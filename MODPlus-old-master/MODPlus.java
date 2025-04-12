import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;

import moda.ThreadPoolManager;
import org.jdom.Document;
import org.jdom.Element;
import org.jdom.JDOMException;
import org.jdom.input.SAXBuilder;

import moda.DPHeap;
import moda.DPPeptide;
import moda.MultiMOD;
import moda.OneMOD;
import modi.AminoAcid;
import modi.AnsPeptide;
import modi.Constants;
import modi.MatchedTagPool;
import modi.PTM;
import modi.PTMDB;
import modi.Spectrum;
import modi.SpectrumAnalyzer;
import modi.TagChainPool;
import modi.TagPool;
import msutil.IsobaricTag;
import msutil.MSMass;
import msutil.PGraph;
import msutil.ProtCutter;
import processedDB.HeatedDB;
import processedDB.PeptideMatchToProtein;
import processedDB.ProtDatabase;
import processedDB.StemTagTrie;
import processedDB.TagTrie;
import scaniter.MSMScan;
import scaniter.ScanIterator;


public class MODPlus {
	static boolean dynamicPMCorrection = false, multiBlind= true;
	
	static int numHeatedPeptides = 50;
	
	private static String[] message = {
		"[Error] Cannot read any MS/MS scan from input dataset.\r\n" +
		"[Error] Check consistency between input file and its format.",
		
		"[Error] Cannot read any protein from input database.\r\n" +
		"[Error] Check input fasta format.",
		
		"[Error] One fixed modification per amino acid can be allowed.\r\n"+
		"[Error] Check specfied fixed modifications.",
		
		"[Error] Unsupported character set in your search parameter",
		
		"[Error] Required field is empty.\r\n"+
		"[Error] Required fields : MS/MS Data, Database", 
		
		"[Error] Wrong usage.\r\n"+
		"[Error] Re-confirm it.",
		
		"[Error] Not defined"
	};
	
	public static void main( String[] args ) throws Exception {
		Constants.engine = "modplus";
		Constants.engineVersion= "hyu";
		
		System.out.println("************************************************************************************");
		System.out.println("Modplus (version "+Constants.engineVersion+") - Identification of post-translational modifications");
		System.out.println("Release Date: Apr 27, 2015");
		System.out.println("************************************************************************************");
		System.out.println();
		
		run(args[0]);
	}
	
	public static void run( String arg  ) throws Exception {
		try {			
			if( set_parameter( arg ) != 0 ) return;			
		} catch (Exception e) {
			e.printStackTrace();
			return;
		}		
		
		try {
			File analPath = new File( Constants.SPECTRUM_LOCAL_PATH );		
			if( analPath.isDirectory() ) {
				String type = Constants.SPECTRA_FILE_TYPE.toString().toLowerCase();
				for( File file : analPath.listFiles() ){
						
					if( file.getName().endsWith(type) ){
						Constants.SPECTRUM_LOCAL_PATH = file.getPath();
						System.out.println("Input datasest : "+Constants.SPECTRUM_LOCAL_PATH);
						modplus_mod_search();
					}
				}	
				System.out.println("End of process");
			}
			else modplus_mod_search();
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	protected static int set_parameter(String Prixparam) throws Exception {

		System.out.println("Reading parameter.....");

		// 최대 슬롯 수를 ThreadPoolManager에서 받아옵니다.
		int numSlots = ThreadPoolManager.numSlots; // 최대 36개 가능

		// 전역 변수들을 슬롯 배열로 초기화합니다.
		Constants.precursorTolerance = new double[numSlots];
		Constants.precursorAccuracy = new double[numSlots];
		Constants.gapTolerance = new double[numSlots];
		Constants.gapAccuracy = new double[numSlots];
		Constants.nonModifiedDelta = new double[numSlots];
		Constants.maxNoOfC13 = new int[numSlots];

		for (int i = 0; i < numSlots; i++) {
			Constants.precursorTolerance[i] = 0.5;
			Constants.precursorAccuracy[i] = 0.5;
			Constants.gapTolerance[i] = 0.6;
			Constants.gapAccuracy[i] = 1.6;
			Constants.nonModifiedDelta[i] = Constants.massToleranceForDenovo;
			Constants.maxNoOfC13[i] = 0;
		}

		// XML 문서를 읽어옵니다.
		Document doc;
		try {
			doc = new SAXBuilder().build(Prixparam);
		} catch (JDOMException e) {
			System.out.println(message[3]);
			return 1;
		} catch (IOException e) {
			System.out.println(message[5]);
			return 5;
		}

		Element search = doc.getRootElement();
		Constants.runDate = DateFormat.getDateInstance().format(new Date());
		if (search.getAttributeValue("user") != null) {
			Constants.runUser = search.getAttributeValue("user");
		}
		if (search.getAttributeValue("title") != null) {
			Constants.runTitle = search.getAttributeValue("title");
		} else {
			Constants.runTitle = String.valueOf(System.currentTimeMillis());
		}

		Element dataset = search.getChild("dataset");
		if (dataset != null) {
			Constants.SPECTRUM_LOCAL_PATH = dataset.getAttributeValue("local_path");
			if (Constants.SPECTRUM_LOCAL_PATH.equals("")) {
				System.out.println(message[4]);
				return 4;
			}

			String type = dataset.getAttributeValue("format");
			if (type.compareToIgnoreCase("mgf") == 0)
				Constants.SPECTRA_FILE_TYPE = Constants.spectra_format.MGF;
			else if (type.compareToIgnoreCase("pkl") == 0)
				Constants.SPECTRA_FILE_TYPE = Constants.spectra_format.PKL;
			else if (type.compareToIgnoreCase("ms2") == 0)
				Constants.SPECTRA_FILE_TYPE = Constants.spectra_format.MS2;
			else if (type.compareToIgnoreCase("dta") == 0)
				Constants.SPECTRA_FILE_TYPE = Constants.spectra_format.DTA;
			else if (type.compareToIgnoreCase("mzxml") == 0)
				Constants.SPECTRA_FILE_TYPE = Constants.spectra_format.MZXML;
			else if (type.compareToIgnoreCase("zip") == 0)
				Constants.SPECTRA_FILE_TYPE = Constants.spectra_format.ZIPDTA;

			Constants.INSTRUMENT_NAME = dataset.getAttributeValue("instrument");
			if (Constants.INSTRUMENT_NAME.equals("QTOF"))
				Constants.INSTRUMENT_TYPE = Constants.msms_type.QTOF;
			else
				Constants.INSTRUMENT_TYPE = Constants.msms_type.TRAP;
		}
		System.out.print("Input datasest : " + Constants.SPECTRUM_LOCAL_PATH);
		System.out.println(" (" + Constants.SPECTRA_FILE_TYPE + " type)");

		Element database = search.getChild("database");
		if (database != null) {
			Constants.PROTEIN_DB_LOCAL_PATH = database.getAttributeValue("local_path");
			if (Constants.PROTEIN_DB_LOCAL_PATH.equals("")) {
				System.out.println(message[4]);
				return 4;
			}
		}
		System.out.println("Input database : " + Constants.PROTEIN_DB_LOCAL_PATH);

		Element enzyme = search.getChild("enzyme"); // DEPRECATED
		if (enzyme != null) {
			String enzymeName = enzyme.getAttributeValue("name");
			String cut = enzyme.getAttributeValue("cut");
			String sence = enzyme.getAttributeValue("sence");
			Constants.protease = new ProtCutter(enzymeName, cut, sence);
		}

		Element com_enzyme = search.getChild("combined_enzyme");
		if (com_enzyme != null) {
			String enzymeName = com_enzyme.getAttributeValue("name");
			String nn = com_enzyme.getAttributeValue("nterm_cleave");
			String cc = com_enzyme.getAttributeValue("cterm_cleave");
			Constants.protease = new ProtCutter(enzymeName, nn, cc, true);
		}

		Constants.alkylatedToCys = 0; // DEPRECATED
		Element cys_alkylated = search.getChild("cys_alkylated");
		if (cys_alkylated != null) {
			Constants.alkylationMethod = cys_alkylated.getAttributeValue("name");
			Constants.alkylatedToCys = Double.valueOf(cys_alkylated.getAttributeValue("massdiff"));
			AminoAcid.modifiedAminoAcidMass('C', Constants.alkylatedToCys);
			MSMass.modifiedAminoAcidMass('C', Constants.alkylatedToCys);
		}

		Element instrument_resolution = search.getChild("instrument_resolution");
		if (instrument_resolution != null) {
			Constants.MSResolution = ("high".compareToIgnoreCase(instrument_resolution.getAttributeValue("ms")) == 0) ? 1 : 0;
			if (Constants.MSResolution == 1)
				System.out.println("High resolution MS!!");
			Constants.MSMSResolution = ("high".compareToIgnoreCase(instrument_resolution.getAttributeValue("msms")) == 0) ? 1 : 0;
			if (Constants.MSMSResolution == 1)
				System.out.println("High resolution MS/MS!!");
		}

		Element parameters = search.getChild("parameters");
		if (parameters != null) {
			Element param;
			param = parameters.getChild("enzyme_constraint");
			if (param != null) {
				Constants.missCleavages = Integer.valueOf(param.getAttributeValue("max_miss_cleavages"));
				Constants.numberOfEnzymaticTermini = Integer.valueOf(param.getAttributeValue("min_number_termini"));
				if (Constants.numberOfEnzymaticTermini > 2)
					Constants.numberOfEnzymaticTermini = 2;
			}

			param = parameters.getChild("isotope_error");
			if (param != null) {
				if (param.getAttributeValue("min_C13_number") != null)
					Constants.minNoOfC13 = Integer.valueOf(param.getAttributeValue("min_C13_number"));

				if (param.getAttributeValue("max_C13_number") != null)
				{
					for(int i=0; i<numSlots; i++){
						Constants.maxNoOfC13[i] = Integer.valueOf(param.getAttributeValue("max_C13_number"));
						if (Constants.maxNoOfC13[i] == 0 && param.getAttributeValue("increment_per_dalton") != null)
							Constants.rangeForIsotopeIncrement = Integer.valueOf(param.getAttributeValue("increment_per_dalton"));
					}
				}
			}

			param = parameters.getChild("peptide_mass_tol");
			if (param != null) {
				if (param.getAttributeValue("unit").compareToIgnoreCase("ppm") == 0) {
					Constants.PPMTolerance = Double.valueOf(param.getAttributeValue("value"));
				} else {
					// 여기서는 precursorTolerance와 precursorAccuracy를 읽은 값을 slot으로 초기화된 배열에 저장 후 사용
					double tolVal = Double.valueOf(param.getAttributeValue("value"));
					for (int i = 0; i < numSlots; i++) {
						Constants.precursorTolerance[i] = tolVal;
						Constants.precursorAccuracy[i] = tolVal;
					}
				}
			}

			param = parameters.getChild("fragment_ion_tol");
			if (param != null) {
				Constants.fragmentTolerance = Double.valueOf(param.getAttributeValue("value"));
			}
			param = parameters.getChild("modified_mass_range");
			if (param != null) {
				Constants.minModifiedMass = Double.valueOf(param.getAttributeValue("min_value"));
				Constants.maxModifiedMass = Double.valueOf(param.getAttributeValue("max_value"));
			}
		}

		Element protocol = search.getChild("protocol");
		if (protocol != null) {
			System.out.print("Protocol Description: ");
			Element isobaric = protocol.getChild("isobaric_labeling");
			if (isobaric != null) {
				if (isobaric.getAttributeValue("reagent") != null) {
					Constants.isobaricTag = isobaric.getAttributeValue("reagent");
					Constants.reporterMassOfIsobaricTag = IsobaricTag.getReporterMasses(isobaric.getAttributeValue("reagent"));
					if (!Constants.isobaricTag.equals(""))
						System.out.print(Constants.isobaricTag + " Labelled" + ((Constants.reporterMassOfIsobaricTag == null) ? " (NOT Supported)" : " (Supported)"));
				}
			}
			Element modEnrich = protocol.getChild("modification_enrichment");
			if (modEnrich != null) {
				if (modEnrich.getAttributeValue("mod") != null) {
					Constants.enrichedModification = modEnrich.getAttributeValue("mod");
					if (!Constants.enrichedModification.equals(""))
						System.out.print(" & " + Constants.enrichedModification + " Enriched" +
								(("Acetyl".compareToIgnoreCase(Constants.enrichedModification) == 0 || "Phospho".compareToIgnoreCase(Constants.enrichedModification) == 0) ? " (Supported)" : " (NOT Supported)"));
				}
			}
			System.out.println();
		}

		Element modifications = search.getChild("modifications");

		// PTMDB 슬롯 초기화
		Constants.variableModifications = new PTMDB[numSlots];
		Constants.fixedModifications = new PTMDB[numSlots];
		for (int i = 0; i < numSlots; i++) {
			// 새로운 PTMDB 객체 생성 (필요시 clone이나 deep copy 적용)
			Constants.variableModifications[i] = new PTMDB();
			Constants.fixedModifications[i] = new PTMDB();
		}

		if (modifications != null) {
			double[] fixedAA = new double[26];

			Element fixed = modifications.getChild("fixed");
			if (fixed != null) {
				if (Constants.fixedModifications[0].setFixedModificatinos(fixed, fixedAA) == 0) {
					System.out.println(message[2]);
					return 2;
				}
			}
			if (Constants.fixedModifications[0].size() > 0)
				System.out.println("Fixed modifications : " + Constants.fixedModifications[0].size() + " selected");

			Element variable = modifications.getChild("variable");
			if (variable != null) {
				Constants.PTM_FILE_NAME = variable.getAttributeValue("local_path");
				boolean canBeModifiedOnFixedAA = variable.getAttributeValue("canBeModifiedOnFixedAA").equals("1");
				Constants.canBeModifiedOnFixedAA = canBeModifiedOnFixedAA;
				if (Constants.PTM_FILE_NAME != null) {
					Constants.variableModifications[0].setVariableModificatinos(Constants.PTM_FILE_NAME, fixedAA, canBeModifiedOnFixedAA);
				}
				Constants.variableModifications[0].setVariableModificatinos(variable, fixedAA, canBeModifiedOnFixedAA);

				if (canBeModifiedOnFixedAA) {
					for (PTM p : Constants.fixedModifications[0]) {
						Constants.variableModifications[0].add(
								new PTM(Constants.variableModifications[0].size(), "De-" + p.getName(), "",
										-p.getMassDifference(), 0, p.getResidue(), p.getPTMPosition(), (p.getAbbAA() == 'C') ? 1 : 0));
					}
				}
				if (variable.getAttributeValue("multi_mods") != null && variable.getAttributeValue("multi_mods").equals("0")) {
					Constants.maxPTMPerGap = Constants.maxPTMPerPeptide = 1;
				}
			}
			if (Constants.variableModifications[0].size() > 0) {
				System.out.print("Variable modifications : " + Constants.variableModifications[0].size() + " selected (");
				Constants.variableModifications[0].setPTMDiagnosticIon();
				if (Constants.maxPTMPerPeptide == 1)
					System.out.println("one modification per peptide)");
				else
					System.out.println("multiple modifications per peptide)");
			}
		}
		Constants.variableModifications[0].constructPTMLookupTable();

		Element decoy_search = search.getChild("decoy_search");
		if (decoy_search != null) {
			if (decoy_search.getAttributeValue("checked") != null) {
				if ("1".compareTo(decoy_search.getAttributeValue("checked")) == 0) {
					Constants.targetDecoy = 1;
					System.out.println("Decoy search checked");
				}
			}
		}

		Element multistages_search = search.getChild("multistages_search");
		if (multistages_search != null) {
			if (multistages_search.getAttributeValue("checked") != null) {
				if ("1".compareTo(multistages_search.getAttributeValue("checked")) == 0) {
					Constants.multiStagesSearch = 1;
					Constants.firstSearchProgram = multistages_search.getAttributeValue("program");
					System.out.println("MultiStages Search checked " + Constants.firstSearchProgram);
				}
			}
		}

		Element mod_map = search.getChild("mod_map");
		if (mod_map != null) {
			if (mod_map.getAttributeValue("checked") != null) {
				if ("1".compareTo(mod_map.getAttributeValue("checked")) == 0) {
					Constants.runMODmap = 1;
					System.out.println("MODMap checked");
				}
			}
		}

		Constants.adjustParameters();

		System.out.println("---------print parameter----------");
		// Constants.printAllConstantsState(1);
		System.out.println();
		return 0;
	}

	public static void printConstantsState(int index) {
		System.out.println("--------------------" + index + "-------------------");

		System.out.println("precursorTolerance        = " + Constants.precursorTolerance);
		System.out.println("precursorAccuracy         = " + Constants.precursorAccuracy);
		System.out.println("gapTolerance              = " + Constants.gapTolerance);
		System.out.println("gapAccuracy               = " + Constants.gapAccuracy);
		System.out.println("nonModifiedDelta          = " + Constants.nonModifiedDelta);
		System.out.println("maxNoOfC13                = " + Constants.maxNoOfC13);

		System.out.println("minModifiedMass           = " + Constants.minModifiedMass);
		System.out.println("maxModifiedMass           = " + Constants.maxModifiedMass);
		System.out.println("fragmentTolerance         = " + Constants.fragmentTolerance);
		System.out.println("PPMTolerance              = " + Constants.PPMTolerance);
		System.out.println("rangeForIsotopeIncrement  = " + Constants.rangeForIsotopeIncrement);

		System.out.println("canBeModifiedOnFixedAA     = " + Constants.canBeModifiedOnFixedAA);
		System.out.println("maxPTMPerGap               = " + Constants.maxPTMPerGap);
		System.out.println("maxPTMPerPeptide           = " + Constants.maxPTMPerPeptide);
		System.out.println("===================================================");
	}

	static int modplus_mod_search() throws Exception{
		System.out.println("Starting MODPlus for modification search!");
		
		Constants.MAX_TAG_SIZE = 100;
		Constants.minTagLength = 2;

		Constants.minTagLengthPeptideShouldContain = 3;
		Constants.tagChainPruningRate = 0.4;
		
		String identifier = Constants.SPECTRUM_LOCAL_PATH;
		identifier = identifier.substring(0, identifier.lastIndexOf('.'));
		
		// SpectrumContainer generation
		ScanIterator scaniter = ScanIterator.get( Constants.SPECTRUM_LOCAL_PATH, Constants.SPECTRA_FILE_TYPE );
		if( scaniter == null || scaniter.size() == 0 ){
			System.out.println("Failed to read msms spectra file" );
			return 1;
		}
		System.out.println( scaniter.size() + " scans"  );		
		
		// Protein DB generation
		System.out.print( "Reading protein database.....  " );		
		StemTagTrie ixPDB= new StemTagTrie( Constants.PROTEIN_DB_LOCAL_PATH );
		if( ixPDB.getSizeOfEntries() == 0 ){
			System.out.println( "Failed to read protein fasta file" );
			return 1;
		}
		System.out.println();
				
		long startTime= System.currentTimeMillis();		
		PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(identifier+".modplus.txt")));

		boolean[] considerIsotopeErr = new boolean[ThreadPoolManager.numSlots];
		for(int i=0; i<ThreadPoolManager.numSlots; i++) {
			considerIsotopeErr[i] =
					(Constants.maxNoOfC13[i] != 0 || Constants.precursorTolerance[i] > 0.50001)
							? true : false;
		}

		int index = 1;
		int iterSize= scaniter.size();
		while( scaniter.hasNext() ){
			
			ArrayList<MSMScan> chargedSpectra = scaniter.getNext();
			System.out.println("MODPlus | " + (index++) + "/" + iterSize);

			//Constants.printAllConstantsState(slot); // 슬롯별 Constants 값 확인

			int selected = -1;
			ArrayList<AnsPeptide> candidates = null;

			for(int i=0; i<chargedSpectra.size(); i++){
				Spectrum spectrum = chargedSpectra.get(i).getSpectrum();
				System.out.println(" └─ Spectrum Charge " + spectrum.getCharge() + " | ObservedMW: " + spectrum.getObservedMW());

				if (spectrum.getObservedMW() > Constants.maxPeptideMass) {
					System.out.println("    → Skipped: Over maxPeptideMass");
					continue;
				}

				PGraph graph = spectrum.getPeakGraph();
				double correctedMW = graph.correctMW(dynamicPMCorrection);
				spectrum.setCorrectedParentMW(correctedMW);
				System.out.println("    → Corrected MW: " + correctedMW);

				TagPool tPool = SpectrumAnalyzer.buildTagPool( spectrum );

				DPHeap heatedPepts = OneMOD.getHeatedPeptides( ixPDB, graph, tPool, considerIsotopeErr[ThreadPoolManager.getSlotIndex()] );
				DPHeap tepidPepts  = null;
				if( Constants.maxPTMPerPeptide > 1 )
					if( heatedPepts == null || !heatedPepts.isConfident() ) {
						tepidPepts = heatedPepts;
						heatedPepts = MultiMOD.getHeatedPeptides( ixPDB, graph, tPool, dynamicPMCorrection );
					}

				if( heatedPepts == null ) continue;

				HeatedDB bitDB = getHeatedDB(ixPDB, heatedPepts, tepidPepts);
				TagTrie bitTrie = bitDB.getPartialDB(ixPDB);
				System.out.println("bitTrie의 전체 태그 개수 = " + bitTrie.getTotalTagCount());
				System.out.println("bitTrie에서 TAG 'A-R-G'의 등장 수 = " + bitTrie.getTagHitCount('A', 'R', 'G'));


				ArrayList<AnsPeptide> tp = dynamicMODeye(bitTrie, graph, tPool);
				System.out.println("    → dynamicMODeye results: " + tp.size() + " candidate(s)");

				if( tp.size() > 0 ) { //decision Of CS
					if( candidates == null || candidates.get(0).compareTo( tp.get(0) ) == 1 ) {
						candidates = tp;
						selected = i;
						System.out.println("    → Selected candidate updated");
					}
				}
			}
			
			if( selected != -1 ) {
				MSMScan scan = chargedSpectra.get(selected);
				
				HashMap<String, ArrayList<PeptideMatchToProtein>> seqToProtMap = new HashMap<String, ArrayList<PeptideMatchToProtein>>();
				out.println( ">>"+scaniter.getFileName()+"\t"+scan.getHeader() );

				for( int k=0; k<candidates.size(); k++ ){
					String tpSeq = candidates.get(k).getPeptideSequence();
					ArrayList<PeptideMatchToProtein> matchedProteins = seqToProtMap.get(tpSeq);
										
					if( matchedProteins == null ){
						matchedProteins = ixPDB.getMatchProteins(tpSeq);
						seqToProtMap.put(tpSeq, matchedProteins);
					}				
					out.println( candidates.get(k).toMODPlus(scan.getObservedMW(), matchedProteins) );	
				}
				out.println();	
			}

			System.out.println("----------------------");
			//Constants.printAllConstantsState(1);

		}
		out.close();
		
		System.out.println("[MOD-Plus] Elapsed Time : " + (System.currentTimeMillis()-startTime)/1000 + " Sec" );
		return 0;
	}
	
	private static HeatedDB getHeatedDB( StemTagTrie stemDB, DPHeap candidates, DPHeap tepids ) {
		
		HeatedDB matchedBits= new HeatedDB();
		int count = 0;
		for( DPPeptide dp : candidates ){
			if( dp.getScore() < 1 ) break;		
			String modapept = dp.getPeptide();		
			int pro_start = dp.getProtein();		
			ProtDatabase proDB = stemDB.get(dp.getStem());
			matchedBits.add( proDB.getProteinIdentity(pro_start), dp.getStem(), pro_start, pro_start+modapept.length() );				
			if( ++count == numHeatedPeptides ) break;
		}

		count = 0;
		if( tepids != null ){
			for( DPPeptide dp : tepids ){
				if( dp.getScore() < 1 ) break;		
				String modapept = dp.getPeptide();		
				int pro_start = dp.getProtein();		
				ProtDatabase proDB = stemDB.get(dp.getStem());
				matchedBits.add( proDB.getProteinIdentity(pro_start), dp.getStem(), pro_start, pro_start+modapept.length() );			
				if( ++count == 10 ) break;
			}
		}
		return matchedBits;
	}

	private static ArrayList<AnsPeptide> dynamicMODeye(TagTrie dynamicDB, PGraph graph, TagPool tPool) throws Exception {
		double parentMW = graph.getCorrectedMW();

		System.out.printf("    → Corrected MW: %.9f%n", parentMW);

		MatchedTagPool matchedList = SpectrumAnalyzer.extendedBuildMatchedTagPool(
				tPool, parentMW, dynamicDB, Constants.protease, Constants.numberOfEnzymaticTermini);

		System.out.println("matchedList size: " + matchedList.size());

		TagChainPool tcPool = new TagChainPool();
		tcPool.putAll(SpectrumAnalyzer.buildTagChain(matchedList));
		System.out.println("TC Pool Size Before discard: " + tcPool.size());

		tcPool.discardPoorTagChain();
		System.out.println("TC Pool Size After discard: " + tcPool.size());

		boolean specAnnotated = false;
		if (tcPool.size() > 0) {
			System.out.println("Calling interpretTagChain()...");
			specAnnotated = SpectrumAnalyzer.interpretTagChain( Constants.variableModifications[ThreadPoolManager.getSlotIndex()], tcPool, graph );
			System.out.println("interpretTagChain() result: " + specAnnotated);
		}

		ArrayList<AnsPeptide> cands = new ArrayList<>();
		if (tcPool.size() > 0 && specAnnotated) {
			cands = tcPool.getAnswerPeptides(graph);
		}

		System.out.println("TC Pool Size After : " + tcPool.size());
		System.out.println("PTM Combination 개수: " + cands.size());

		// 상세 TagChain 내부 내용 출력 (상위 몇 개만 보기)


		return cands;
	}

	
}
