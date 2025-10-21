import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DateFormat;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;

import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Vector;
import java.util.concurrent.Callable;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import moda.ThreadPoolManager;
import modi.MatchedTag;
import modi.Peptide;
import modi.TagChain;
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
						int dot = Constants.SPECTRUM_LOCAL_PATH.lastIndexOf('.');
						Constants.SPECTRUM_FILE_NAME = (dot > 0) ? Constants.SPECTRUM_LOCAL_PATH.substring(0, dot) : Constants.SPECTRUM_LOCAL_PATH;
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
		final int numSlots = ThreadPoolManager.numSlots; // 최대 36개 가능

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

		// dataset
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

		// database
		Element database = search.getChild("database");
		if (database != null) {
			Constants.PROTEIN_DB_LOCAL_PATH = database.getAttributeValue("local_path");
			if (Constants.PROTEIN_DB_LOCAL_PATH.equals("")) {
				System.out.println(message[4]);
				return 4;
			}
		}
		System.out.println("Input database : " + Constants.PROTEIN_DB_LOCAL_PATH);

		Element output = search.getChild("output");
		if (output != null) {
			String path = output.getAttributeValue("local_path");
			if (path != null && !path.trim().isEmpty()) {
				Constants.OUTPUT_LOCAL_PATH = path.trim();
				new File(Constants.OUTPUT_LOCAL_PATH).mkdirs(); // 디렉터리 보장
				System.out.println("Output directory : " + Constants.OUTPUT_LOCAL_PATH);
			}
		}

		System.out.println("Output file path : " + Constants.OUTPUT_LOCAL_PATH);

		// enzyme (deprecated) / combined_enzyme
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

		// cysteine alkylation
		Constants.alkylatedToCys = 0; // DEPRECATED
		Element cys_alkylated = search.getChild("cys_alkylated");
		if (cys_alkylated != null) {
			Constants.alkylationMethod = cys_alkylated.getAttributeValue("name");
			Constants.alkylatedToCys = Double.valueOf(cys_alkylated.getAttributeValue("massdiff"));
			AminoAcid.modifiedAminoAcidMass('C', Constants.alkylatedToCys);
			MSMass.modifiedAminoAcidMass('C', Constants.alkylatedToCys);
		}

		// instrument resolution
		Element instrument_resolution = search.getChild("instrument_resolution");
		if (instrument_resolution != null) {
			Constants.MSResolution = ("high".compareToIgnoreCase(instrument_resolution.getAttributeValue("ms")) == 0) ? 1 : 0;
			if (Constants.MSResolution == 1) System.out.println("High resolution MS!!");
			Constants.MSMSResolution = ("high".compareToIgnoreCase(instrument_resolution.getAttributeValue("msms")) == 0) ? 1 : 0;
			if (Constants.MSMSResolution == 1) System.out.println("High resolution MS/MS!!");
		}

		// ===========================================
		// 2) parameters → DEFAULT_* 값 먼저 세팅  // [CHANGED]
		// ===========================================
		Element parameters = search.getChild("parameters");
		if (parameters != null) {
			Element param;

			// enzyme_constraint
			param = parameters.getChild("enzyme_constraint");
			if (param != null) {
				Constants.missCleavages = Integer.valueOf(param.getAttributeValue("max_miss_cleavages"));
				Constants.numberOfEnzymaticTermini = Integer.valueOf(param.getAttributeValue("min_number_termini"));
				if (Constants.numberOfEnzymaticTermini > 2)
					Constants.numberOfEnzymaticTermini = 2;
			}

			// isotope_error → DEFAULT_maxNoOfC13 / rangeForIsotopeIncrement
			param = parameters.getChild("isotope_error");
			if (param != null) {
				if (param.getAttributeValue("min_C13_number") != null)
					Constants.minNoOfC13 = Integer.valueOf(param.getAttributeValue("min_C13_number"));

				if (param.getAttributeValue("max_C13_number") != null) {
					Constants.DEFAULT_maxNoOfC13 = Integer.valueOf(param.getAttributeValue("max_C13_number"));
					if (Constants.DEFAULT_maxNoOfC13 == 0 && param.getAttributeValue("increment_per_dalton") != null)
						Constants.rangeForIsotopeIncrement = Integer.valueOf(param.getAttributeValue("increment_per_dalton"));
				}
			}

			// peptide_mass_tol → PPM or Dalton
			param = parameters.getChild("peptide_mass_tol");
			if (param != null) {
				if ("ppm".equalsIgnoreCase(param.getAttributeValue("unit"))) {
					// ppm: per-scan에서만 환산 적용. DEFAULT는 참고값(0)으로 두어 오용 탐지
					Constants.PPMTolerance = Double.valueOf(param.getAttributeValue("value"));
					Constants.DEFAULT_precursorTolerance = 0.0;  // [CHANGED]
					Constants.DEFAULT_precursorAccuracy  = 0.0;  // [CHANGED]
				} else {
					double tolVal = Double.valueOf(param.getAttributeValue("value"));
					Constants.PPMTolerance = 0.0;                // [CHANGED] da 모드 명시
					Constants.DEFAULT_precursorTolerance = tolVal;
					Constants.DEFAULT_precursorAccuracy  = tolVal;
				}
			}

			// fragment_ion_tol
			param = parameters.getChild("fragment_ion_tol");
			if (param != null) {
				Constants.fragmentTolerance = Double.valueOf(param.getAttributeValue("value"));
				// gapTolerance/gapAccuracy의 per-run 기본은 그대로 둔다(0.6/1.6). per-scan 보정은 스캔 단계에서.
			}

			// modified_mass_range
			param = parameters.getChild("modified_mass_range");
			if (param != null) {
				Constants.minModifiedMass = Double.valueOf(param.getAttributeValue("min_value"));
				Constants.maxModifiedMass = Double.valueOf(param.getAttributeValue("max_value"));
			}
		}

		// DEFAULT_nonModifiedDelta는 런 파라미터 기준으로 설정  // [CHANGED]
		Constants.DEFAULT_nonModifiedDelta = Constants.massToleranceForDenovo;

		// ===========================================
		// 3) 슬롯 배열 생성 + DEFAULT_*로 일괄 리셋   // [CHANGED]
		// ===========================================
		Constants.initSlotArrays(numSlots);

		// ==== protocol(출력용 부가 정보) ====
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

		// ==== modifications / PTMDB 구성 ====
		Element modifications = search.getChild("modifications");

		PTMDB tempFixed = new PTMDB();
		PTMDB tempVar = new PTMDB();

		if (modifications != null) {
			double[] fixedAA = new double[26];
			Element fixed = modifications.getChild("fixed");
			if (fixed != null) {
				if (tempFixed.setFixedModificatinos(fixed, fixedAA) == 0) {
					System.out.println(message[2]);
					return 2;
				}
			}
			if (tempFixed.size() > 0)
				System.out.println("Fixed modifications : " + tempFixed.size() + " selected");

			Element variable = modifications.getChild("variable");
			if (variable != null) {
				Constants.PTM_FILE_NAME = variable.getAttributeValue("local_path");
				boolean canBeModifiedOnFixedAA = variable.getAttributeValue("canBeModifiedOnFixedAA").equals("1");
				Constants.canBeModifiedOnFixedAA = canBeModifiedOnFixedAA;
				if (Constants.PTM_FILE_NAME != null) {
					tempVar.setVariableModificatinos(Constants.PTM_FILE_NAME, fixedAA, canBeModifiedOnFixedAA);
				}
				tempVar.setVariableModificatinos(variable, fixedAA, canBeModifiedOnFixedAA);

				if (canBeModifiedOnFixedAA) {
					for (PTM p : tempFixed) {
						tempVar.add(new PTM(tempVar.size(), "De-" + p.getName(), "",
							-p.getMassDifference(), 0, p.getResidue(), p.getPTMPosition(), (p.getAbbAA() == 'C') ? 1 : 0));
					}
				}
				if (variable.getAttributeValue("multi_mods") != null && variable.getAttributeValue("multi_mods").equals("0")) {
					Constants.maxPTMPerGap = Constants.maxPTMPerPeptide = 1;
				}
				if (tempVar.size() > 0) {
					System.out.print("Variable modifications : " + tempVar.size() + " selected (");
					tempVar.setPTMDiagnosticIon();
					if (Constants.maxPTMPerPeptide == 1)
						System.out.println("one modification per peptide)");
					else
						System.out.println("multiple modifications per peptide)");
				}
			}
		}

		// 모든 슬롯에 대해 PTMDB 객체를 동일한 내용으로 초기화 및 lookup table 구성
		Constants.fixedModifications = new PTMDB[numSlots];
		Constants.variableModifications = new PTMDB[numSlots];
		for (int i = 0; i < numSlots; i++) {
			Constants.fixedModifications[i] = tempFixed.deepCopy();
			Constants.variableModifications[i] = tempVar.deepCopy();
			Constants.variableModifications[i].constructPTMLookupTable();
		}

		// ==== decoy / multistage / mod_map ====
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

		// 파라미터 보정(기존 로직)
		Constants.adjustParameters();

		System.out.println("---------print parameter----------");
		// Constants.printAllConstantsState(1);
		System.out.println();
		return 0;
	}

	// javac -cp ".:lib/*" -d out $(find . -name "*.java")
	// java -cp "out:lib/*" MODPlus param.xml
	// nohup java -cp "out:lib/*" MODPlus param.xml > 09B_QE3.txt 2>&1 &  -> 백그라운드 실행
	static int modplus_mod_search() throws Exception {
		System.out.println("Thread Number : " + ThreadPoolManager.numSlots);
		System.out.println("Starting MODPlus for modification search!");

		Constants.MAX_TAG_SIZE = 100;
		Constants.minTagLength = 2;
		Constants.minTagLengthPeptideShouldContain = 3;
		Constants.tagChainPruningRate = 0.4;

		//String identifier = Constants.SPECTRUM_LOCAL_PATH;
		//identifier = identifier.substring(0, identifier.lastIndexOf('.'));
		String identifier = Constants.OUTPUT_LOCAL_PATH;
		//identifier = identifier.substring(0, identifier.lastIndexOf('.'));

		// 1. SpectrumContainer 생성
		ScanIterator scaniter = ScanIterator.get(Constants.SPECTRUM_LOCAL_PATH, Constants.SPECTRA_FILE_TYPE);
		if (scaniter == null || scaniter.size() == 0) {
			System.out.println("Failed to read msms spectra file");
			return 1;
		}
		System.out.println(scaniter.size() + " scans");

		// 2. Protein DB 생성
		System.out.print("Reading protein database.....  ");
		StemTagTrie ixPDB = new StemTagTrie(Constants.PROTEIN_DB_LOCAL_PATH);
		if (ixPDB.getSizeOfEntries() == 0) {
			System.out.println("Failed to read protein fasta file");
			return 1;
		}
		System.out.println();
		long startTime = System.currentTimeMillis();

		// 3. 스캔 데이터 미리 수집
		List<ArrayList<MSMScan>> allScanBlocks = new ArrayList<>();
		while (scaniter.hasNext()) {
			allScanBlocks.add(scaniter.getNext());
		}

		// 4. 병렬 처리 실행
		List<Future<String>> futureList = runParallelSearch(allScanBlocks, ixPDB, scaniter.getFileName());

		// 5. 결과 수집 및 출력
		try (PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(identifier+ Constants.SPECTRUM_FILE_NAME +"_P.txt")))) {
			for (int i = 0; i < futureList.size(); i++) {
				try {
					out.print(futureList.get(i).get());
				} catch (ExecutionException e) {
					System.err.println("⚠ Error in block " + (i + 1) + ": " + e.getCause());
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		}

		System.out.println("[MOD-Plus] Elapsed Time : " + (System.currentTimeMillis() - startTime) / 1000 + " Sec");
		return 0;
	}


	private static List<Future<String>> runParallelSearch(List<ArrayList<MSMScan>> scanBlocks,
			StemTagTrie ixPDB, String fileName) throws InterruptedException {
		final int numThreads = ThreadPoolManager.numSlots;
		ThreadPoolManager.enterParallelPhase();
		ExecutorService executor = java.util.concurrent.Executors.newFixedThreadPool(
			numThreads,
			r -> {
				Thread t = new Thread(() -> {
					ThreadPoolManager.bindSlotForWorker();
					r.run();
				});
				t.setName("modplus-worker-" + t.getId());
				t.setDaemon(false);
				return t;
			}
		);

		try {
			List<Future<String>> futures = new ArrayList<>(scanBlocks.size());
			for (int blockIdx = 0; blockIdx < scanBlocks.size(); blockIdx++) {
				final int idx = blockIdx;
				final ArrayList<MSMScan> block = scanBlocks.get(blockIdx);
				futures.add(executor.submit(() -> {
					final int slot = ThreadPoolManager.getSlotIndex();
					Constants.resetSlot(slot); // 슬롯 경계에서 초기화
					try {
						return processBlock(idx+1, scanBlocks.size(), block, ixPDB, fileName);
					} finally {
						Constants.resetSlot(slot);
					}
				}));
			}

			return futures;
		} finally {
			executor.shutdown();
			ThreadPoolManager.exitParallelPhase();
		}
	}

	private static String processBlock(int blockIndex, int totalBlocks, ArrayList<MSMScan> chargedSpectra, StemTagTrie ixPDB, String fileName)
            throws Exception {
		int slotIdx = ThreadPoolManager.getSlotIndex();
		StringBuilder resultSB = new StringBuilder();

		System.out.println("MODPlus | " + blockIndex + "/" + totalBlocks);

		int selected = -1;
		ArrayList<AnsPeptide> candidates = null;

		for (int i = 0; i < chargedSpectra.size(); i++) {
			Spectrum spectrum = chargedSpectra.get(i).getSpectrum();

			if (spectrum.getObservedMW() > Constants.maxPeptideMass)
				continue;

			PGraph graph = spectrum.getPeakGraph();
			double correctedMW = graph.correctMW(dynamicPMCorrection);
			spectrum.setCorrectedParentMW(correctedMW);

			TagPool tPool = SpectrumAnalyzer.buildTagPool(spectrum);
			DPHeap heatedPepts = OneMOD.getHeatedPeptides(ixPDB, graph, tPool,
					(Constants.maxNoOfC13[slotIdx] != 0 || Constants.precursorTolerance[slotIdx] > 0.50001));

			DPHeap tepidPepts = null;
			if (Constants.maxPTMPerPeptide > 1 && (heatedPepts == null || !heatedPepts.isConfident())) {
				tepidPepts = heatedPepts;
				heatedPepts = MultiMOD.getHeatedPeptides(ixPDB, graph, tPool, dynamicPMCorrection);
				//System.out.println("MultiMOD : "+blockIndex);
			}

			if (heatedPepts == null)
				continue;

			HeatedDB bitDB = getHeatedDB(ixPDB, heatedPepts, tepidPepts);
			TagTrie bitTrie = bitDB.getPartialDB(ixPDB);
			ArrayList<AnsPeptide> tp = dynamicMODeye(bitTrie, graph, tPool);

			if (!tp.isEmpty()) {
				if (candidates == null || candidates.get(0).compareTo(tp.get(0)) == 1) {
					candidates = tp;
					selected = i;
				}
			}
		}

		if (selected != -1) {
			MSMScan scan = chargedSpectra.get(selected);
			HashMap<String, ArrayList<PeptideMatchToProtein>> seqToProtMap = new HashMap<>();

			resultSB.append(">> ").append(scan.getHeader()).append("\n");
			for (AnsPeptide peptide : candidates) {
				String tpSeq = peptide.getPeptideSequence();
				ArrayList<PeptideMatchToProtein> matchedProteins = seqToProtMap.computeIfAbsent(tpSeq, ixPDB::getMatchProteins);
				resultSB.append(peptide.toMODPlus(scan.getObservedMW(), matchedProteins)).append("\n");
			}
			resultSB.append("\n");
		}

		return resultSB.toString();
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

		// Step 1: MatchedTagPool 생성
		MatchedTagPool matchedList = SpectrumAnalyzer.extendedBuildMatchedTagPool(
				tPool, parentMW, dynamicDB, Constants.protease, Constants.numberOfEnzymaticTermini);

		if (matchedList == null || matchedList.size() == 0) return new ArrayList<>();

		// Step 2: TagChainPool 안전하게 생성
		TagChainPool tcPool = TagChainPool.buildTagChainPoolSafely(matchedList);
		tcPool.discardPoorTagChain();

		// Step 3: PTM 해석 적용
		boolean specAnnotated = false;
		if (!tcPool.isEmpty()) {
			specAnnotated = SpectrumAnalyzer.interpretTagChain(
					Constants.variableModifications[ThreadPoolManager.getSlotIndex()], tcPool, graph);
		}

		// Step 4: 최종 후보 펩타이드 추출
		ArrayList<AnsPeptide> cands = new ArrayList<>();
		if (!tcPool.isEmpty() && specAnnotated) {
			cands = tcPool.getAnswerPeptides(graph);
		}

		return cands;
	}
}
