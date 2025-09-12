package modi;

import java.text.DecimalFormat;

import moda.ThreadPoolManager;

import msutil.MSMass;
import msutil.ProtCutter;

public class Constants {
	
	public static String engine;
	public static String engineVersion;

	public static String runDate;
	public static String runUser= "anonymous";
	public static String runTitle;
	
	public static String 			SPECTRUM_LOCAL_PATH;
	public static String 			OUTPUT_LOCAL_PATH;
	public static String 			SPECTRUM_FILE_NAME;
	public static String 			INSTRUMENT_NAME = "TRAP";
	
	public static msms_type	 		INSTRUMENT_TYPE = msms_type.TRAP; //TOF(0), LOW_TRAP(1), HIGH_TRAP(2)
	public static spectra_format 	SPECTRA_FILE_TYPE = spectra_format.MGF;	
	
	public enum spectra_format {
		PKL,		// read spectrums in SPECTRUM_FILE_NAME
		DTA,		// read all dta file from SPECTRUM_FILE_NAME(compressed file)
		MGF,		// read spectrums in SPECTRUM_FILE_NAME
		MS2,
		MZXML,		// read spectrums in SPECTRUM_FILE_NAME
		ZIPDTA, 
	}
	
	public enum msms_type {
		QTOF,
		TRAP,
	}

	public enum instrument_resolution {
		HIGH,
		LOW
	}

	public enum experimental_protocol {
		iTRAQ,
	}
	
	public static String 		PROTEIN_DB_LOCAL_PATH;
	public static String 		PROTEIN_DB_NAME;
	public static String 		DECOY_LABEL="dec_";
	public static String 		DECOY_DB_LOCAL_PATH;
	public static String 		DECOY_DB_NAME;
	
	public static int			targetDecoy=0;
	public static int			runMODmap = 0;
	
	public static int 			multiStagesSearch = 0;
	public static String 		firstSearchProgram = "";

	public static final double	UNIT_MASS = 1.;
	
	public static final double	Electron = 0.000549;
	public static final double	Hydrogen = 1.007825035;
	public static final double	Oxygen = 15.99491463;
	public static final double	Nitrogen = 14.003074;
	public static final double	Proton = Hydrogen-Electron;
	public static final double	HO = Hydrogen + Oxygen;
	public static final double	H2O = Hydrogen*2 + Oxygen;	
	public static final double	NH3 = Hydrogen*3 + Nitrogen;		
	public static final double	IsotopeSpace = 1.00235;
	
	public static double		NTERM_FIX_MOD = 0;
	public static double		CTERM_FIX_MOD = 0;
	public static final double	B_ION_OFFSET = Proton;
	public static final double	Y_ION_OFFSET = H2O + Proton;
	public static final double	A_ION_OFFSET = Oxygen + 12.;
	public static final double  IMM_OFFSET = -A_ION_OFFSET + Proton;

	public static double		minPeptideMass = 300.;
	public static double		maxPeptideMass = 5000.;//
	
	public static ProtCutter 	protease = ProtCutter.getCutter("Trypsin");
	public static int			numberOfEnzymaticTermini = 2;
	public static int			missCleavages = 2;
	
	public static int			minNoOfC13 = 0;
	public static int			rangeForIsotopeIncrement = 0;
	
	public static double		alkylatedToCys = 0;
	public static String		alkylationMethod;

	// 초기화를 위해
	public static double DEFAULT_precursorTolerance  = 0.5;
	public static double DEFAULT_precursorAccuracy   = 0.5;
	public static double DEFAULT_gapTolerance        = 0.6;
	public static double DEFAULT_gapAccuracy         = 1.6;
	public static double DEFAULT_nonModifiedDelta    = 0.3; // 너희 기본에 맞춰 조정
	public static int    DEFAULT_maxNoOfC13          = 0;


	// 병렬처리 slot 방식
	public static double[]		precursorTolerance; // 0.5
	public static double[] 		precursorAccuracy; // 0.5
	public static double[] 		gapTolerance; // 0.6
	public static double[] 		gapAccuracy;  //  1.6
	public static double[] 		nonModifiedDelta; // massToleranceForDenovo
	public static int[]    		maxNoOfC13; // 0

	public static void resetSlot(int idx) {
		precursorTolerance[idx] = DEFAULT_precursorTolerance;
		precursorAccuracy[idx]  = DEFAULT_precursorAccuracy;
		gapTolerance[idx]       = DEFAULT_gapTolerance;
		gapAccuracy[idx]        = DEFAULT_gapAccuracy;
		nonModifiedDelta[idx]   = DEFAULT_nonModifiedDelta;
		maxNoOfC13[idx]         = DEFAULT_maxNoOfC13;
	}

	public static void initSlotArrays(int numSlots) {
		precursorTolerance = new double[numSlots];
		precursorAccuracy  = new double[numSlots];
		gapTolerance       = new double[numSlots];
		gapAccuracy        = new double[numSlots];
		nonModifiedDelta   = new double[numSlots];
		maxNoOfC13         = new int[numSlots];
		for (int i = 0; i < numSlots; i++) resetSlot(i);
	}

	public static double		PPMTolerance = 0;
	public static double		fragmentTolerance = 0.6;
	public static double		minNormIntensity = 0.00;

	// 근데 이거 두개도 문제 생길 수 있음.
	public static PTMDB[] 		variableModifications;
	public static PTMDB[] 		fixedModifications;
	public static double		minModifiedMass = -0.5;
	public static double		maxModifiedMass = 0.5;
	public static boolean		canBeModifiedOnFixedAA = false;
	public static boolean		isInModifiedRange( double v ){
		if( minModifiedMass-gapTolerance[ThreadPoolManager.getSlotIndex()] < v && v < maxModifiedMass+gapTolerance[ThreadPoolManager.getSlotIndex()] ) return true;
		else if( Math.abs(v) <= gapTolerance[ThreadPoolManager.getSlotIndex()] ) return true;
		else return false;
	}
	
	public static int			MSResolution 	= 0; // if 1, high (FT, OrbiTrap)
	public static int			MSMSResolution 	= 0; // if 1, high (FT, OrbiTrap)
	
	//for De novo sequencing
	public static double		massToleranceForDenovo = 0.3;
	public static int 			MAX_TAG_SIZE = 50;
	public static double		selectionWindowSize   = 70;
	public static int			minNumOfPeaksInWindow = 4;
	public static int			minTagLength = 3;
	public static int			minTagLengthPeptideShouldContain = 3;
	public static boolean		Leu_indistinguishable_Ile = true;
	public static boolean		Lys_indistinguishable_Qln = true;

	public static double		tagChainPruningRate = 0.5;
	public static int			maxTagPerPept     	= 12;
	public static int			maxTagChainPerPept  = 30;	
	public static int 			maxInterpretationPerGap	= 10;	
	public static int			maxPTMPerGap		= 2;
	public static int 			maxPTMPerPeptide	= 4;
	
	public static int getMaxPTMOccurrence( int seqLength ){		
		if( seqLength > 10 ) return 1;
		return maxPTMPerGap;
	}
	
	public static String		PTM_FILE_NAME = "PTMDB.xml";
	
	public static String		isobaricTag = "";
	public static double[]		reporterMassOfIsobaricTag = null;
	
	public static String		enrichedModification = "";
	
	public static void	adjustParameters(){
		if( INSTRUMENT_TYPE == msms_type.QTOF ) { // TOF
			massToleranceForDenovo = ( MSMSResolution == 0 )? 0.2 : 0.2;	
			minNumOfPeaksInWindow = 4;
			rNorm[0]= 6;
		}
		else {
			massToleranceForDenovo = ( MSMSResolution == 0 )? 0.3 : 0.03;	
			minNumOfPeaksInWindow = 4;
			rNorm[0]= 6;
		}
		if( massToleranceForDenovo > fragmentTolerance/2 ) massToleranceForDenovo = fragmentTolerance/2;
		if( fragmentTolerance < 0.1 ) MSMSResolution = 1;
		Constants.Lys_indistinguishable_Qln = MSMass.isIndistinguishableAA('K', 'Q');
		Constants.Leu_indistinguishable_Ile = MSMass.isIndistinguishableAA('L', 'I');
		
		if( canBeModifiedOnFixedAA ){			
			double fixedOff = -20;
			if( fixedModifications[ThreadPoolManager.getSlotIndex()].size() > 0 ){
				for( PTM p : fixedModifications[ThreadPoolManager.getSlotIndex()] ){
					fixedOff -= p.getMassDifference();
				}
				if( fixedOff < minModifiedMass ) minModifiedMass = fixedOff;
			}
		}
	}
	
	public static boolean	fEqual(double v1, double v2){
		if( Math.abs(v1-v2) <= fragmentTolerance ) return true;
		else return false;
	}
	
	public static boolean	pEqual(double v1, double v2){
		if( Math.abs(v1-v2) <= precursorTolerance[ThreadPoolManager.getSlotIndex()] ) return true;
		else return false;
	}	
	
	public static final double[] rNorm= {6,
		2.928968, 1.928968, 1.428968, 1.095635, 0.845635,
		0.645635, 0.478968, 0.336111, 0.211111, 0.100000};

	public static String	getString(double value){
		return new DecimalFormat("#.###").format(value).toString();
	}

	public static boolean isWithinTolerance(double calc, double obsv, double tol){

		if( minNoOfC13 ==0 && maxNoOfC13[ThreadPoolManager.getSlotIndex()] == 0 ) {
			if( Math.abs(calc-obsv) > tol ) return false;
		}
		else {
			double tempError = obsv - calc;		
			int isoerr = round( tempError / IsotopeSpace );		
			if( isoerr < minNoOfC13 || maxNoOfC13[ThreadPoolManager.getSlotIndex()] < isoerr ) return false;
			if(	Math.abs( tempError - isoerr*IsotopeSpace ) > precursorAccuracy[ThreadPoolManager.getSlotIndex()] ) return false;
		}
		return true;
	}
	public static boolean isWithinAccuracy(double err){		
		if( gapAccuracy[ThreadPoolManager.getSlotIndex()] > 0.5 ) return true;
		int isoerr = round( err / IsotopeSpace );		
		if(	Math.abs( err - isoerr*IsotopeSpace ) > gapAccuracy[ThreadPoolManager.getSlotIndex()] ) return false;
		return true;
	}
	
	public static double PPMtoDalton(double mass, double ppm){	
		return mass/1000000*ppm;
	}

	public static int round(double a){
		if( a > 0 ) return (int)(a + 0.5);
		else return (int)(a - 0.5);
	}

	public static void printAllConstantsState(int slot) {
		System.out.println("======= Constants 상태 (slot = " + slot + ") =======");

		// Slot-based
		System.out.println("precursorTolerance        = " + precursorTolerance);
		System.out.println("precursorAccuracy         = " + precursorAccuracy);
		System.out.println("gapTolerance              = " + gapTolerance);
		System.out.println("gapAccuracy               = " + gapAccuracy);
		System.out.println("nonModifiedDelta          = " + nonModifiedDelta);
		System.out.println("maxNoOfC13                = " + maxNoOfC13);

		// Shared (non-slot)
		System.out.println("minModifiedMass           = " + minModifiedMass);
		System.out.println("maxModifiedMass           = " + maxModifiedMass);
		System.out.println("fragmentTolerance         = " + fragmentTolerance);
		System.out.println("PPMTolerance              = " + PPMTolerance);
		System.out.println("rangeForIsotopeIncrement  = " + rangeForIsotopeIncrement);
		System.out.println("massToleranceForDenovo    = " + massToleranceForDenovo);
		System.out.println("minNumOfPeaksInWindow     = " + minNumOfPeaksInWindow);
		System.out.println("minTagLength              = " + minTagLength);
		System.out.println("minTagLengthPeptideShouldContain = " + minTagLengthPeptideShouldContain);

		System.out.println("maxPTMPerGap              = " + maxPTMPerGap);
		System.out.println("maxPTMPerPeptide          = " + maxPTMPerPeptide);
		System.out.println("canBeModifiedOnFixedAA    = " + canBeModifiedOnFixedAA);
		System.out.println("alkylatedToCys            = " + alkylatedToCys);
		System.out.println("alkylationMethod          = " + alkylationMethod);

		System.out.println("MSResolution              = " + MSResolution);
		System.out.println("MSMSResolution            = " + MSMSResolution);
		System.out.println("maxPeptideMass            = " + maxPeptideMass);

		System.out.println("engine                    = " + engine);
		System.out.println("engineVersion             = " + engineVersion);
		System.out.println("runDate                   = " + runDate);
		System.out.println("runUser                   = " + runUser);
		System.out.println("runTitle                  = " + runTitle);

		System.out.println("SPECTRUM_LOCAL_PATH       = " + SPECTRUM_LOCAL_PATH);
		System.out.println("PROTEIN_DB_LOCAL_PATH     = " + PROTEIN_DB_LOCAL_PATH);
		System.out.println("INSTRUMENT_NAME           = " + INSTRUMENT_NAME);
		System.out.println("INSTRUMENT_TYPE           = " + INSTRUMENT_TYPE);
		System.out.println("SPECTRA_FILE_TYPE         = " + SPECTRA_FILE_TYPE);

		System.out.println("targetDecoy               = " + targetDecoy);
		System.out.println("runMODmap                 = " + runMODmap);
		System.out.println("multiStagesSearch         = " + multiStagesSearch);
		System.out.println("firstSearchProgram        = " + firstSearchProgram);

		System.out.println("PTM_FILE_NAME             = " + PTM_FILE_NAME);
		System.out.println("isobaricTag               = " + isobaricTag);
		System.out.println("enrichedModification      = " + enrichedModification);

		System.out.println("Leu_indistinguishable_Ile = " + Leu_indistinguishable_Ile);
		System.out.println("Lys_indistinguishable_Qln = " + Lys_indistinguishable_Qln);

		System.out.println("===================================================");
	}
}
