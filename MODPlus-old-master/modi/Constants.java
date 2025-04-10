package modi;

import java.text.DecimalFormat;

import scaniter.MSMScan;

import msutil.MSMass;
import msutil.ProtCutter;

public class Constants {
	
	public static String engine;
	public static String engineVersion;

	public static String runDate;
	public static String runUser= "anonymous";
	public static String runTitle;
	
	public static String 			SPECTRUM_LOCAL_PATH;
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
	public static int			maxNoOfC13 = 0;
	public static int			rangeForIsotopeIncrement = 0;
	
	public static double		alkylatedToCys = 0;
	public static String		alkylationMethod;
	
	public static double		precursorAccuracy = 0.5;
	public static double		precursorTolerance = 0.5;
	public static double		PPMTolerance = 0;
	public static double		fragmentTolerance = 0.6;
	public static double		gapTolerance = 0.6;	
	public static double		gapAccuracy = 1.6;	
	public static double		minNormIntensity = 0.00;
	
	public static PTMDB 		variableModifications;
	public static PTMDB 		fixedModifications;
	public static double		minModifiedMass = -precursorTolerance;
	public static double		maxModifiedMass = precursorTolerance;
	public static boolean		canBeModifiedOnFixedAA = false;
	public static boolean		isInModifiedRange( double v ){
		if( minModifiedMass-gapTolerance < v && v < maxModifiedMass+gapTolerance ) return true;
		else if( Math.abs(v) <= gapTolerance ) return true;
		else return false;
	}
	
	public static int			MSResolution 	= 0; // if 1, high (FT, OrbiTrap)
	public static int			MSMSResolution 	= 0; // if 1, high (FT, OrbiTrap)
	
	//for De novo sequencing
	public static double		massToleranceForDenovo = 0.3;
	public static int 			MAX_TAG_SIZE = 50;
	public static double		selectionWindowSize   = 70;
	public static int			minNumOfPeaksInWindow = 4;
	public static int			maxNumOfPeaksInWindow = 8;
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
	
	// for Peptide DB
	public static final int		proteinIDModeSeqLength	= 3;
	public static final String 	SOURCE_PROTEIN_FILE_NAME = "sourceProtein.mprot";	
	
	public static final String	UNIMOD_FILE_NAME = "unimod.xml";
	
	// for mother mass correction for LTQ/LCQ
	public static final double	MINIMUM_PRECURSOR_MASS_ERROR = -1.5;
	public static final double	MAXIMIM_PRECURSOR_MASS_ERROR = 1.5;

	// if true, write unidrawing only tag chains whose all gaps are annotated
	public static final boolean	writeAnnotatedTagChainOnly = false;
	
	public static final int		MINIMUM_SHARED_PEAK_COUNT = 2; 
	
	// for offset
	public static final int newLineCharSize = new String("\r\n").getBytes().length;
	
	public static double		nonModifiedDelta = massToleranceForDenovo;
	public static int[] 		maxPTMOccurrence = {1, 1, 2, 2, 3, 3, 3, 2};// = new int[7];
	
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
			if( fixedModifications.size() > 0 ){
				for( PTM p : fixedModifications ){
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
		if( Math.abs(v1-v2) <= precursorTolerance ) return true;
		else return false;
	}	
	
	public static final double[] rNorm= {6,
		2.928968, 1.928968, 1.428968, 1.095635, 0.845635,
		0.645635, 0.478968, 0.336111, 0.211111, 0.100000};
	
	public static double[] coEfft= {0.3159, -34.6288, 1.3209, -8.7609, 0., - 5.0206};		
	public static double getMODScore( double a, double b, double c, double d, double e){
		return coEfft[0]*a + coEfft[1]*b + coEfft[2]*c + coEfft[3]*d + coEfft[4]*e + coEfft[5]; 		
	}
	
	public static String	getString(double value){
		return new DecimalFormat("#.###").format(value).toString();
	}
	
	public static double	MASS_CAL_STD_THRESHOLD = 0.1;
	public static double	PTM_ADD_PENALTY = 0.2;
	
	// Can use Wysocki paper results?
	public static double	getMissingPenaltyWeight(PeakProperty property){
		if(property == PeakProperty.Y_ION)
			return 0.4;	
		else if(property == PeakProperty.Y_MINUS_NH3_ION)
			return 0;
		else if(property == PeakProperty.Y_MINUS_H2O_ION)
			return 0;
		else if(property == PeakProperty.B_ION)
			return 0;
		else if(property == PeakProperty.A_ION)
			return 0;
		else
			return 0;
	}
	public static double	getNotExplainedPenaltyWeight(){
		return 0.15;
	}
	public static final	double	ANALYSIS_VERSION = 0.8;
	
	public static boolean isWithinTolerance(double calc, double obsv, double tol){

		if( minNoOfC13 ==0 && maxNoOfC13 == 0 ) {
			if( Math.abs(calc-obsv) > tol ) return false;
		}
		else {
			double tempError = obsv - calc;		
			int isoerr = round( tempError / IsotopeSpace );		
			if( isoerr < minNoOfC13 || maxNoOfC13 < isoerr ) return false;		
			if(	Math.abs( tempError - isoerr*IsotopeSpace ) > precursorAccuracy ) return false;
		}
		return true;
	}
	public static boolean isWithinAccuracy(double err){		
		if( gapAccuracy > 0.5 ) return true;
		int isoerr = round( err / IsotopeSpace );		
		if(	Math.abs( err - isoerr*IsotopeSpace ) > gapAccuracy ) return false;
		return true;
	}
	
	public static double PPMtoDalton(double mass, double ppm){	
		return mass/1000000*ppm;
	}

	public static int round(double a){
		if( a > 0 ) return (int)(a + 0.5);
		else return (int)(a - 0.5);
	}
}
