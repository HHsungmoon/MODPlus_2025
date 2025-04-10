package moda;

import modi.Constants;
import msutil.MSMass;

public class MODaConst {
	
	static int maxTagPoolSize = 50; 
	
	static double 	baseScore = -10000.;
	static double 	isotopeUnit = ( Constants.MSMSResolution==1 )? Constants.IsotopeSpace : 1.0;
	static PtmMassUnit ptmUnit = ( Constants.MSMSResolution==1 )? new PtmMassUnitForHighRES() : new PtmMassUnit();
	
	static double 	maxIsotopeError =  -1*isotopeUnit;
	static int	 	maxIntIsotopeError =  -1;
	static int 		isotopePointsToBeCorrected = ( Constants.INSTRUMENT_TYPE == Constants.msms_type.QTOF )? 2 : 3;
	
	static double 	TERMINALMOD = Constants.NTERM_FIX_MOD + Constants.CTERM_FIX_MOD;
	static double 	minimumDistance = MSMass.getMinAAMass() - Constants.fragmentTolerance;
	 
	static double 	minCutOfProb = 0.9;
	static int		minConsecutiveAALength = 5;	
}
