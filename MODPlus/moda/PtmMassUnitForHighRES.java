package moda;

import java.text.DecimalFormat;

public class PtmMassUnitForHighRES extends PtmMassUnit {
	
	public PtmMassUnitForHighRES() {
		df = new DecimalFormat("#.###");
	} 
	
	public double getPtmMass(double deltaMass) { return deltaMass; }
	public void setPtmMasses(double[] actptm, int[] intptm) {}
	public String toString(double m){ return df.format(m); }
	
}
