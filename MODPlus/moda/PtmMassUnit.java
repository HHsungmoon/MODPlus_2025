package moda;

import java.text.DecimalFormat;

import modi.Constants;

public class PtmMassUnit { //for integer ptm
	
	protected DecimalFormat df;

	public PtmMassUnit() {
		df = new DecimalFormat("#");
	}
	
	public double getPtmMass(double deltaMass) { return Constants.round(deltaMass); }
	public void setPtmMasses(double[] actptm, int[] intptm ) { 
		for(int i=0; i<actptm.length; i++) {
			actptm[i] = intptm[i];
		}
	}
	public String toString( double m ){ return df.format(m); }
	
}
