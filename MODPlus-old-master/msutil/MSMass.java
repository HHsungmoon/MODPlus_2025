package msutil;

import modi.Constants;

public class MSMass {
	private static double minAAMass = 57.02146;
	private static double maxAAMass = 186.07931;
	private static double aaMass[]={
		71.03711, 0, 103.00919, 115.02694, 129.04259, 
		147.06841, 57.02146, 137.05891, 113.08406, 0, 
		128.09496, 113.08406, 131.04049, 114.04293, 0, 
		97.05276, 128.05858, 156.10111, 87.03203, 101.04768,
		0, 99.06841, 186.07931, 0, 163.06333, 0};
	
	public static double getAAMass(char aa){
		return aaMass[aa-'A'];
	}
	public static double getAAMass(byte aa){
		return aaMass[aa-'A'];
	}
	
	public static double getPepMass(String pept){
		double mass=0;
		for(int i=0; i<pept.length(); i++){			
			if( aaMass[pept.charAt(i)-'A'] == 0 )
				return -1;
			mass += aaMass[pept.charAt(i)-'A'];
		}
		return mass;
	}
	public static double getMinAAMass(){
		return minAAMass;
	}
	public static void modifiedAminoAcidMass(char AA, double fixedModification){
		aaMass[AA-'A'] += fixedModification;
		if( aaMass[AA-'A'] > maxAAMass ) maxAAMass = aaMass[AA-'A'];
		if( aaMass[AA-'A'] < minAAMass ) minAAMass = aaMass[AA-'A'];
	}
	
	public static boolean isIndistinguishableAA(char a, char b){
		if( Math.abs(aaMass[a-'A']-aaMass[b-'A']) <= Constants.massToleranceForDenovo ) return true;
		return false;
	}
}














