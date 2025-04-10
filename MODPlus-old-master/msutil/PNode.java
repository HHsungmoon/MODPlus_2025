package msutil;

import java.util.Comparator;

public class PNode implements Comparable<PNode> {

	double 	mass;
	double 	intensity;
	double 	norm;
	double	prmScore;
	double	b_prm;
	double	y_prm;
	int 	localRank;
	String 	annotation;
	boolean assigned = false;
	
	public PNode(double m, double i, int r, double rn){
		mass = m;
		intensity = i;
		localRank = r;
		norm = rn;
	}
	public double 	getMass() { return mass; }
	public double 	getIntensity() { return intensity; }
	public double 	getNorm() { return norm; }
	public double 	getBPRMScore() { return b_prm; }
	public double 	getYPRMScore() { return y_prm; }
	public boolean 	isAssigned() { return assigned; }
	
	public void assign(boolean a){
		assigned = a;
	}

	public String toString(){
		StringBuffer a = new StringBuffer(annotation);
		a.append(String.format(" %f %f", mass, intensity));
		return a.toString();
	}
	public int compareTo(PNode o) {
		if( this.mass > o.mass ) return 1;
		else if( this.mass == o.mass ) return 0;
		else return -1;
	}	
}

class PNodeIntComparator implements Comparator<PNode>{
	public int compare(PNode x1, PNode x2){
		if( x1.intensity > x2.intensity ) return 1;
		else if( x1.intensity == x2.intensity ) return 0;
		else return -1;
	}	
	public boolean equals(PNode x1, PNode x2){
		return x1 == x2;
	}
}










