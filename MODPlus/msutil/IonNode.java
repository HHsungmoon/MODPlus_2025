package msutil;

import modi.Constants;

public class IonNode implements Comparable<IonNode> {
	int 	index;
	int 	type; // nterm=1, cterm=2
	int 	charge;
	double 	mz; // theoretical
	int 	priority;
	double 	observed;
	double 	intensity;
	
	public IonNode(int i, int t, int cs, double m, int p){
		index = i;
		type = t;
		mz = m;
		charge =cs;
		priority = p;
	}
	
	public int getIndex(){ return index; }
	public int getCharge(){ return charge; }
	
	public int compareTo(IonNode o) {
		if( this.priority > o.priority ) return 1;
		else if( this.priority == o.priority ){
			if( this.mz > o.mz ) return 1;
			else if( this.mz < o.mz ) return -1;
			else return 0;
		}
		else return -1;
	}	
}
