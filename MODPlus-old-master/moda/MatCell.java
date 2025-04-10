package moda;

public class MatCell {
	double 		mass;
	double 		score;
	int 		isInsideTag; // before 0, inside 1, after 2 
	int 		isAAJump = -1;
	double 		delta;
	int 		nominalDelta;
	MatCell		parent = null;
	
	MatCell(int a){	
		isInsideTag = a; 			
	}

	void refresh(){
		isAAJump= -1;
		score	= 0;
	}
	
	void setMass(double baseMass, double deltaMass, int deltaInt){
		mass		 = baseMass + deltaMass;
		delta		 = deltaMass;
		nominalDelta = deltaInt;
	}
	
	void correctMass(double unit){
		mass += unit;
		delta+= unit;	
		nominalDelta++;
	}
	
	void setParent(MatCell mc){
		parent = mc;
	}
	
	boolean canObliqueJumpFrom( MatCell a ){		
		if( a.isAAJump < 1 )
			return false;
		if( a.isInsideTag > 0 && this.isInsideTag < 2 )
			return true; 
		return false;
	}
}
