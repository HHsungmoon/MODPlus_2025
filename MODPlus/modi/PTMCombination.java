package modi;

public class PTMCombination {
	String 		ptmComb = "";
	double[] 	ptms;
	PTM[] 		ptmList;
	
	public PTMCombination(int len){
		ptms = new double[len];
		ptmList = new PTM[len];
	}
	
	public int hashCode(){ 
		return ptmComb.hashCode();
	}
	
	public boolean equals(Object o){		
		if(!(o instanceof PTMCombination))
			   return false;		
		PTMCombination p= (PTMCombination)o;
		if( this.ptmComb.equals(p.ptmComb) )
			return true;
		return false;
	}
}
