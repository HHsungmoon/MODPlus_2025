package processedDB;

public class MODPeptide implements Comparable<MODPeptide> {
	int pStart, pLeft, pRight, pEnd;

	public MODPeptide(int s, int e){
		pStart = s;
		pEnd   = e;
	}
	
	public int getStart()	{ return pStart; }
	public int getEnd()		{ return pEnd; }
	public int getLeft()	{ return pLeft; }
	public int getRight()	{ return pRight; }
	public void setConservedRegion(int a, int b){
		pLeft = a;
		pRight = b;
	}
	public String getPeptide( ProtDatabase trie ){
		return trie.getPeptide(pStart, pEnd);
	}
	
	public boolean extend( MODPeptide xp ){ //one-mod	
		if( this.pStart == xp.pStart && xp.pEnd == this.pEnd ) {
			this.pLeft= Math.max( this.pLeft, xp.pLeft );
			this.pRight= Math.min( this.pRight, xp.pRight );
			return true;
		}
		return false;
	}
	
	public int compareTo( MODPeptide x ) { 
		if( pStart > x.pStart  ) return 1;
		else if( pStart < x.pStart ) return -1;
		else {
			if( pEnd > x.pEnd ) return 1;
			else if( pEnd < x.pEnd ) return -1;
			else return 0;
		}
	}	
}







