package processedDB;

public class PeptideMatchToProtein implements Comparable<PeptideMatchToProtein> {

	String protein;
	int protIndex;
	int start, end;
	char prevAA, nextAA;
	int ntt;
	
	public PeptideMatchToProtein(String pname, int prot, char pa, int si, int ei, char na, int n){
		protein = pname;
		protIndex = prot; //protein index in fasta
		start 	= si;     //position in protein
		end 	= ei;
		prevAA 	= pa;
		nextAA 	= na;
		ntt 	= n;
	}

	public String getProtName() { return protein; }
	public int getProtIndex() { return protIndex; }
	public int getStartSite() { return start; }

	public String getWrappingAA() { return prevAA+""+nextAA; }
	
	public String toString(){
		return String.format("%s[%d:%c.%d~%d.%c]", protein, protIndex, prevAA, start, end, nextAA);
	}
	
	public int compareTo(PeptideMatchToProtein x) {
		if( ntt < x.ntt ) return 1;
		else if( ntt > x.ntt ) return -1;		

		if( start < x.start ) return -1;
		else if( start > x.start ) return 1;
		
		if( protIndex < x.protIndex ) return -1;
		else if( protIndex > x.protIndex ) return 1;
		
		return 0;
	}
	

}
