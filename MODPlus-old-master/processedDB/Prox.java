package processedDB;

public class Prox implements Comparable<Prox> {
	String accession = "";
	String sequence = "";
	String description = "";
	
	public Prox(){}	
	public Prox(String h, String s){
		accession = h;
		sequence = s;
	}
	
	public String getSequence() { return sequence; }
	public String getAccession() {return accession; }
	
	public void setHeader(String h){ 
		
		int i, cut= 0;
		for( i=0; i<h.length(); i++ ){
			if( h.charAt(i) == '|' || h.charAt(i) == ':' ) cut++;			
			if( cut == 2 || h.charAt(i) == ' ') break;
		}
				
		accession = h.substring(0, i);
		
		if( i != h.length() )
			description = h.substring(i+1);
	}
	
	public void setSequence(String s){ sequence= s; }	
	public String getReverseSequence(){
		StringBuffer seq = new StringBuffer( this.sequence );		
		return seq.reverse().toString();
	}
	
	public int compareTo(Prox p) {
		if( accession.compareTo(p.accession) > 0 ) return 1;
		else if( accession.compareTo(p.accession) == 0 ) return 0;
		else return -1;
	}
	
}



















