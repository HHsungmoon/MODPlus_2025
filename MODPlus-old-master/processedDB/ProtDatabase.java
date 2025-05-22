package processedDB;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import java.nio.charset.StandardCharsets;
import msutil.ProtCutter;
import modi.Constants;
import modi.Peptide;

public class ProtDatabase {
	static byte delimeter = '-';
	protected byte[] 	sequence;
	protected int 		at_end_of_line;
	
	protected String[] 	proteins;
	protected int[] 	protPos;
	private int 		sizeOfEntries;
	private int 		sizeOfResidues;	
	
	public ProtDatabase( String fileName ) throws Exception {
		construct( fileName );
	}
	
	public ProtDatabase( ProxDB proxDB ){	
		construct(proxDB);
	}
	
	private void construct(ProxDB proxDB){
		
		sizeOfResidues = proxDB.getSizeOfResidues();
		sizeOfEntries = proxDB.size();
		
		proteins = new String[sizeOfEntries];		
		protPos  = new int[sizeOfEntries];

		int res = 0;
		sequence = new byte[sizeOfResidues+sizeOfEntries+1];

		sequence[res++] = delimeter;
		
		int id=0;
		for( Prox pro : proxDB ){	
			String s = pro.getSequence();		
			for(int k=0; k<s.length(); k++){
				sequence[res++] = (byte)s.charAt(k);
			}
			
			proteins[id] = pro.getAccession();
			protPos[id++] = res;
			
			sequence[res++] = delimeter;
		}
		at_end_of_line = res - 1;
	}

	public int getSizeOfEntries(){ return sizeOfEntries; }
	
	public String getProteinIdentity(int position) { return proteins[searchProtein(position)]; }

	/*
	public String getPeptide(int s, int e){ //inclusive, exclusive like substring
		StringBuffer x= new StringBuffer();
		for(int i=s; i<e; i++){
			x.append( (char)sequence[i] );			
		}	
		return x.toString();
	}
	*/

	public String getPeptide(int s, int e) {
		return new String(sequence, s, e - s, StandardCharsets.ISO_8859_1);
	}
	
	public String getSequenceAroundPeptide(int s, int e, int around){ //inclusive, exclusive like substring
		StringBuffer x= new StringBuffer();
		int i = s-around;
		if( i < 1 ) i = 1;
		if( sequence[i-1] != delimeter ) {
		//	x.append( 'X' );
			if( sequence[i-2] != delimeter ) x.append( 'X' );
			else x.append( (char)sequence[i-1] );
		}
		for( ; i<s; i++){
			x.append( (char)sequence[i] );
			if( sequence[i] == delimeter ) x = new StringBuffer();
		}		
		for(i=s; i<e+around; i++){
			if( sequence[i] == delimeter ) break;
			x.append( (char)sequence[i] );
		}
		if( sequence[i] != delimeter ) x.append( 'X' );
		return x.toString();
	}
	
	public int	getProteinIndex(int position) { 
		int index = -1;
		if( ( index = searchProtein(position) ) < sizeOfEntries )
			return index;
		return -1; 
	}
	public int getRealPositionInProtein(int position) { 
		int protein = searchProtein( position ); 		
		if( protein == 0 ) return position;
		else return ( position - protPos[protein-1] );	
	}
	
	protected boolean isPotentialProteinNterm(int pos) {
		if( sequence[pos-1] == 'M' && sequence[pos-2] == delimeter ) return true;
		else return false;
	}
	public int getNTTOfPeptide( int start, int end, ProtCutter cutter ) { //in moda dynamic programming
		if( !cutter.isSpecified() ) return 0;
		
		int NTT = 0;	
		if( cutter.isCleavage(sequence[start-1], sequence[start]) || isPotentialProteinNterm(start) ) NTT++;		
		if( cutter.isCleavage(sequence[end-1], sequence[end]) ) NTT++;	
		return NTT;
	}

	public char[] getTeminalResiduesOfPeptide( Peptide pept ){
		char[] nc = new char[2];
		nc[0] = (char)sequence[pept.getStartCoordinate()-1];
		nc[1] = (char)sequence[pept.getEndCoordinate()+1];
		
		return nc;
	}
	
	protected int searchProtein( int pos ){
		int index;	
		if( pos <= protPos[0] )
			index= 0;
		else if( pos > protPos[protPos.length-1] )
			index= protPos.length-1;
		else{
			int M, L= 0, R= protPos.length-1;
			while( R - L > 1 )
			{
				M= ( L + R ) /2;

				if( pos <= protPos[M] )
					R= M;
				else
					L= M;
			}
			index= R;
		}	
		return index;
	}
	
	private int getMemorySizeForSequence(String fileName)
	{
		sizeOfEntries = sizeOfResidues = 0;
		try {							
			BufferedReader in = new BufferedReader(new FileReader(fileName));
			String buf;
			while( (buf=in.readLine()) != null ){
				if( buf.startsWith(">") ) sizeOfEntries++;
				else sizeOfResidues += buf.length();
			}
			in.close();		
		}		
		catch (FileNotFoundException e) {
			System.out.println( "Cannot find the protein file, "+fileName );
			e.printStackTrace();
			System.exit(1);
		} 
		catch (Exception e) {
			System.out.println( "Cannot read the protein file" );
			e.printStackTrace();
			System.exit(1);
		}		
	
		return (sizeOfEntries+sizeOfResidues+1);
	}
	
	private void construct( String fileName ) throws IOException  {		
		int memSize = getMemorySizeForSequence( fileName );
		if( Constants.targetDecoy == 1 ) {
			memSize *= 2;
			sizeOfEntries *= 2;
		}
		
		sequence = new byte[memSize];
		proteins = new String[sizeOfEntries];		
		protPos  = new int[sizeOfEntries];
		sizeOfEntries = sizeOfResidues = 0;
		
		int ent = 0, res = 0;
		sequence[res++] = delimeter;
		
		BufferedReader in = new BufferedReader(new FileReader(fileName));
		String s= "<";
		StringBuffer buffer;	
		
		while( s.startsWith(">") == false ) {
			s = in.readLine();
			if( s == null ) break;
		}
	
		while( s != null ) {
			proteins[ent] = getAccession(s);//s.substring(1);
			
			buffer = new StringBuffer();
			while( (s = in.readLine()) != null ){
				if( s.startsWith(">") ) break;
				for( int aa=0; aa<s.length(); aa++ ){
					if( Character.isLetter( s.charAt(aa) ) )
						buffer.append( Character.toUpperCase(s.charAt(aa)) );
				}
			}
			if( buffer.length() < 3 ) continue; // check sequence

			for( int aa=0; aa<buffer.length(); aa++ ){
				sequence[res++] = (byte)buffer.charAt(aa);
			}
			protPos[ent++] = res;
			sequence[res++] = delimeter;
			
			sizeOfEntries++;
			sizeOfResidues += buffer.length();

		}
		in.close();			
		System.out.println( sizeOfEntries+" proteins / "  + sizeOfResidues+" residues" );
		
		if( Constants.targetDecoy == 1 ) {
			System.out.print("Decoy search checked.....  ");
			Constants.DECOY_DB_NAME = String.valueOf( sizeOfEntries )+" reversed sequences";
			int border = res;
			StringBuffer forward = new StringBuffer();
			for(int p=1; p<border; p++){
				if( sequence[p] == delimeter ){
					proteins[ent] = Constants.DECOY_LABEL + proteins[ent-sizeOfEntries];
					for( int aa=forward.length()-1; aa>-1; aa-- ){
						sequence[res++] = (byte)forward.charAt(aa);
					}
					protPos[ent++] = res;
					sequence[res++] = delimeter;
					forward = new StringBuffer();
				}
				else forward.append( (char)sequence[p] );
			}
			System.out.println( "Generated " + sizeOfEntries+" decoy proteins");
		}
		at_end_of_line = res - 1;
		while( ent < proteins.length  ){
			proteins[ent] = "";
			protPos[ent++] = res++;
		}
	}
	
	private String getAccession(String h){ 
		int i, cut= 0;
		for( i=1; i<h.length(); i++ ){
			if( h.charAt(i) == '|' || h.charAt(i) == ':' ) cut++;			
			if( cut == 2 || h.charAt(i) == ' ') break;
		}				
		return h.substring(1, i);
	}
	
}




