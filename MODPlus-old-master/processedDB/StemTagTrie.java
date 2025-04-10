package processedDB;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Collections;

public class StemTagTrie extends ArrayList<TagTrie> {

	static private int 	capacityPerStem = 50 * 1024 * 1024;
	
	private int 		sizeOfEntries;
	private int 		sizeOfResidues;	
	
	public StemTagTrie( String fileName ) throws Exception{
		construct(fileName);
	}

	public int getSizeOfEntries(){ return sizeOfEntries; }
	
	
	private void construct(String fileName){
		
		sizeOfEntries = sizeOfResidues = 0;
		
		int stemAA = 0;
		ProxDB stem= new ProxDB();
		
		try {									
			BufferedReader in = new BufferedReader(new FileReader(fileName));
			String s= "<";
			StringBuffer buffer;	
			
			while( !s.startsWith(">") ) {
				s = in.readLine();
				if( s == null ) break;
			}
		
			while( s != null ) {
				Prox protein= new Prox();				
				protein.setHeader( s.substring(1) );			
				
				buffer = new StringBuffer();
				while( (s = in.readLine()) != null ){
					if( s.startsWith(">") ) break;
					for( int aa=0; aa<s.length(); aa++ ){
						if( Character.isLetter( s.charAt(aa) ) )
							buffer.append( Character.toUpperCase(s.charAt(aa)) );
					}
				}
				if( buffer.length() < 3 ) continue; // check sequence
				
				protein.setSequence( buffer.toString() );
				stemAA += buffer.length();
				stem.add( protein );
				
				if( capacityPerStem < stemAA ){
					stem.setSizeOfEntries(stem.size());
					stem.setSizeOfResidues(stemAA);
					TagTrie tt = new TagTrie(stem);
					tt.setStemNo( this.size() );
					this.add( tt );
					sizeOfEntries += stem.size();
					sizeOfResidues += stemAA;
					stemAA = 0;
					stem= new ProxDB();
				}
					
			}
			if( 0 < stemAA ){
				stem.setSizeOfEntries(stem.size());
				stem.setSizeOfResidues(stemAA);
				TagTrie tt = new TagTrie(stem);
				tt.setStemNo( this.size() );
				this.add( tt );
				sizeOfEntries += stem.size();
				sizeOfResidues += stemAA;
				stemAA = 0;
				stem= null;
			}
			in.close();						
			System.out.println( sizeOfEntries+" proteins / "  + sizeOfResidues+" residues" + " ("+this.size()+")" );
			
		}		
		catch (FileNotFoundException e) {
			System.out.println( "Cannot find the protein file, "+fileName );
			e.printStackTrace();
		} 
		catch (Exception e) {
			System.out.println( "Cannot read the protein file" );
			e.printStackTrace();
		}
	}
	
	public ArrayList<PeptideMatchToProtein> getMatchProteins(String peptide){
		ArrayList<PeptideMatchToProtein> protMatch = new ArrayList<PeptideMatchToProtein>();
		int prior = 0;
		for(int i=0; i<this.size(); i++){
			protMatch.addAll( this.get(i).getMatchProteins(peptide, prior) );
			prior += this.get(i).getSizeOfEntries();
		}
		Collections.sort(protMatch);
		return protMatch;
	}
	
}
