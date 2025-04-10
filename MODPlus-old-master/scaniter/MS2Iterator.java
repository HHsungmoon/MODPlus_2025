package scaniter;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.StringTokenizer;

public class MS2Iterator extends ScanIterator {
	
	private String reversedStartLine=""; 
		
	public MS2Iterator( String fileName ) throws IOException {		
		super(fileName);
		
		BufferedReader in = new BufferedReader( new FileReader(fileName) );		
		String buf;
		int size = 0;			
		while( (buf = in.readLine()) != null ){
			if( buf.startsWith("S") )// begining of new spectrum
				size++;
		}	
		in.close();
		sizeOfScans = size;
		scanIndex = 1;
		
		fin = new BufferedReader( new FileReader(fileName) );
		while( (buf = fin.readLine()) != null ){
			 if( buf.startsWith("S") ) {
				 reversedStartLine = buf;
				 break;
			 }
		}
	}
	
	public ArrayList<MSMScan> getNext() throws IOException {
		
		ArrayList<MSMScan> scanlist = new ArrayList<MSMScan>();
		MSMScan curScan = null;
		
		String buf = reversedStartLine;
		if( buf.startsWith("S") ){// begining of new spectrum				
			StringTokenizer token = new StringTokenizer(buf);
			if( token.countTokens() != 4 ) {
				 while( (buf = fin.readLine()) != null ){
					 if( buf.startsWith("S") ) {
						 reversedStartLine = buf;
						 scanIndex++;
						 return scanlist;
					 }
				 }
			}
			
			token.nextToken();//S
			String startScan = token.nextToken();
			String endScan   = token.nextToken();
			String title= baseName+"."+startScan+"."+endScan;
			int scanNo = Integer.parseInt(startScan);
			
			double pmz = Double.parseDouble( token.nextToken() );
			ArrayList<Integer> cslist = new ArrayList<Integer>();	

			while( (buf = fin.readLine()) != null ) {	
				
				if( !Character.isLetter(buf.charAt(0)) ) break;
				
				if( buf.startsWith("Z") ){
					token = new StringTokenizer(buf);
					if( token.countTokens() != 3 ) break;
					token.nextToken();//Z
					cslist.add(Integer.parseInt(token.nextToken()));
				}
			}
			
			ArrayList<RawPeak> rawPL = new ArrayList<RawPeak>();				
			do {
				if( buf.startsWith("S") ) {
					reversedStartLine = buf;
					break;
				}
				token = new StringTokenizer(buf);
				if( token.countTokens() > 1 )
					rawPL.add( new RawPeak(Double.parseDouble(token.nextToken()), Double.parseDouble(token.nextToken())) );
			} while( (buf = fin.readLine()) != null );
			Collections.sort(rawPL);

			if( pmz > 0 ){
				if( cslist.size() != 0 ){
					for( Integer cs : cslist ){
						curScan = new MSMScan(title+"."+cs, scanIndex, scanNo, pmz, cs);
						if( curScan.setSpectrum(rawPL) ) scanlist.add(curScan);
					}
				}
				else{
					for(int cs=MIN_ASSUMED_CHARGE; cs<=MAX_ASSUMED_CHARGE; cs++){						
						curScan = new MSMScan(title+"."+cs, scanIndex, scanNo, pmz, cs);
						if( curScan.setSpectrum(rawPL) ) scanlist.add(curScan);
					}
				}
			}
			scanIndex++;
		}	
		return scanlist;
	}
}
