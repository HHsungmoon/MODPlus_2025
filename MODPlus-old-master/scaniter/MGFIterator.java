package scaniter;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.StringTokenizer;

public class MGFIterator extends ScanIterator {
		
	public MGFIterator( String fileName ) throws IOException {	
		super(fileName);
		
		BufferedReader in = new BufferedReader( new FileReader(fileName) );
		
		String buf;
		int size = 0;			
		while( (buf = in.readLine()) != null ){
			if( buf.startsWith("BEGIN IONS") )	// begining of new spectrum
				size++;
		}	
		in.close();
		sizeOfScans = size;
		scanIndex = 1;
		fin = new BufferedReader( new FileReader(fileName) );
	}
	
	public ArrayList<MSMScan> getNext() throws IOException {
		
		ArrayList<MSMScan> scanlist = new ArrayList<MSMScan>();
		MSMScan curScan = null;
		
		String buf;
		while( (buf = fin.readLine()) != null ) {

			if( buf.startsWith("BEGIN") ){	// begining of new spectrum	
				String title="";
				int scanNo = 0;
				int charge= 0;
				double pmz= -1;

				while( (buf = fin.readLine()) != null ) {	
					if( !buf.contains("=") ) break;				

					if( buf.startsWith("TITLE=") ){
						title = buf.substring(buf.indexOf("=")+1).trim();	
						title = title.replaceAll("[/,:*?\"<>|\\\\]", "_");
					}
					else if( buf.startsWith("CHARGE=") ){												
						int st = buf.lastIndexOf('=')+1;
						int ed = st+1;
						for(int i=st; i<buf.length(); i++){
							if( Character.isDigit( buf.charAt(i) ) ){
								st = i;
								ed = st+1;
								break;
							}
						}
						for(int i=ed; i<buf.length(); i++){
							if( !Character.isDigit( buf.charAt(i) ) ){
								ed = i;
								break;
							}
						}
						try {
							charge= Integer.parseInt( buf.substring(st,ed) );
						} catch (NumberFormatException e) {
							charge = 0;
						}
					}
					else if( buf.startsWith("PEPMASS=") ){						
						StringTokenizer token = new StringTokenizer(buf.substring(buf.indexOf("=")+1));
						if( token.countTokens() != 0 ) pmz = Double.parseDouble( token.nextToken() );
					}
					else if( buf.startsWith("SCANS=") ){						
						StringTokenizer token = new StringTokenizer(buf.substring(buf.indexOf("=")+1));
						if( token.countTokens() != 0 ) scanNo = Integer.parseInt( token.nextToken() );
					}
				}
				
				ArrayList<RawPeak> rawPL = new ArrayList<RawPeak>();				
				do {
					if( buf.startsWith("END") ) break;
					StringTokenizer token = new StringTokenizer(buf);
					if( token.countTokens() > 1 )
						rawPL.add( new RawPeak(Double.parseDouble(token.nextToken()), Double.parseDouble(token.nextToken())) );
				} while( (buf = fin.readLine()) != null );
				Collections.sort(rawPL);

				if( pmz > 0 ){
					if( charge > 0 ){			
						curScan = new MSMScan(title, scanIndex, scanNo, pmz, charge);
						if( curScan.setSpectrum(rawPL) ) scanlist.add(curScan);
					}
					else{
						for(int cs=MIN_ASSUMED_CHARGE; cs<=MAX_ASSUMED_CHARGE; cs++){						
							curScan = new MSMScan(title, scanIndex, scanNo, pmz, cs);
							if( curScan.setSpectrum(rawPL) ) scanlist.add(curScan);
						}
					}
				}
				scanIndex++;
				break;
			}	
		}
		return scanlist;
	}
	
}

