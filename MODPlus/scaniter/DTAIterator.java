package scaniter;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.StringTokenizer;

import modi.Constants;

public class DTAIterator extends ScanIterator {

	public DTAIterator( String fileName ) throws IOException {
		
		super(fileName);
		
		BufferedReader in = new BufferedReader( new FileReader(fileName) );	
		String buf;
		int size = 0;			
		boolean newDTA = true;
		while( (buf = in.readLine()) != null ) {
			StringTokenizer token = new StringTokenizer(buf);

			if( token.countTokens() != 2 ){
				newDTA = true;
				continue;
			}
			
			if( newDTA && token.countTokens() == 2 ) {
				size++;
				newDTA = false;
			}
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
		boolean newDTA = true;
		while( (buf = fin.readLine()) != null ) {			
			StringTokenizer token = new StringTokenizer(buf);
			if( token.countTokens() != 2 ){
				newDTA = true;
				continue;
			}
			
			if( newDTA && token.countTokens() == 2 ){
				double precursor = Double.parseDouble(token.nextToken());
				int charge = (int)Double.parseDouble(token.nextToken());
				precursor = (precursor+(charge-1)*Constants.Proton)/charge;
				
				ArrayList<RawPeak> rawPL = new ArrayList<RawPeak>();				
				while( (buf = fin.readLine()) != null ) {
					token = new StringTokenizer(buf);
					if( token.countTokens() < 2 ) break;
					rawPL.add( new RawPeak(Double.parseDouble(token.nextToken()), Double.parseDouble(token.nextToken())) );
				}
				Collections.sort(rawPL);
				
				curScan = new MSMScan( scanIndex++, precursor, charge );				
				if( curScan.setSpectrum(rawPL) ) scanlist.add(curScan);
				break;
			}
		}	
		return scanlist;
	}
}





















