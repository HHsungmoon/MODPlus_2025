package scaniter;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;

import modi.Constants;

public abstract class ScanIterator {
	
	public static int MIN_ASSUMED_CHARGE = 2, MAX_ASSUMED_CHARGE = 3;
	public int sizeOfScans;
	public int scanIndex;
	public String fileName;
	public String baseName;
	public BufferedReader fin;
	
	public ScanIterator(String fName) throws IOException{
		fileName = fName;
		baseName = fileName.replace('\\', '/');
		baseName = baseName.substring(baseName.lastIndexOf('/')+1, baseName.lastIndexOf('.'));
		sizeOfScans = 0;
		scanIndex = 0;
		fin = null;
	}
	
	public String getFileName(){ return fileName; }
	public int size(){ return sizeOfScans; }
	public int getIndex(){ return scanIndex; }
	
	public boolean hasNext() throws IOException {
		if( scanIndex <= sizeOfScans ) return true;
		if( fin != null ) fin.close();		
		return false;
	}
	
	public abstract ArrayList<MSMScan> getNext() throws IOException; 	
	
	public static ScanIterator get(String specFile) throws IOException{
		System.out.print( "Reading MS/MS spectra.....  " );
		
		ScanIterator scaniter = null;
		if( specFile.toLowerCase().endsWith(".pkl") ){
			Constants.SPECTRA_FILE_TYPE = Constants.spectra_format.PKL;
			scaniter = new PKLIterator( specFile );
		}
		else if( specFile.toLowerCase().endsWith(".mgf") ){
			Constants.SPECTRA_FILE_TYPE = Constants.spectra_format.MGF;
			scaniter = new MGFIterator( specFile );
		}
		else if( specFile.toLowerCase().endsWith(".ms2") ){
			Constants.SPECTRA_FILE_TYPE = Constants.spectra_format.MS2;
			scaniter = new MS2Iterator( specFile );
		}
		else if( specFile.toLowerCase().endsWith(".dta") ){
			Constants.SPECTRA_FILE_TYPE = Constants.spectra_format.DTA;
			scaniter = new DTAIterator( specFile );
		}
		else if( specFile.toLowerCase().endsWith(".mzxml") ){
			Constants.SPECTRA_FILE_TYPE = Constants.spectra_format.MZXML;
			scaniter = new MZXMLIterator( specFile );
		}		
		return scaniter;
	}
	
	public static ScanIterator get(String specFile, Constants.spectra_format type) throws IOException{
		System.out.print( "Reading MS/MS spectra.....  " );
		
		ScanIterator scaniter = null;
		if( Constants.spectra_format.PKL == type ){
			scaniter = new PKLIterator( specFile );
		}
		else if( Constants.spectra_format.MGF == type ){
			scaniter = new MGFIterator( specFile );
		}
		else if( Constants.spectra_format.MS2 == type ){
			scaniter = new MS2Iterator( specFile );
		}
		else if( Constants.spectra_format.DTA == type ){
			scaniter = new DTAIterator( specFile );
		}
		else if( Constants.spectra_format.MZXML == type ){
			scaniter = new MZXMLIterator( specFile );
		}		
		return scaniter;
	}
}
