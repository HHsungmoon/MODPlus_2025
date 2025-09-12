package modi;

import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.Collections;
import java.util.StringTokenizer;

public class ScanCap implements Comparable<ScanCap> {
	private String 	title;
	private double 	pmz;
	private double 	neutralMW;
	private int 	charge;
	private int 	scanNo;
	private long 	offset;

	public ScanCap(String title, double pmz, int charge){
	
		this.title 		= title;	
		this.pmz 		= pmz;
		this.charge 	= charge;
		this.neutralMW 	= (pmz - Constants.Proton)*charge;
	}
	
	public ScanCap(String title, int sn, double pmz, int charge){		
		this.title 		= title;	
		this.scanNo		= sn;
		this.pmz 		= pmz;
		this.charge 	= charge;
		this.neutralMW 	= (pmz - Constants.Proton)*charge;
	}

	public void setOffset(long offset){ this.offset = offset; }
	public String getTitle(){ return title; }
	public double getObservedMW(){ return neutralMW; }
	public double getPMZ(){ return pmz; }
	public int getCharge(){ return charge; }
	public long getOffset(){ return offset; }
	
	public int compareTo(ScanCap s) 
	{
		if( this.neutralMW > s.neutralMW ) return 1;
		else if( this.neutralMW < s.neutralMW ) return -1;
		
		if( this.charge > s.charge ) return 1;
		else if( this.charge < s.charge ) return -1;
		else return 0;
	}

}
