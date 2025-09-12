package processedDB;

import modi.Constants;

public class SequenceTag implements Comparable<SequenceTag> {
	private double 	nGap;
	private short	mStart, mEnd;
	private byte	mType;//0:B, 1:Y
	
	public SequenceTag(double n, int start, int end, int type){
		nGap= n;
		mStart= (short)start;
		mEnd  = (short)end;
		mType = (byte)type;
	}
	
	public double getNGap(){ return nGap; }
	public int getStart(){ return mStart; }
	public int getEnd(){ return mEnd; }
	public int getType(){ return mType; }

	public void shiftTag(double gap, int offset){
		nGap -= gap;
		mStart += offset;
		mEnd += offset;
	}
	
	public boolean extendable(SequenceTag x){
		if( this.mType == x.mType ) {
			if( Math.abs(this.nGap-x.nGap) <= Constants.massToleranceForDenovo ){
					return true;
			}		
		}
		return false;
	}
	
	public boolean extend(SequenceTag x){
		if( extendable(x) ){
			this.mEnd= (short)Math.max( this.mEnd, x.mEnd );
			return true;
		}	
		return false;
	}
	
	public int compareTo(SequenceTag x) {
		if( x.mStart < this.mStart ) return 1;
		else if( x.mStart > this.mStart ) return -1;
		else return 0;
	}
	
}


