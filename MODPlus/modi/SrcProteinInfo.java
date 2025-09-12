package modi;

public class SrcProteinInfo implements Comparable<SrcProteinInfo> {
	private int srcProteinID; // 0:nterm, 1, 2: cterm
	private int startPos, endPos;
	private int ntt;
	
	public SrcProteinInfo( int srcProteinID, int startPos, int endPos ) {
		this.srcProteinID	= srcProteinID;
		this.startPos		= startPos;
		this.endPos			= endPos;
	}
	
	public SrcProteinInfo( int srcProteinID, int startPos, int endPos, int ntt_v ) {
		this(srcProteinID, startPos, endPos);
		ntt = ntt_v;
	}
	
	public int			getSrcProteinID()	{ return srcProteinID; }
	public int			getStartPos()		{ return startPos; }
	public int			getEndPos()			{ return endPos; }
	public int			getNTT()			{ return ntt; }

	public String		toString() {
		return "(SRC:\""+srcProteinID+"\""+startPos+"~"+endPos+")";
	}

	public int compareTo(SrcProteinInfo x) {
		if( this.ntt > x.ntt ) return -1;
		else if( this.ntt < x.ntt ) return 1;		
		else {
			if( this.srcProteinID > x.srcProteinID ) return 1;
			else if( this.srcProteinID < x.srcProteinID ) return -1;
			return 0;
		}
	}
	
}
