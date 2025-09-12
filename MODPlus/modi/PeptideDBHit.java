package modi;

import java.io.Serializable;

public class PeptideDBHit implements Serializable {
	private Peptide 		hitPeptide;
	private int 			start;
	private int 			end;
	private double			nTermOffset;
	private double			cTermOffset;
	
	public int			getStart()			{ return start; }
	public int			getEnd()			{ return end; }
	public Peptide		getHitPeptide()		{ return hitPeptide; }
	
	public void			setStart(int start)	{ this.start = start; }
	public void			setEnd(int end)		{ this.end = end; }
	public boolean		isNtermIon()		{ return start == 0; }
	public boolean		isCtermIon()		{ return end == hitPeptide.size()-1; }
	
	public PeptideDBHit( Peptide hitPeptide, int start, int end, double nTermOffset, double cTermOffset )
	{
		this.hitPeptide		= hitPeptide;
		this.start			= start;
		this.end			= end;
		this.nTermOffset	= nTermOffset;
		this.cTermOffset	= cTermOffset;		
	}

}
