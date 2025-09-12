package modi;

public class PTMOccurrence {
	int	position;	// position by index in peptide
	PTM	ptm;
	
	public PTMOccurrence( int position, PTM ptm )
	{
		this.position = position;
		this.ptm = ptm;
	}

	public	int		getPosition()	{ return position; }
	public	PTM		getPTM()		{ return ptm; }
	
	public String toString()
	{
		return ptm + "(" + Constants.getString(ptm.getMassDifference()) + ")," + position;
	}
	
	public String toString(int basePos)
	{
		return ptm + "(" + Constants.getString(ptm.getMassDifference()) + ")," + (basePos+position);
	}
}
