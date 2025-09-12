package modi;

public class DoublePair {
	public double x;		// ! Must x<y
	public double y;
	
	public DoublePair(double x, double y)
	{
		this.x = x;
		this.y = y;
	}
	boolean inRange( double offset )
	{
		return offset>x && offset<y;
	}
	
	double getRange()		{ return y-x; }
}
