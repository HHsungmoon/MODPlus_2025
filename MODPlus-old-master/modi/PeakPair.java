package modi;

public class PeakPair {
	public Peak first;
	public Peak second;
	
	public PeakPair(Peak first, Peak second)
	{
		this.first = first;
		this.second = second;
	}
	public	Peak	getFirst()	{ return first; }
	public	Peak	getSecond()	{ return second; }

	public	String	toString()
	{
		return new String("(" + first + "," + second + ")");
	}
}
