package modi;

import java.util.ArrayList;

public class Sequence extends ArrayList<AminoAcid> implements Comparable<Sequence> {
	// override of ArrayList<AminoAcid> get method
	public AminoAcid get(int i)
	{
		if(i <= -1)	// N-terminal
			return AminoAcid.getAminoAcid(21);
		else if(i >= this.size())	// C-terminal
			return AminoAcid.getAminoAcid(22);
		else
			return super.get(i);
	}
	public String toString()
	{
		StringBuffer output = new StringBuffer();
		for (AminoAcid aa : this)
		{
			output.append( aa.getResidue() );
		}
		return output.toString();
	}

	public String toString(IonDirection direction)
	{
		StringBuffer output = new StringBuffer();
		if(direction == IonDirection.B_DIRECTION)
		{
			return toString();
		}
		else
		{
			for (AminoAcid aa : this)
			{
				output.append(Character.toLowerCase(aa.getResidue()));
			}
			return output.toString();
		}
	}
	public int compareTo(Sequence o) {
		for (int i=0; i<size() && i<o.size(); i++)
		{
			int r = get(i).compareTo(o.get(i));
			if (r<0) return -1;
			else if (r>0) return 1;
		}
		int r = size()-o.size();
		if (r<0) return -1;
		else if (r>0) return 1;
		else return 0;
	}

	
	public double getMonoMass()
	{
		double sum = 0.;
		for(int i=0; i<this.size(); i++)
			sum += this.get(i).getMonoMass();
		return sum;
	}
	
	public double getMonoMass(int from, int to)	// from : inclusive, to : exclusive
	{
		if( from < 0 || to > this.size() ) return 0.;
		double sum = 0.;
		for(int i=from; i<to; i++)
			sum += this.get(i).getMonoMass();
		
		return sum;
	}
	
	public Sequence subSequence(int from, int to)	// from : inclusive, to : exclusive
	{
		if(from < 0 || to > this.size())	// error
			return null;
		Sequence newSeq = new Sequence();
		for(int i=from; i<to; i++)
			newSeq.add(this.get(i));
		return newSeq;
	}
	
	public static Sequence getSequence(String seq)
	{
		Sequence retSeq = new Sequence();
		for(int i=0; i<seq.length(); i++)
		{
			AminoAcid aa = AminoAcid.getAminoAcid(seq.charAt(i));
			if(aa == null)
				return null;
			else
				retSeq.add(aa);
		}
		return retSeq;
	}
}
