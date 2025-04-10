package modi;

import java.util.ArrayList;
import java.util.Comparator;

public interface SpecInterpretation extends Comparable<SpecInterpretation> {
	public	Sequence	sequence();
	public	Peptide		getMatchedPeptide();
	public	int			getStart();
	public	int			getEnd();
	public	String		getSequenceStr();
	public	double		getScore();
	public	int			compareTo(SpecInterpretation t);
	public ArrayList<Peak> getTheoreticalPeaks();

}

class SpecInterpretationComparator implements Comparator<SpecInterpretation> {
	public	int	compare(SpecInterpretation t1, SpecInterpretation t2)
	{
		return t1.getStart() - t2.getStart();
	}
	
	public boolean equals(SpecInterpretation t1, SpecInterpretation t2)
	{
		return (t1.getStart() == t2.getStart());
	}
}