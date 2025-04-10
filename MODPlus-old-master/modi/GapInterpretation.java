package modi;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

import msutil.PGraph;

public class GapInterpretation extends ArrayList<PTMRun> {
	private		double gapScore[] = null;
	
	public double	getGapScore(int index)
	{
		if (gapScore!=null && gapScore.length>index) return gapScore[index];
		else return 0.0;
	}
	
	public class PTMRunScorePair {
		public PTMRun ptmRun;
		public double score;
	}
	
	public class PTMRunScorePairComparator implements Comparator<PTMRunScorePair>
	{
		public int compare(PTMRunScorePair run1, PTMRunScorePair run2)
		{
			double score1 = run1.score;
			double score2 = run2.score;
			if(score1 > score2)
				return -1;
		//	else if(score1 == score2)
		//		return 0;
			else if(score1 < score2)
				return 1;
			
			double delta1=0, delta2=0;
			for( PTMOccurrence pc : run1.ptmRun ){
				delta1 += Math.abs(pc.getPTM().getMassDifference());
			}
			for( PTMOccurrence pc : run2.ptmRun ){
				delta2 += Math.abs(pc.getPTM().getMassDifference());
			}
			if(delta1 > delta2) return 1;
			else if(delta1 < delta2) return -1;
			else return 0;
		}
		
		public boolean equals(PTMRun run1, PTMRun run2)
		{
			return run1.score == run2.score;
		}
	}	

	
	public GapInterpretation( Gap gap, ArrayList<PTMRun> interpretation, PGraph graph )
	{
		ArrayList<PTMRunScorePair> list = new ArrayList<PTMRunScorePair>();
		
		for(PTMRun run : interpretation)
		{
			PTMRunScorePair pair = new PTMRunScorePair();
			pair.ptmRun = run;
			pair.score = gap.getScore(run, graph);
			list.add(pair);
		}
		Collections.sort(list, new PTMRunScorePairComparator());
		
		this.gapScore = new double[interpretation.size()];	
		for (int i=0; i<interpretation.size() && i<Constants.maxInterpretationPerGap; i++){
			this.add(list.get(i).ptmRun);
			this.gapScore[i] = list.get(i).score;
		}
	}
	
}
