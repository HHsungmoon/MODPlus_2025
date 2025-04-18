package modi;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.ListIterator;
import java.util.Map;
import java.util.TreeMap;

import msutil.PGraph;
import msutil.Scoring;

public class TagChainPool extends TreeMap<Peptide, LinkedList<TagChain>> {

	// 추가: 스레드별 TagChainPool 생성 도우미
	public static TagChainPool buildTagChainPoolSafely(MatchedTagPool matchedTags) {
		TagChainPool localPool = new TagChainPool();

		if (matchedTags == null || matchedTags.size() == 0)
			return localPool;

		for (Map.Entry<Peptide, LinkedList<MatchedTag>> entry : matchedTags.entrySet()) {
			LinkedList<TagChain> tcList = TagChain.buildTagChainList(entry);
			if (tcList.size() == 0) continue;
			Peptide matchedPeptide = new Peptide(entry.getKey());
			localPool.put(matchedPeptide, tcList);
		}

		return localPool;
	}

	public void mergeWith(TagChainPool other) {
		for (Map.Entry<Peptide, LinkedList<TagChain>> entry : other.entrySet()) {
			Peptide pep = entry.getKey();
			if (!this.containsKey(pep)) {
				this.put(pep, entry.getValue());
			} else {
				this.get(pep).addAll(entry.getValue()); // 합치기
			}
		}
	}


	public void discardPoorTagChain()	// should be modified
	{
		Iterator<Map.Entry<Peptide, LinkedList<TagChain>>> it = this.entrySet().iterator();

		// calculating best score
		double bestScore = 0;
		TagChain bet = null;
		while( it.hasNext() ) {
			for(TagChain tc : it.next().getValue())
			{				
				double curScore = tc.getScore();
				if( bestScore < curScore ){
					bestScore = curScore;
					bet = tc;
				}
			}
		}
	
		// discarding tag chain		
		it = this.entrySet().iterator();
		ArrayList<Peptide> deletePepList = new ArrayList<Peptide>();
		
		while(it.hasNext())
		{
			LinkedList<TagChain> tcList = it.next().getValue();
			Peptide pep = tcList.getFirst().getMatchedPeptide();
			ListIterator<TagChain> listIt = tcList.listIterator();
			while(listIt.hasNext()) {
				TagChain tc = listIt.next();
				if( tc.getScore() <= bestScore * Constants.tagChainPruningRate ){
					listIt.remove();
				}	
			}
			
			if( tcList.size() == 0 )
				deletePepList.add(pep);
		}
		
		for( Peptide pep : deletePepList )
			this.remove(pep);	
	}
	
	public String toString()
	{
		StringBuffer output = new StringBuffer();
		output.append("Matched Tag Pool");
		Iterator<Map.Entry<Peptide, LinkedList<TagChain>>> it = this.entrySet().iterator();
		while(it.hasNext())
		{
			Map.Entry<Peptide, LinkedList<TagChain>> entry = it.next();
			output.append("\n").append(entry.getKey()).append("\n");
			for(TagChain tc : entry.getValue())
			{
				output.append(tc).append("\n");
			}
		}
		return output.toString();
	}
	
	public ArrayList<AnsPeptide> getAnswerPeptides( PGraph graph ){		
		AnsHeap answerPepts = new AnsHeap();
		
		Iterator<Map.Entry<Peptide, LinkedList<TagChain>>> entries = this.entrySet().iterator();
		while(entries.hasNext()){
			Map.Entry<Peptide, LinkedList<TagChain>> entry = entries.next();
			String pept = entry.getKey().toString();
			LinkedList<TagChain> alignedTagChainList = entry.getValue();
			if( alignedTagChainList.size() < 1 ) continue;
			
			HashSet<PTMCombination> ptmComb = new HashSet<PTMCombination>();
			for( TagChain tc : alignedTagChainList ){
				if( tc.allGapAnnotated ) ptmComb.addAll( tc.getPTMCombination() );
			}
			
			Iterator<PTMCombination> iter = ptmComb.iterator();
			while( iter.hasNext() ){
				PTMCombination p = iter.next();
	
				int s = Scoring.getModEyeRankScore(pept, p.ptms, graph);
				if( s < 0 ) continue;
				AnsPeptide candidate = new AnsPeptide(entry.getKey(), p.ptmComb, p.ptms, p.ptmList, s);
				answerPepts.add( candidate ); 
			}
		}
		
		return answerPepts.getFinalList(graph);
	}
	
}
	
