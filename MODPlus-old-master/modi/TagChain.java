package modi;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.ListIterator;
import java.util.Map;
import java.util.TreeSet;

public class TagChain extends TreeSet<SpecInterpretation> implements Comparable<TagChain> {
	Peptide		matchedPeptide;
	Spectrum	sourceSpectrum;
	double 		score = 0;
	boolean		allGapAnnotated = true;
	double 		mostAbundantBYPeakIntensity;
	int 		longestTagLen= 0; 
	int 		tagHit= 0;
	int 		tagCoverage= 0;

	public Peptide	getMatchedPeptide()			{ return matchedPeptide; }
	
	public	TagChain(Peptide matchedPeptide, Spectrum sourceSpectrum)
	{
		this.matchedPeptide = matchedPeptide; 
		this.sourceSpectrum = sourceSpectrum;
	}
	
	public	void	setAllGapAnnotated(boolean annotated)	{ allGapAnnotated = annotated; } 

	public	int	 compareTo(TagChain tc)	
	{ 
		double score1 = this.getScore();
		double score2 = tc.getScore();
		
		if(score1 < score2)
			return 1;
		else if(score1 == score2)
			return 0;
		else
			return -1;
	}
	
	private void setMostAbundantBYPeakIntensity()
	{
		ArrayList<Peak> theoPeaks = new ArrayList<Peak>();
		for(SpecInterpretation t : this)
			theoPeaks.addAll(t.getTheoreticalPeaks());
		PeakListComparator comparator = new PeakListComparator(sourceSpectrum, theoPeaks);
		for(PeakPair pp : comparator.getSharedPeaks())
		{
			if(pp.getSecond().getPeakProperty() == PeakProperty.B_ION ||
					pp.getSecond().getPeakProperty() == PeakProperty.Y_ION)
			{
				if(this.mostAbundantBYPeakIntensity <  pp.getFirst().getIntensity())
					this.mostAbundantBYPeakIntensity = pp.getFirst().getIntensity();
			}
		}
	}
	
	public boolean makeGap() {
		int peptLastIndex = matchedPeptide.size()-1;

		setMostAbundantBYPeakIntensity();
		
		// check if matchedPeptide is Protein N-Term/C-term : should move to Peptide
		boolean protNTerm = false, protCTerm = false;
		
		if( matchedPeptide.getSrcProteinInfo().get(0).getSrcProteinID() == 0 ) protNTerm = true;
		else if( matchedPeptide.getSrcProteinInfo().get(0).getSrcProteinID() == 2 ) protCTerm = true;

		int start = 0, end = peptLastIndex;
		ArrayList<Gap> gapList = new ArrayList<Gap>();
		PTMPosition position;
		double prevOffset = Constants.NTERM_FIX_MOD;
		IonDirection prevDirection = IonDirection.B_DIRECTION;
		double bStart = Constants.B_ION_OFFSET;
		double yStart;
		boolean hasQTag = false;
		for(SpecInterpretation element : this)
		{
			assert(element instanceof MatchedTag);
			MatchedTag curTag = (MatchedTag)element;
			if( curTag.size() > Constants.minTagLengthPeptideShouldContain ) hasQTag = true;
			double motherMass = curTag.getSourceSpectrum().getCorrectedMW();
			end = element.getStart()-1;
			
			if(curTag.getDirection() == IonDirection.B_DIRECTION)
				yStart = curTag.get(0).getComplementMass(motherMass);
			else
				yStart = curTag.get(curTag.size()-1).getMass();
			
			if(start <= end)
			{
				// determine gap position
				if(start == 0 && protNTerm )
					position = PTMPosition.PROTEIN_N_TERM;
				else if(start == 0 && !protNTerm )
					position = PTMPosition.ANY_N_TERM;
				else
					position = PTMPosition.ANYWHERE;
				
				double offset = curTag.getNTermOffset() - prevOffset;
				
				Gap tpGap = new Gap(matchedPeptide, start, end, bStart, yStart, offset, 
						position, sourceSpectrum, this.mostAbundantBYPeakIntensity);
				gapList.add(tpGap);
				
				prevOffset = curTag.getNTermOffsetByLastPeak();
				prevDirection = curTag.getDirection();
			}
			if(curTag.getDirection() == IonDirection.B_DIRECTION)
				bStart = curTag.get(curTag.size()-1).getMass();
			else
				bStart = curTag.get(0).getComplementMass(curTag.getSourceSpectrum().getCorrectedMW());
			start = element.getEnd() + 1;
		}
		if( !hasQTag ) return false;
		
		// make C-term gap 
		if(start <= peptLastIndex)
		{
			end = peptLastIndex;
			if(protCTerm)
				position = PTMPosition.PROTEIN_C_TERM;
			else
				position = PTMPosition.ANY_C_TERM;
			
			Gap tpGap = new Gap(matchedPeptide, start, end, bStart, Constants.Y_ION_OFFSET, 
					((MatchedTag)this.last()).getCTermOffset()-Constants.CTERM_FIX_MOD, position, sourceSpectrum, this.mostAbundantBYPeakIntensity);
			gapList.add( tpGap );
		}
		
		for(Gap gap : gapList){	
			if( !gap.isReasonable() )
				return false;
			this.add(gap);
		}
		return true;
	}

	public void setTagChainScore() {
		
		HashSet<Integer> modtype = new HashSet<Integer>();
		
		double gapScore= 0;
		this.score= 0;
		for( SpecInterpretation t : this ){
			if( t instanceof Gap ) {
				Gap gap= (Gap)t;
				if( Math.abs( gap.getOffset() ) > Constants.gapTolerance ) {
					modtype.add( Constants.round(gap.getOffset()) );
				}
			}
			if( t instanceof MatchedTag ){
				this.score -= 1;	
				
				MatchedTag tag= (MatchedTag)t;
				tagCoverage += tag.sequence().size();
							
				for( Peak p: tag ) this.score += p.getProbability();			
			}
		}
		
		this.score -= modtype.size();
		
		if( this.tagHit == 1 && this.longestTagLen < Constants.minTagLengthPeptideShouldContain ){
			if( modtype.size() == 0 ) {			
				this.score = 0; //because this is unmodified
			}
		}//*/
	}
	
	public double getScore(){		// should be optimized later
		return score;
	}
	
	public static LinkedList<TagChain> buildTagChainList(Map.Entry<Peptide, LinkedList<MatchedTag>> tagPool){
		initTagMerge(tagPool);
		return enumTagChain(tagPool);
	}
	
	public String getTagChainSequence()
	{
		StringBuffer output = new StringBuffer();
		for(SpecInterpretation t : this)
			output.append(t.getSequenceStr());

		return output.toString();
	}
	
	private static void initTagMerge( Map.Entry<Peptide, LinkedList<MatchedTag>> entry )
	{
		LinkedList<MatchedTag> tagList = entry.getValue();
		Collections.sort(tagList, new SpecInterpretationComparator());	
		Spectrum ccspec = tagList.getFirst().getSourceSpectrum();
	
		if( tagList.getFirst().start != 0 ){
			Peptide pept= entry.getKey();		
			Tag gto= ccspec.getB2Tag(pept.subSequence(0, 2));
			if( gto != null ){
				MatchedTag b2mTag= MatchedTag.getMatchedB2Tag(gto, pept);
				tagList.addFirst( b2mTag );
			}		
		}//*/	

		int initSize  = tagList.size();
		for(int i=0; i<initSize-1; i++)
		{	
			MatchedTag seed = tagList.get(i);
			for(int j=i+1; j<initSize; j++)
			{
				MatchedTag tag = tagList.get(j);		
				if( seed.getRelativePosition(tag) == RelativePosition.SEPERATED )
					break;
				
				if( seed.getDirection() == tag.getDirection() && 
						Constants.fEqual( seed.getOffset(), tag.getOffset() ) ){
					MatchedTag temp= new MatchedTag(seed);
					temp.extend(tag);
					if( temp.size() != temp.sequence().size()+1 )
						continue;
					seed.extend(tag);	// extend seed from extendable TagList
					tagList.remove(j);		
					j--;
					initSize--;
					continue;
				}
				
				if( seed.isComplementarySame(tag) ){
					tagList.remove(j);		// remove seed from matchedTagList
					j--;
					initSize--;		
					continue;
				}						
			}
		}	
		
		int countOfLongTag= 0;
		MatchedTag sLongTag= null;
		for(int k=0; k<tagList.size(); k++){
			tagList.get(k).setScore();
			if( tagList.get(k).sequence().size() > 2 ){
				countOfLongTag++;
				sLongTag= tagList.get(k);
			}
		}
		
//		topept++;
		if( tagList.size() > Constants.maxTagPerPept ){
		//	maxHitperPept++;
			Collections.sort(tagList, Collections.reverseOrder(new TagComparator()) );
			for(int k=Constants.maxTagPerPept; k<tagList.size(); k++){
				tagList.remove(k--);
			}
			Collections.sort(tagList, new SpecInterpretationComparator());	
		}//*/
		
		if( countOfLongTag == 1 && tagList.size() < 4 ){
			int rev= 0;
			MatchedTag newTag= new MatchedTag(sLongTag);
			if( sLongTag.getDirection() == IonDirection.B_DIRECTION ){				
			//	if( !Constants.fEqual(sLongTag.getNTermOffset(), 0) )
				{
					rev++;
					int s;
					for(s=1; sLongTag.sequence().size()-rev > 2 ;s++){
						if( ccspec.isConfidentPeak(sLongTag.get(s), 1) ) break;
						rev++;
					}
					newTag.trim(newTag.getStart(), s);
				}
			//	if( !Constants.fEqual(sLongTag.getCTermOffset(), 0) )
			
				{
					rev++;
					int s;
					for(s=1; sLongTag.sequence().size()-rev > 2 ;s++){
						if( ccspec.isConfidentPeak(sLongTag.get(sLongTag.size()-s-1), 1) ) break;
						rev++;
					}
					newTag.trim(newTag.getEnd(), s);
				}
			}
			else{
			//	if( !Constants.fEqual(sLongTag.getCTermOffset(), 0) )
				{
					rev++;
					int s;
					for(s=1; sLongTag.sequence().size()-rev > 2 ;s++){
						if( ccspec.isConfidentPeak(sLongTag.get(s), 1) ) break;
						rev++;
					}
					newTag.trim(newTag.getEnd(), s);
				}
			//	if( !Constants.fEqual(sLongTag.getNTermOffset(), 0) )
				{
					rev++;
					int s;
					for(s=1; sLongTag.sequence().size()-rev > 2 ;s++){
						if( ccspec.isConfidentPeak(sLongTag.get(sLongTag.size()-s-1), 1) ) break;
						rev++;
					}
					newTag.trim(newTag.getStart(), s);
				}
			}
			if( rev > 0 ){
				if( newTag.size() > 2 ) tagList.add(newTag);		
			}
		}
	}
	
	private static LinkedList<TagChain> enumTagChain(Map.Entry<Peptide, LinkedList<MatchedTag>> entry)
	{				
		LinkedList<MatchedTag> tagList = entry.getValue();
		Collections.sort(tagList, new SpecInterpretationComparator());	// sort matched tag list by start position
				
		LinkedList<TagChain> tagChainList = new LinkedList<TagChain>();		
		
	//	System.out.println("TAGCHAIN.JAVA : " + entry.getKey()+ " " + tagList.size() );		
		for( MatchedTag tag : tagList ) {	
			TagChain t = new TagChain(entry.getKey(), tag.getSourceSpectrum());
			t.add(tag);
			tagChainList.add(t);
	//		System.out.println(tag);
		}//generate base tag_chain		
	//	System.out.println(entry.getKey() + " " + tagList.size());	
		
		// tag chain enumeration
		ArrayList<TagChain> addedTagChain = new ArrayList<TagChain>();
		int start = 0, end = tagChainList.size(), addedCount=0;
		
		while(true)
		{
			ListIterator<TagChain> listIt = tagChainList.listIterator(start);
			while(listIt.hasNext())
			{
				TagChain curTC = listIt.next();
				if( !(curTC.last() instanceof MatchedTag) )
					continue;
				MatchedTag seed= (MatchedTag)curTC.last();
				
				for( MatchedTag tag : tagList ) {
					if( seed.compareTo(tag) > 0 || curTC.useSamePeak(tag) ) continue;	
		
					combineTagChains(addedTagChain, curTC, tag);
				}
			}
			addedCount = addedTagChain.size();
			tagChainList.addAll(addedTagChain);
			addedTagChain.clear();
			start = end;
			end += addedCount;
			if( addedCount == 0 ) break;
		}
		
		// make gap for all tag chain list
		ListIterator<TagChain> listIt = tagChainList.listIterator();
		double topSCore = 0;
		while( listIt.hasNext() ){
			TagChain curTC = listIt.next();		
			if( curTC.makeGap() ){ 
				curTC.setTagChainScore(); 
				if( curTC.score < 0 ) listIt.remove();
				if( curTC.score > topSCore )
					topSCore = curTC.score;
			}
			else{ listIt.remove(); }
		}
		
		listIt = tagChainList.listIterator();
		while( listIt.hasNext() ){
			TagChain curTC = listIt.next();
			if( curTC.score < topSCore * Constants.tagChainPruningRate ){
				listIt.remove();				
			}
		}
		
		if( tagChainList.size() > Constants.maxTagChainPerPept ){
			Collections.sort( tagChainList );
			for(int i=Constants.maxTagChainPerPept; i<tagChainList.size(); i++)
				tagChainList.remove(i--);
		}
		
		return tagChainList;
	}
	
	private static void  combineTagChains (ArrayList<TagChain> newTCList, TagChain baseTC, MatchedTag tag) {	
		if( !(baseTC.last() instanceof MatchedTag) )
			return;
		
		TagChain newTC = (TagChain)baseTC.clone();
		MatchedTag seed= (MatchedTag)baseTC.last();
		
		RelativePosition ir= seed.getRelativePosition(tag);
		if( ir == RelativePosition.SEPERATED ){			
			newTC.add(tag);
			newTCList.add(newTC);
		}
		
		else if( ir == RelativePosition.ADJACENT ){
			
			if(	seed.size() < 3 )
				return;

			if(	seed.sequence().size() < Constants.minTagLengthPeptideShouldContain &&
					tag.sequence().size()< Constants.minTagLengthPeptideShouldContain )
				return;

			newTC.remove(newTC.last());	

			if( LinkableTags(seed, tag) ){			
				if( seed.size() < tag.size() )
					newTCList.add( extendTagChain(newTC, seed, tag.getTrimedTag( tag.getStart(), 1 )) );
				else
					newTCList.add( extendTagChain(newTC, seed.getTrimedTag( seed.getEnd(), 1 ), tag) );
			}
			else{				
				if( tag.size() > 2 ) newTCList.add( extendTagChain(newTC, seed, tag.getTrimedTag( tag.getStart(), 1 )) );
				newTCList.add( extendTagChain(newTC, seed.getTrimedTag( seed.getEnd(), 1 ), tag) );
			}
		}
		
		else if( ir == RelativePosition.OVERLAP ){ 

			if( seed.sequence().size()< Constants.minTagLengthPeptideShouldContain &&
						tag.sequence().size()< Constants.minTagLengthPeptideShouldContain )
				return;

			newTC.remove(newTC.last());	
			
			int overLabPart= seed.getEnd() - tag.getStart() + 1;
				
			newTCList.add( extendTagChain(newTC, seed.getTrimedTag( seed.getEnd(), overLabPart ),
					tag.getTrimedTag( tag.getStart(), overLabPart )) );
			if( LinkableTags(seed, tag) ){ return; }
			
			if( tag.sequence().size() > overLabPart+1 )
				newTCList.add( extendTagChain(newTC, seed,
					tag.getTrimedTag( tag.getStart(), overLabPart+1 )) );
			if( seed.sequence().size() > overLabPart+1 )
				newTCList.add( extendTagChain(newTC, seed.getTrimedTag( seed.getEnd(), overLabPart+1 ),
					tag) );	
		}
		
		else if( ir == RelativePosition.INCLUDING ){
			
			if(	seed.sequence().size() < Constants.minTagLengthPeptideShouldContain &&
					tag.sequence().size() < Constants.minTagLengthPeptideShouldContain )
				return;
			if( Constants.fEqual( seed.getOffset(), tag.getOffset()) ) return;
			if( !seed.isLikelyChild(tag) ) return;
			
			int rev= 0;
			MatchedTag newTag;
			if( seed.score > tag.score )
				newTag= new MatchedTag(seed);
			else
				newTag= new MatchedTag(tag);

			if( seed.getDirection() == IonDirection.B_DIRECTION ){
				while( newTag.sequence().size() > 2 ) {
					if( newTag.sourceSpectrum.isConfidentPeak(newTag.get(0), 1) ) break;
					newTag.trim(newTag.getStart(), 1);
					rev++;
				}

				while( newTag.sequence().size() > 2 ) {
					if( newTag.sourceSpectrum.isConfidentPeak(newTag.get(newTag.size()-1), 1) ) break;
					newTag.trim(newTag.getEnd(), 1);
					rev++;
				}
			}
			else {
				while( newTag.sequence().size() > 2 ) {
					if( newTag.sourceSpectrum.isConfidentPeak(newTag.get(0), 1) ) break;
					newTag.trim(newTag.getEnd(), 1);
					rev++;
				}

				while( newTag.sequence().size() > 2 ) {
					if( newTag.sourceSpectrum.isConfidentPeak(newTag.get(newTag.size()-1), 1) ) break;
					newTag.trim(newTag.getStart(), 1);
					rev++;
				}
			}

			if( rev > 0 ){									
				newTC.remove(newTC.last());					
				newTC.add(newTag);
				newTCList.add(newTC);
			}
		}//INCLUDING*/
	}
	
	private static boolean LinkableTags(MatchedTag one, MatchedTag two){
		//ADJACENT tags citation
		if( one.getDirection() == two.getDirection() ) return false;
		if( Constants.fEqual( one.getNTermOffset(), two.getNTermOffset()) ) return true;
		return false;
	}

	private boolean useSamePeak(MatchedTag t) {
		for(SpecInterpretation element : this)
		{
			if(!(element instanceof MatchedTag)) continue;
			if( ((MatchedTag)element).useSamePeak(t) ) {
				return true;
			}
		}
		return false;
	}
	
	public boolean add(SpecInterpretation si){
	//	System.out.println(si);
		if( si instanceof MatchedTag ){
			tagHit++;
			if( longestTagLen < ((MatchedTag)si).sequence().size() )
				longestTagLen= ((MatchedTag)si).sequence().size();
		}			
		return super.add(si);
	}
	
	public boolean remove(SpecInterpretation si){
		if( si instanceof MatchedTag ){	tagHit--; }
		return super.remove(si);
	}
	
	private static TagChain extendTagChain (TagChain baseTC, MatchedTag one, MatchedTag two){
		TagChain newTC = (TagChain)baseTC.clone();
		newTC.add(one);
		newTC.add(two);		
		return newTC;
	}
	
	
	public ArrayList<PTMCombination> getPTMCombination(){
		
		ArrayList<PTMCombination> answers = new ArrayList<PTMCombination>();
		ArrayList<Gap> gapList = new ArrayList<Gap>();
		
		int madeSize = 1;
		for( SpecInterpretation t : this ){		
			if(t instanceof Gap){
				Gap gap = (Gap)t;
				if( gap.getGapInterpretation().size() == 0 ) continue;
				gapList.add(gap);
				madeSize *= gap.getGapInterpretation().size();
			}
		}
		
		if( gapList.size() == 0 ){
			answers.add( new PTMCombination( matchedPeptide.size() ) );
			return answers;
		}
		
		int poolSize = 1;
		int[] indexTrace = new int[gapList.size()];
		while( true ){

			PTMCombination ptmComb = new PTMCombination( matchedPeptide.size() );
			StringBuffer tempComb = new StringBuffer();
			int ix = 0, modCount = 0;
			for( Gap gap : gapList ){
				PTMRun prun = gap.getGapInterpretation().get(indexTrace[ix]);
				modCount += prun.size();
				tempComb.append( prun.toString(matchedPeptide, gap.getStart()) );
				prun.setPTMMass(ptmComb.ptms, gap.getStart());
				prun.setPTMs(ptmComb.ptmList, gap.getStart());
				ix++;
			}
			if( modCount <= Constants.maxPTMPerPeptide ) {
				ptmComb.ptmComb = tempComb.toString();
				answers.add(ptmComb);
			}
			if( poolSize++ ==  madeSize ) break;	
	
			indexTrace[--ix]++;
			while( ix > -1 && indexTrace[ix] == gapList.get(ix).getGapInterpretation().size() ){
				indexTrace[ix--] = 0;
				indexTrace[ix]++;
			}
		}
	
		return answers;
	}
	
	public String toString()
	{
		StringBuffer output = new StringBuffer(this.getTagChainSequence());
	
		output.append("(");
		for(SpecInterpretation t : this)
		{
			if(t instanceof Gap)
				output.append(new DecimalFormat("#.###").format(((Gap)t).getOffset())+ ",");
		}
		output.deleteCharAt(output.length()-1);
		output.append(") " + getScore());
		
		for (SpecInterpretation si : this)
		{
			if (!(si instanceof Gap)) continue;
			Gap gap = (Gap)si;
			if (gap.getGapInterpretation()==null || gap.getGapInterpretation().size()==0) continue;
			output.append("\nGap interpretation for "+gap.getStart()+"~"+gap.getEnd()+":");
			output.append(gap.getGapInterpretation());
		}
		
		return output.toString();
	}
}
