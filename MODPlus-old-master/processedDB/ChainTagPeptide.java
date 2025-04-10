package processedDB;

import java.util.ArrayList;
import java.util.Collections;

import modi.Constants;

public class ChainTagPeptide extends MODPeptide {
	
	private static double	a = ( Constants.maxModifiedMass < 0 )?  Constants.gapTolerance : Constants.maxModifiedMass+Constants.gapTolerance; 
	private static double	b = ( Constants.minModifiedMass > 0 )? -Constants.gapTolerance : Constants.minModifiedMass-Constants.gapTolerance; 
	private static double	shiftWindow = (a - b);
	
	ArrayList<SequenceTag> mTags;

	public ChainTagPeptide(int s, int e, SequenceTag tag){
		super(s, e);
		mTags = new ArrayList<SequenceTag>();
		mTags.add(tag);
	}
	
	public ArrayList<SequenceTag> getMatchedTags( ){ return mTags; }

	public boolean extend( TagPeptide merged, TagPeptide toMerge, TagTrie trie ){ //multi-mod
		
		if( merged.pStart <= toMerge.pStart && toMerge.pStart < merged.pEnd )
		{
			double offset = msutil.MSMass.getPepMass( trie.getPeptide(merged.pStart, toMerge.pStart) );			
			if( offset <= shiftWindow && mTags.size() < 100 ) {
				this.pLeft = Math.max( this.pLeft,  toMerge.pLeft );
				this.pRight= Math.min( this.pRight, toMerge.pRight );
				this.pEnd  = toMerge.pEnd;
				
				SequenceTag tag = toMerge.mTag;
				tag.shiftTag( getNOffset(merged, trie)+offset, (toMerge.pStart - this.pStart) );
				this.mTags.add(tag);	
				
				return true;
			}
		}
		return false;
	}
	
	public double getNOffset( TagPeptide xp, TagTrie trie ){
		return msutil.MSMass.getPepMass( trie.getPeptide(this.pStart, xp.pStart) );		
	}
	
	public void arrangeTags(){
		Collections.sort( mTags );	
		int initSize = mTags.size();
		for(int i=0; i<initSize-1; i++){					
			for(int j=i+1; j<initSize; j++){
				if( mTags.get(i).getEnd()+1 < mTags.get(j).getStart() ) break;
				if( mTags.get(i).extend( mTags.get(j) ) ){
					mTags.remove(j);		
					j--;
					initSize--;
					continue;
				}
			}
		}			
	}
}
