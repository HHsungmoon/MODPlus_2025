package processedDB;
import java.util.Comparator;

public class TagPeptide extends MODPeptide {
	
	SequenceTag 	mTag;	
	public TagPeptide(int s, int e, SequenceTag tag){
		super(s, e);
		mTag   = tag;
	}
}

class TagPeptComparator implements Comparator<TagPeptide>{
	public int compare(TagPeptide x1, TagPeptide x2){
		if( x1.pStart > x2.pStart  ) return 1;
		else if( x1.pStart < x2.pStart ) return -1;

		if( x1.pEnd > x2.pEnd ) return 1;
		else if( x1.pEnd < x2.pEnd ) return -1;
			
		if( x1.mTag.getNGap() < x2.mTag.getNGap() ) return 1;
		else if( x1.mTag.getNGap() > x2.mTag.getNGap() ) return -1;
		
		return 0;
	}	
	public boolean equals(TagPeptide x1, TagPeptide x2){
		return x1 == x2;
	}
}






