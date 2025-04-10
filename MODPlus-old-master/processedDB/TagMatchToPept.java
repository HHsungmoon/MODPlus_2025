package processedDB;

import modi.IonDirection;
import modi.Tag;
import modi.Constants;

public class TagMatchToPept implements Comparable<TagMatchToPept> {
	
	Tag 			matchedTag;
	double 			nGap, cGap;
	int 			staSite, endSite;	
	IonDirection 	ir;
	
	public TagMatchToPept(Tag tag, double n, double c, int start, int end, IonDirection x){
		matchedTag= tag;
		nGap= n;
		cGap= c;
		staSite= start;
		endSite= end;
		ir= x;
	}

	public int compareTo(TagMatchToPept x){	// default comparator : mass
		if( x.staSite < this.staSite ) return 1;
		else if( x.staSite > this.staSite ) return -1;
		else return 0;
	}	
	
}


