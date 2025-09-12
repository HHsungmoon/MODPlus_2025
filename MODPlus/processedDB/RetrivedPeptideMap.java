package processedDB;

import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import modi.Constants;
import modi.IonDirection;
import modi.MatchedTag;
import modi.MatchedTagPool;
import modi.Peptide;
import modi.PeptideDBHit;
import modi.TagPool;

public class RetrivedPeptideMap extends HashMap<String, RetrivedPeptide> {

	public void add(RetrivedPeptide pept){
		RetrivedPeptide rp = this.get(pept.getSequence());
		
		if ( rp == null ) this.put(pept.getSequence(), pept);
		else rp.addProtein( pept.getProteinList().get(0) );	
	}
	
	public void combine(RetrivedPeptideMap rpMap){		
		Iterator<Map.Entry<String, RetrivedPeptide>> it = rpMap.entrySet().iterator();
		while(it.hasNext())
		{
			Map.Entry<String, RetrivedPeptide> entry = it.next();			
			RetrivedPeptide rp = this.get(entry.getKey());			
			if ( rp == null ) this.put( entry.getKey(), entry.getValue() );
			else rp.addTag( entry.getValue().getTagList().get(0) );			
		}		
	}
	
	public MatchedTagPool convertToMatchedTagPool(TagPool shortTags){
				
		MatchedTagPool matchedTagMap = new MatchedTagPool();
		Iterator<Map.Entry<String, RetrivedPeptide>> it = this.entrySet().iterator();
		while(it.hasNext())
		{
			Map.Entry<String, RetrivedPeptide> entry = it.next();	
			Collections.sort( entry.getValue().getProteinList() );
			Peptide pept= new Peptide( entry.getKey(), entry.getValue().getProteinList() );
			
			if( shortTags.size() > 0 ){
				entry.getValue().addExtraTags(shortTags, Constants.maxModifiedMass, Constants.minModifiedMass);
			}
			
			for( TagMatchToPept mt : entry.getValue().getTagList() ){
				PeptideDBHit hit= new PeptideDBHit( pept, mt.staSite, mt.endSite, mt.nGap, mt.cGap );
				
				if( mt.ir == IonDirection.B_DIRECTION ){
					if(!(mt.matchedTag.isCtermOnly() && !hit.isCtermIon() || mt.matchedTag.isNtermOnly() && !hit.isNtermIon())){
						MatchedTag newTag = new MatchedTag(mt.matchedTag, hit, IonDirection.B_DIRECTION);
						matchedTagMap.add(newTag);						
					}
				}
				else if( mt.ir == IonDirection.Y_DIRECTION ){
					if(!(mt.matchedTag.isCtermOnly() && !hit.isCtermIon() || mt.matchedTag.isNtermOnly() && !hit.isNtermIon()))
					{
						MatchedTag newTag = new MatchedTag(mt.matchedTag, hit, IonDirection.Y_DIRECTION);
						matchedTagMap.add(newTag);
					}
				}
			}
		}
		return matchedTagMap;
	}
}










