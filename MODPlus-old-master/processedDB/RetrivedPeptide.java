package processedDB;

import java.util.ArrayList;

import msutil.MSMass;

import modi.Constants;
import modi.IonDirection;
import modi.SrcProteinInfo;
import modi.Tag;
import modi.TagPool;

public class RetrivedPeptide implements Comparable<RetrivedPeptide> {
	
	private String 	peptide;
	private ArrayList<SrcProteinInfo> mProteinList = new ArrayList<SrcProteinInfo>();
	private ArrayList<TagMatchToPept> mTagList = new ArrayList<TagMatchToPept>();
	
	public RetrivedPeptide(String seq, int pro, int start, int end, int ntt){
		peptide= seq;
		mProteinList.add( new SrcProteinInfo( pro, start, end, ntt ) );
	}
	public String getSequence(){ return peptide; }
	public ArrayList<SrcProteinInfo> getProteinList(){ return mProteinList; }
	public ArrayList<TagMatchToPept> getTagList(){ return mTagList; }
	
	public void addTag(Tag tag, double n, double c, int site, int end, IonDirection ir){
		mTagList.add( new TagMatchToPept(tag, n, c, site, end, ir) );
	}	
	
	public void addProtein( SrcProteinInfo sc ){
		mProteinList.add( sc );
	}
	
	public void addTag( TagMatchToPept sc ){
		mTagList.add( sc );
	}
	
	public void addExtraTags( TagPool shortTags, double maxModified,  double minModified ){
		
		for(Tag tag : shortTags){	
			boolean included = false;
			for( TagMatchToPept ttp : mTagList ){
				if( ttp.matchedTag.containsAll(tag) ) {
					included = true;
					break;
				}
			}
			if( included ) {
				continue;
			}
			
			ArrayList<TagMatchToPept> addedTags = new ArrayList<TagMatchToPept>();
			
			int index = 0, site = -1;
			String tSeq = tag.sequence().toString();
			double ntFlank=tag.getBIonNtermOffset(), ctFlank=tag.getBIonCtermOffset();
			double preMass = Constants.NTERM_FIX_MOD;
			
			while( ( site= peptide.indexOf(tSeq, index) ) != -1 ){
				preMass += MSMass.getPepMass(peptide.substring(index, site));
				double nGap= ntFlank - preMass;
				if( ( nGap < maxModified ) && ( nGap > minModified ) ){
					double cGap= ctFlank - Constants.CTERM_FIX_MOD - MSMass.getPepMass(peptide.substring(site+2));
					if( ( cGap < maxModified ) && ( cGap > minModified ) ){
						addedTags.add(new TagMatchToPept(tag, nGap, cGap, site, site+1, IonDirection.B_DIRECTION));
						if( addedTags.size() > 2 ) break;
					}
				}
				index= site+1;
				preMass += MSMass.getAAMass(peptide.charAt(site));
			}
			if( addedTags.size() > 2 ) continue;
			
			Tag reverseTag= tag.reverseTag();
			index = 0; site = -1;
			tSeq = reverseTag.sequence().toString();
			ntFlank= reverseTag.getYIonNtermOffset();
			ctFlank= reverseTag.getYIonCtermOffset();
			preMass = Constants.NTERM_FIX_MOD;
			
			while( ( site= peptide.indexOf(tSeq, index) ) != -1 ){
				preMass += MSMass.getPepMass(peptide.substring(index, site));
				double cGap= ntFlank - preMass;
				if( ( cGap < maxModified ) && ( cGap > minModified ) ){
					double nGap= ctFlank - Constants.CTERM_FIX_MOD - MSMass.getPepMass(peptide.substring(site+2));
					if( ( nGap < maxModified ) && ( nGap > minModified ) ){
						addedTags.add(new TagMatchToPept(reverseTag, cGap, nGap, site, site+1, IonDirection.Y_DIRECTION));
						if( addedTags.size() > 2 ) break;
					}
				}
				index= site+1;
				preMass += MSMass.getAAMass(peptide.charAt(site));
			}
			
			if( addedTags.size() < 3 ){
				for( TagMatchToPept ttp: addedTags ){
					mTagList.add( ttp );
				}
			}
		}
	}
		
	public int compareTo(RetrivedPeptide rp) {	// default comparator : mass
		if( peptide.compareTo(rp.peptide) > 0 ) return 1;
		else if( peptide.compareTo(rp.peptide) < 0 ) return -1;
		else return 0;
	}	

}
