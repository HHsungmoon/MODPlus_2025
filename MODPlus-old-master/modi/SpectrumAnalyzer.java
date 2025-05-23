package modi;

import java.util.LinkedList;

import java.util.ListIterator;
import moda.ThreadPoolManager;
import msutil.PGraph;
import msutil.ProtCutter;

import processedDB.RetrivedPeptideMap;
import processedDB.TagTrie;

public class SpectrumAnalyzer {
	
	private	SpectrumAnalyzer() {}
	public	static TagPool buildTagPool( Spectrum sourceSpec ) {
		if( sourceSpec == null ) return null;

		for(int i=0; i<sourceSpec.size(); i++){
			if( Math.abs(sourceSpec.get(i).getMass() - sourceSpec.getPrecursor() ) < 2 ) {
				sourceSpec.remove(i);
				i--;			
				continue;
			}
			sourceSpec.get(i).setIndex(i);
		}//temporal
	
		sourceSpec.normalizeIntensityLocally();
		
		int extra = ( sourceSpec.getCharge() > 2 && Constants.INSTRUMENT_TYPE != Constants.msms_type.QTOF )? 2 : 0; 		
		sourceSpec.peakSelection(Constants.selectionWindowSize, Constants.minNumOfPeaksInWindow+extra );	
		TagPool primitiveTags = sourceSpec.generateTags(Constants.minTagLength, Constants.minTagLengthPeptideShouldContain, Constants.massToleranceForDenovo);
	
		return primitiveTags;
	}

	public	static TagChainPool buildTagChain(MatchedTagPool matchedTags)
	{
		TagChainPool tagChainPool = new TagChainPool();
		tagChainPool.buildTagChainPoolSafely(matchedTags);
		return tagChainPool;
	}

	public static boolean interpretTagChain(PTMDB ptmDB, TagChainPool tcPool, PGraph graph) {
		Spectrum sourceSpectrum = null;
		boolean specAnnotated = false;

		for (LinkedList<TagChain> tagChainList : tcPool.values()) {
			ListIterator<TagChain> listIt = tagChainList.listIterator();

			while (listIt.hasNext()) {
				TagChain tc = listIt.next();

				boolean allGapAnnotated = true;
				if (sourceSpectrum == null) {
					sourceSpectrum = tc.sourceSpectrum;
				}

				Peptide pep = tc.getMatchedPeptide();
				for (SpecInterpretation si : tc) {
					if (!(si instanceof Gap)) continue;

					Gap gap = (Gap) si;
					Sequence subSeq = pep.subSequence(gap.getStart(), gap.getEnd() + 1);
					PTMSearchResult interpretation = ptmDB.searchPTM(subSeq, gap.getOffset(), gap.getPosition());

					if (!interpretation.isInterpreted()) {
						gap.setInterpreted(false);
						allGapAnnotated = false;
						tc.setAllGapAnnotated(false);
						break; // PTM 검색 실패 → 이 TagChain 제거 대상
					} else {
						gap.setInterpreted(true);
						gap.setInterpretation(interpretation, graph);
					}
				}

				if (allGapAnnotated) {
					tc.setAllGapAnnotated(true);
					specAnnotated = true;
				} else {
					listIt.remove(); // 안전한 제거 방식
				}
			}
		}

		return specAnnotated;
	}


	public static MatchedTagPool extendedBuildMatchedTagPool(TagPool primitiveTags, double motherMass,
			TagTrie ixPDB, ProtCutter enzyme, int NTT)
	{
		if( primitiveTags == null || ixPDB == null )
			return null;
		int slotIdx = ThreadPoolManager.getSlotIndex();
		double minDelta = (Constants.minModifiedMass < 0)? Constants.minModifiedMass - Constants.gapTolerance[slotIdx] : - Constants.gapTolerance[slotIdx];
		double maxDelta = (Constants.maxModifiedMass > 0)? Constants.maxModifiedMass + Constants.gapTolerance[slotIdx] : + Constants.gapTolerance[slotIdx];
		TagPool longTags = primitiveTags.extractAbove(Constants.minTagLengthPeptideShouldContain);

		int realTag = 0;
		double orbMass= motherMass - Constants.H2O;		
		RetrivedPeptideMap searchResults= new RetrivedPeptideMap();
		for(Tag tag : longTags){
			
			RetrivedPeptideMap bRes= ixPDB.getRetrivedPeptides(orbMass, enzyme, NTT, tag.getBIonNtermOffset()-Constants.NTERM_FIX_MOD, tag, 
					tag.getBIonCtermOffset()-Constants.CTERM_FIX_MOD, IonDirection.B_DIRECTION, minDelta, maxDelta, Constants.gapTolerance[slotIdx]);
			searchResults.combine(bRes);
			
			Tag reverseTag= tag.reverseTag();
			RetrivedPeptideMap yRes= ixPDB.getRetrivedPeptides(orbMass, enzyme, NTT, reverseTag.getYIonNtermOffset()-Constants.NTERM_FIX_MOD, reverseTag, 
					reverseTag.getYIonCtermOffset()-Constants.CTERM_FIX_MOD, IonDirection.Y_DIRECTION, minDelta, maxDelta, Constants.gapTolerance[slotIdx]);
			searchResults.combine(yRes);
			realTag++;
			
			if( realTag > Constants.MAX_TAG_SIZE*2 ) break;
		}
		return searchResults.convertToMatchedTagPool( primitiveTags.extract(Constants.minTagLength, Constants.minTagLengthPeptideShouldContain) );
	}

}





