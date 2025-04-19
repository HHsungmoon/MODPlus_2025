package processedDB;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

public class CandidateContainer {
	
	int size = 0;
	MODPeptide[] modlist;
	
	public int size() { return size; }
	public MODPeptide[] getList() { return modlist; }

	/*
	public CandidateContainer( LinkedList<TagPeptide> hmap, TagTrie trie ){ //multi-mod
		// 질문: 정렬 필요한 이유
		Collections.sort( hmap, new TagPeptComparator() );

		size  = 0;
		modlist = new ChainTagPeptide[hmap.size()];

		while( hmap.size() != 0 ){
			TagPeptide parent = hmap.getFirst();
			ChainTagPeptide ctp = new ChainTagPeptide(parent.pStart, parent.pEnd, parent.mTag);
			ctp.setConservedRegion(parent.pLeft, parent.pRight);
			hmap.removeFirst();
			while( hmap.size() != 0 ){
				TagPeptide entry = hmap.getFirst();
				if( !ctp.extend(parent, entry, trie) ){
					ctp.arrangeTags();
					modlist[size++] = ctp;
					ctp = new ChainTagPeptide(entry.pStart, entry.pEnd, entry.mTag);
					ctp.setConservedRegion(entry.pLeft, entry.pRight);
				}
				parent = entry;
				hmap.removeFirst();
			}
			ctp.arrangeTags();
			modlist[size++] = ctp;
		}
	}
	*/
	  //이게 위에거보다 조금 더 빠르긴 한데... 정확한지 확인이 아직
	public CandidateContainer(List<TagPeptide> hmap, TagTrie trie) {
		if (hmap.isEmpty()) {
			size = 0;
			modlist = new ChainTagPeptide[0];
			return;
		}

		hmap.sort(new TagPeptComparator());

		size = 0;
		modlist = new ChainTagPeptide[hmap.size()];

		TagPeptide parent = hmap.get(0);
		ChainTagPeptide ctp = new ChainTagPeptide(parent.pStart, parent.pEnd, parent.mTag);
		ctp.setConservedRegion(parent.pLeft, parent.pRight);

		for (int i = 1; i < hmap.size(); i++) {
			TagPeptide entry = hmap.get(i);
			if (!ctp.extend(parent, entry, trie)) {
				ctp.arrangeTags();
				modlist[size++] = ctp;

				ctp = new ChainTagPeptide(entry.pStart, entry.pEnd, entry.mTag);
				ctp.setConservedRegion(entry.pLeft, entry.pRight);
			}
			parent = entry;
		}

		ctp.arrangeTags();
		modlist[size++] = ctp;
	}

	public CandidateContainer(List<MODPeptide> hmap) {
		Collections.sort(hmap);  // Comparable 구현되어 있어야 함

		size = 0;
		modlist = new MODPeptide[hmap.size()];

		if (hmap.isEmpty()) return;

		MODPeptide parent = hmap.get(0);
		for (int i = 1; i < hmap.size(); i++) {
			MODPeptide entry = hmap.get(i);
			if (!parent.extend(entry)) {
				modlist[size++] = parent;
				parent = entry;
			}
		}
		modlist[size++] = parent;
	}


	/*
	public CandidateContainer( List<MODPeptide> hmap ){ //one-mod
		Collections.sort( hmap );
	
		size  = 0;
		modlist = new MODPeptide[hmap.size()];
		
		while( hmap.size() != 0 ){
			MODPeptide parent = hmap.getFirst();
			hmap.removeFirst();
			while( hmap.size() != 0 ){
				MODPeptide entry = hmap.getFirst();			
				if( !parent.extend(entry) ){
					modlist[size++] = parent;				
					parent = entry;
				}			
				hmap.removeFirst();
			}	
			modlist[size++] = parent;
		}		
	}
	*/

}
