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

	// Multi-mod: List<TagPeptide> 버전
	public CandidateContainer(List<TagPeptide> hmap, TagTrie trie) { // multi-mod
		ArrayList<TagPeptide> list = new ArrayList<TagPeptide>(hmap);
		list.sort(new TagPeptComparator());

		size = 0;
		modlist = new ChainTagPeptide[list.size()];

		int i = 0, n = list.size();
		while (i < n) {
			TagPeptide parent = list.get(i++);
			ChainTagPeptide ctp = new ChainTagPeptide(parent.pStart, parent.pEnd, parent.mTag);
			ctp.setConservedRegion(parent.pLeft, parent.pRight);

			while (i < n) {
				TagPeptide entry = list.get(i);
				if (!ctp.extend(parent, entry, trie)) {
					ctp.arrangeTags();
					modlist[size++] = ctp;
					ctp = new ChainTagPeptide(entry.pStart, entry.pEnd, entry.mTag);
					ctp.setConservedRegion(entry.pLeft, entry.pRight);
				}
				parent = entry;
				i++;
			}
			ctp.arrangeTags();
			modlist[size++] = ctp;
		}
	}

	// One-mod: List<MODPeptide> 버전
	public CandidateContainer(List<MODPeptide> hmap) { // one-mod
		ArrayList<MODPeptide> list = new ArrayList<MODPeptide>(hmap);
		Collections.sort(list);

		size = 0;
		modlist = new MODPeptide[list.size()];

		int i = 0, n = list.size();
		while (i < n) {
			MODPeptide parent = list.get(i++);
			while (i < n) {
				MODPeptide entry = list.get(i);
				if (!parent.extend(entry)) {
					modlist[size++] = parent;
					parent = entry;
				}
				i++;
			}
			modlist[size++] = parent;
		}
	}

}
