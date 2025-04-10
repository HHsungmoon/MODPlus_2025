package processedDB;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

public class HeatedDB extends HashMap<String, int[]> {

	public void add(String prot, int stem, int start, int end){
		int[] local = this.get(prot);			
		if ( local == null ) {
			local = new int[3]; //0: stem, 1: start site, 2: end site
			local[0] = stem;
			local[1] = start;
			local[2] = end;
			this.put( prot, local );
		}
		else {
			if( start < local[1] ) local[1] = start;
			if( local[2] < end ) local[2] = end;
		}
		
	}
	
	public TagTrie getPartialDB(StemTagTrie ixPDB) throws Exception {		
		int res = 0;
		ProxDB parts = new ProxDB();			
		Iterator<Map.Entry<String, int[]>> it = this.entrySet().iterator();
		while( it.hasNext() ) {
			Map.Entry<String, int[]> entry = it.next();			
			int[] local = entry.getValue();
			ProtDatabase ptdb = ixPDB.get(local[0]);
			String stripSequence = ptdb.getSequenceAroundPeptide(local[1], local[2], 3);
			parts.add( new Prox(entry.getKey(), stripSequence) );
			res += stripSequence.length();
		}
		parts.setSizeOfResidues(res);
		parts.setSizeOfEntries(parts.size());		
		return new TagTrie(parts);
	}

}
