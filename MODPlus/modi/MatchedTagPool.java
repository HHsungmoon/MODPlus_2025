package modi;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map;
import java.util.TreeMap;

public class MatchedTagPool extends TreeMap<Peptide, LinkedList<MatchedTag>>{
	public int add(MatchedTag tag)
	{
		if( !tag.checkCompatibility() )
			return 0;
		// New peptide tag

		Peptide matchedPeptide = tag.getMatchedPeptide();
		LinkedList<MatchedTag> tagList = this.get( matchedPeptide );
		if( tagList == null )
		{
			tagList = new LinkedList<MatchedTag>();
			tagList.add(tag);
			put( tag.getMatchedPeptide(), tagList );
			return 1;
		}

		tagList.add(tag);
		return tagList.size();
	}

	
	public String toString() 
	{
		StringBuffer output = new StringBuffer();
		output.append("Matched Tag Pool");
		Iterator<Map.Entry<Peptide, LinkedList<MatchedTag>>> it = this.entrySet().iterator();
		while(it.hasNext())
		{
			Map.Entry<Peptide, LinkedList<MatchedTag>> entry = it.next();
			output.append(entry.getKey()).append("\n");
			output.append(entry.getValue()).append("\n");
		}
		return output.toString();
	}
	
}
