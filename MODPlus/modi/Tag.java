package modi;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;

public class Tag extends ArrayList<Peak> {
	Sequence	tagSequence = new Sequence();
	Spectrum	sourceSpectrum;
	boolean		bIonOnly = false; 
	boolean		yIonOnly = false;
	boolean		nTermOnly = false;
	boolean		cTermOnly = false;
	private double 		offset= 0;
	double 		tagScore= 0.;
	
	public	Tag()	{}
	public	Tag(Peak p1, Peak p2, AminoAcid residue, Spectrum sourceSpectrum) {
		this.add(p1);
		this.add(p2);
		
		tagSequence.add(residue);
		this.sourceSpectrum = sourceSpectrum;
	}
	
	public	Tag(Tag t)
	{
		this.addAll(t);
		tagSequence.addAll(t.sequence());
		this.tagScore= t.tagScore;
		sourceSpectrum = t.getSourceSpectrum();
		bIonOnly = t.bIonOnly;
		yIonOnly = t.yIonOnly;
		nTermOnly = t.nTermOnly;
		cTermOnly = t.cTermOnly;
	}
	public	double		getOffset() { return offset; }	
	public	Sequence	sequence()	{ return tagSequence; }

	public	Tag			reverseTag()
	{
		Tag revTag = new Tag(this);
		Collections.reverse( revTag.tagSequence );
		return revTag;
	}

	public	void		setTagScore()	
	{ 
		double score = 0.;
		Iterator<Peak> it = this.iterator();
		while(it.hasNext())
		{
			score += it.next().getProbability();
		}
		tagScore= score;
	}
	
	public	double		getScore() { return tagScore; }		
	public	Spectrum	getSourceSpectrum() { return sourceSpectrum; }
	// t1 | t2 extendable : 1,  elas : 0
	public	static	boolean		extendable(Tag t1, Tag t2)
	{
		boolean ex= true;
		if(t1.get(t1.size()-1) != t2.get(0))
			ex = false;			
		return ex;
	}
	public	boolean	contains(Tag t)
	{
		return this.sourceSpectrum == t.sourceSpectrum && this.containsAll(t) ;
	}

	public	boolean	isNtermOnly() 	{ return nTermOnly; }
	public	boolean	isCtermOnly()	{ return cTermOnly; }
	
	// last peak of this tag and first peak of t are same
	public	static Tag merge(Tag t1, Tag t2)
	{
		assert(t1.getSourceSpectrum() == t2.getSourceSpectrum() && Tag.extendable(t1, t2));
		
		Tag merged = new Tag(t1);
		merged.remove( merged.size()-1 );
		merged.addAll(t2);
		merged.sequence().addAll(t2.sequence());
		
		return merged;
	}
	
	public	String toString()
	{		
		StringBuffer str = new StringBuffer(tagSequence.toString());
		for( Peak p : this ){
			str.append( String.format(" | %.2f", p.mass) );
		}
		str.append(" = " + this.tagScore);
		return str.toString();
	}
	
	public	double getBIonNtermOffset()
	{
		double mean = 0;
		for(int i=0; i<this.size(); i++){
			mean += this.get(i).getMass();
		}
		mean -= getMassCumSumFromNTerm();
		mean /= this.size();
		return mean - Constants.B_ION_OFFSET;//*/
	}
	
	public double getBIonCtermOffset()//(double motherMass)
	{
		double mean = 0;
		for(int i=0; i<this.size(); i++){
			mean += this.get(i).getMass();
		}
		mean += getMassCumSumFromCTerm();
		mean /= this.size();
		return sourceSpectrum.getCorrectedMW() - mean - Constants.UNIT_MASS*17;//*/
	}
	
	public double getYIonNtermOffset()//(double motherMass)
	{
		double mean = 0;
		for(int i=0; i<this.size(); i++){
			mean += this.get(i).getMass();
		}
		mean += getMassCumSumFromNTerm();
		mean /= this.size();
		return sourceSpectrum.getCorrectedMW() - mean + Constants.UNIT_MASS;//*/
	}
	
	public double getYIonCtermOffset()
	{
		double mean = 0;
		for(int i=0; i<this.size(); i++){
			mean += this.get(i).getMass();			
		}
		mean -= getMassCumSumFromCTerm();
		mean /= this.size();
		return mean - Constants.Y_ION_OFFSET;//*/
	}
	
	public double getMassCumSumFromNTerm(){
		double sum = 0;
		for(int i=0; i<this.size(); i++){
			sum += tagSequence.getMonoMass(0, i);
		}
		return sum;
	}
	public double getMassCumSumFromCTerm(){
		double sum = 0;
		for(int i=this.size()-1; i>-1; i--){
			sum += tagSequence.getMonoMass(i, this.size()-1);
		}
		return sum;
	}

	public int hashCode()
	{
		return tagSequence.toString().hashCode();
	}
	
	public boolean equals(Object o)
	{
		if (o instanceof Tag){
			Tag x= (Tag)o;
			if( this.sequence().toString().equals(x.sequence().toString()) )
				if ( Constants.fEqual(this.getBIonNtermOffset(), x.getBIonNtermOffset()) )
					if ( Constants.fEqual(this.getBIonCtermOffset(), x.getBIonCtermOffset()) )
						return true;
		}
		return false;
	}
}

class TagComparator implements Comparator<Tag> {
	public int compare(Tag t1, Tag t2) {
		if(t1.getScore() > t2.getScore())
			return 1;
		else if(t1.getScore() == t2.getScore())
			return 0;
		else
			return -1;
	}
}












