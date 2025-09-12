package modi;

import java.util.ArrayList;
import java.util.Collections;

public class Peptide extends Sequence {
		
	private ArrayList<SrcProteinInfo> srcProteinInfo = new ArrayList<SrcProteinInfo>();
	
	public Peptide( int srcProteinID, int start, int end ) {
		srcProteinInfo.add( new SrcProteinInfo( srcProteinID, start, end ));
	}
	
	public Peptide( String sequence, ArrayList<SrcProteinInfo> src) {
		this.addAll( Sequence.getSequence(sequence) );
		this.srcProteinInfo = src;
	}
	
	public Peptide( Peptide peptide ) {
		this.addAll( peptide );
		this.srcProteinInfo = peptide.srcProteinInfo;
	}
	
	public ArrayList<SrcProteinInfo> getSrcProteinInfo() {
		return srcProteinInfo;
	}
	public int getNTT(){
		return srcProteinInfo.get(0).getNTT();
	}
	public int getStartCoordinate(){
		return srcProteinInfo.get(0).getStartPos();
	}
	public int getEndCoordinate(){
		return srcProteinInfo.get(0).getEndPos();
	}

	public double getMolecularMass()
	{
		return this.getMonoMass() + Constants.H2O;
	}

	public	int compareTo( Peptide p ) {
		for (int i=0; i<size() && i<p.size(); i++)
		{
			int r = get(i).compareTo(p.get(i));
			if (r<0) return -1;
			else if (r>0) return 1;
		}
		int r = size()-p.size();
		if (r<0) return -1;
		else if (r>0) return 1;
		else return 0;		
	}
	
	public int hashCode(){ 
		return srcProteinInfo.hashCode();
	}
	
	public boolean equals(Object o){		
		if(!(o instanceof Peptide))
			   return false;		
		Peptide p= (Peptide)o;
		if( this.srcProteinInfo.hashCode() == p.srcProteinInfo.hashCode() )
			return true;
		else return false;
	}
	
}
