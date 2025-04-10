package modi;

import java.util.Comparator;

public class PTM implements Comparable<PTM> {
	private int					id;
	private String				name;
	private String				fullName;
	private double				massDifference;
	private double				avgMassDifference;
	private AminoAcid			residue;
	private PTMPosition			position;
	private boolean				residueNTerm;	// residue must be N-term
	private boolean				residueCTerm;	// residue must be C-term
	private String				classification;	
	private int					modCount;
	private double				penalty;
	private double				diagnosticIon;
	private double				neutralLoss;
	
	public PTM( int id, String name, String fullName, double massDifference, double avgMassDifference, AminoAcid residue, PTMPosition position, String classification ){
		this.id						= id;
		this.name					= name;
		this.fullName				= fullName;
		this.massDifference			= massDifference;
		this.avgMassDifference		= avgMassDifference;
		this.residue				= residue;
		this.position				= position;
		this.residueNTerm			= false;
		this.residueCTerm			= false;
		this.classification 		= classification;
		this.modCount 				= 1;
		this.penalty 				= 1;
		this.diagnosticIon			= -1;
		this.neutralLoss			= 0;
		
		if( position==PTMPosition.ANY_N_TERM || position==PTMPosition.PROTEIN_N_TERM ) this.residueNTerm = true;
		if( position==PTMPosition.ANY_C_TERM || position==PTMPosition.PROTEIN_C_TERM ) this.residueCTerm = true;
	}

	public PTM( int id, String name, String fullName, double massDifference, double avgMassDifference, AminoAcid residue, PTMPosition position ) {
		this(id, name, fullName, massDifference, avgMassDifference, residue, position, "");
	}
	
	public PTM( int id, String name, String fullName, double massDifference, double avgMassDifference, AminoAcid residue, PTMPosition position, int mcount ) {
		this(id, name, fullName, massDifference, avgMassDifference, residue, position, "");
		this.modCount = mcount;
		if( mcount == 0 ) this.penalty = 0;
	}
	
	public PTM( double massDifference ) {
		this(-1, "", "", massDifference, massDifference, null, PTMPosition.ANYWHERE, "");
	}
	
	public int			getID()						{ return id; }
	public String		getName()					{ return name; }
	public String		getFullName()				{ return fullName; }
	public AminoAcid	getResidue()				{ return residue; }
	public double		getMassDifference()			{ return massDifference; }
	public double		getAvgMassDifference()		{ return avgMassDifference; }
	public PTMPosition	getPTMPosition()			{ return position; }
	public boolean		isResidueNTerm()			{ return residueNTerm; }
	public boolean		isResidueCTerm()			{ return residueCTerm; }
	public String		getClassification()			{ return classification; }
	
	public int			getModCount()				{ return modCount; }
	public double		getPenalty()				{ return penalty; }
	public double		getDiagnosticIon()			{ return diagnosticIon; }
	public double		getNeutralLoss()			{ return neutralLoss; }
	public void			setNeutralLoss(double d)	{ neutralLoss = d; }
	public void			setDiagnosticIon(double d)	{ diagnosticIon = d; }
	public void			setPenalty(double d)		{ penalty = d; }
	public void			setModCount(int d)		    { modCount = d; }
	
	public String	getSite(){ 
		if( residue != null ) return residue.toString();
		else {
			if (residueNTerm) return "N-term";
			else if (residueCTerm) return "C-term";
			else return "Xnull";
		}
	}
	
	public char	getAbbAA(){ 
		if( residue != null ) return residue.getResidue();
		else {
			if (residueNTerm) return 'X';
			else return 'Z';
		}
	}

	public String toString()
	{
		return this.name;
	}

	
	public String getDescription(){
		
		StringBuffer pos= new StringBuffer();
		
		if( this.position == PTMPosition.ANY_N_TERM )
			pos.append("N-term");
		else if( this.position == PTMPosition.ANY_C_TERM )
			pos.append("C-term");
		else if( this.position == PTMPosition.PROTEIN_N_TERM )
			pos.append("Protein N-term");
		else if( this.position == PTMPosition.PROTEIN_C_TERM )
			pos.append("Protein C-term");
		
		if( residue != null ){
			if( pos.length() != 0 ) pos.append(' ');
			pos.append( residue.toString() );	
		}
		
		pos.append( String.format("|%.4f|", massDifference) );
		
		return name+"|"+pos.toString();
	}//*/

	public int compareTo(PTM p){	// default comparator : mass
		if(this.id > p.id ) return 1;	
		else if(this.id < p.id ) return -1;	
		return 0;
	}
	
	public int hashCode(){ 
		return name.hashCode();
	}
	public boolean equals(Object o){		
		if(!(o instanceof PTM))
			   return false;		
		PTM p= (PTM)o;
		if( this.name.equals(p.name) && this.penalty==p.penalty && this.modCount == p.modCount )
			return true;
		return false;
	}
}

class PTMPosComparator implements Comparator<PTM>{
	public int compare( PTM x1, PTM x2 ){
		if( x1.getPTMPosition().ordinal() > x2.getPTMPosition().ordinal() )
			return 1;		
		else if( x1.getPTMPosition().ordinal() == x2.getPTMPosition().ordinal() )
			return 0;
		else return -1;
	}	
	public boolean equals(PTM x1, PTM x2){
		return x1.getPTMPosition().ordinal() == x2.getPTMPosition().ordinal();
	}
}









