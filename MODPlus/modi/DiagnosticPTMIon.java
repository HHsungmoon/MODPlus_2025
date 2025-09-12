package modi;

public class DiagnosticPTMIon {
	String		name;
	AminoAcid	residue;
	double		mass;
	double		prob;
	
	public DiagnosticPTMIon(String name, AminoAcid residue, double mass, double prob)
	{
		this.name = name;
		this.residue = residue;
		this.mass = mass;
		this.prob = prob;
	}
	
	public	String		getName()		{ return name; }
	public	double		getMass()		{ return mass; }
	public	double		getProb()		{ return prob; }
}

