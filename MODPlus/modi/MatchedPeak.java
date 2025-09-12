package modi;

public class MatchedPeak extends Peak {
	AminoAcid NTermResidue;
	AminoAcid CTermResidue;
	
	public MatchedPeak(int index, double mass, double intensity, int charge, PeakProperty property, 
			AminoAcid NTermResidue, AminoAcid CTermResidue)
	{
		super(index, mass, intensity, charge, property);
		this.NTermResidue = NTermResidue;
		this.CTermResidue = CTermResidue;
	}

	
}
