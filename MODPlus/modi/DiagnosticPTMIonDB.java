package modi;

import java.util.ArrayList;

public class DiagnosticPTMIonDB extends ArrayList<DiagnosticPTMIon> {

	public	DiagnosticPTMIonDB()
	{//VEMS Table
		add(new DiagnosticPTMIon("Methyl", AminoAcid.getAminoAcid('R'), 143.1291, 0.212));
		add(new DiagnosticPTMIon("Methyl", AminoAcid.getAminoAcid('R'), 115.0866, 0.124));
		add(new DiagnosticPTMIon("Methyl", AminoAcid.getAminoAcid('R'), 112.0869, 0.126));
		add(new DiagnosticPTMIon("Methyl", AminoAcid.getAminoAcid('R'), 74.0713, 0.002));
		add(new DiagnosticPTMIon("Methyl", AminoAcid.getAminoAcid('R'), 70.0651, 0.002));
		add(new DiagnosticPTMIon("di-Methylation", AminoAcid.getAminoAcid('R'), 157.1448, 0.206));
		add(new DiagnosticPTMIon("di-Methylation", AminoAcid.getAminoAcid('R'), 115.0866, 0.124));
		add(new DiagnosticPTMIon("di-Methylation", AminoAcid.getAminoAcid('R'), 112.0869, 0.124));
		add(new DiagnosticPTMIon("di-Methylation", AminoAcid.getAminoAcid('R'), 88.0869, 0.01));
		add(new DiagnosticPTMIon("di-Methylation", AminoAcid.getAminoAcid('R'), 71.0604, 0));
		add(new DiagnosticPTMIon("tri-Methylation", AminoAcid.getAminoAcid('R'), 171.1604, 0.243));
		add(new DiagnosticPTMIon("Methyl", AminoAcid.getAminoAcid('K'), 115.123, 0.118));
		add(new DiagnosticPTMIon("Methyl", AminoAcid.getAminoAcid('K'), 98.0964, 0.008));
		add(new DiagnosticPTMIon("Methyl", AminoAcid.getAminoAcid('K'), 84.0808, 0.089));
		add(new DiagnosticPTMIon("di-Methylation", AminoAcid.getAminoAcid('K'), 129.1386, 0.555));
		add(new DiagnosticPTMIon("di-Methylation", AminoAcid.getAminoAcid('K'), 84.0808, 0.089));
		add(new DiagnosticPTMIon("tri-Methylation", AminoAcid.getAminoAcid('R'), 143.1543, 0.206));
		add(new DiagnosticPTMIon("tri-Methylation", AminoAcid.getAminoAcid('R'), 84.0808, 0.089));
		add(new DiagnosticPTMIon("Acetyl", AminoAcid.getAminoAcid('K'), 143.1179, 0.212));
		add(new DiagnosticPTMIon("Acetyl", AminoAcid.getAminoAcid('K'), 126.0913, 0.02));
		add(new DiagnosticPTMIon("Acetyl", AminoAcid.getAminoAcid('K'), 84.0808, 0.089));
		add(new DiagnosticPTMIon("Phospho", AminoAcid.getAminoAcid('Y'), 216.042, 0.142));
		add(new DiagnosticPTMIon("Oxidation", AminoAcid.getAminoAcid('M'), 120.0478, 0.327));
	}
	
	public	ArrayList<DiagnosticPTMIon> getDiagIon(PTM ptm)
	{
		ArrayList<DiagnosticPTMIon> result = new ArrayList<DiagnosticPTMIon>();
		for(DiagnosticPTMIon ion : this)
		{
			if(ion.name.equalsIgnoreCase(ptm.getName()) && ion.residue == ptm.getResidue())
				result.add(ion);
		}
		return result;
	}
}
