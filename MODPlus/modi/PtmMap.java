
package modi;

import java.util.ArrayList;

import msutil.IonGraph;
import msutil.PGraph;
import msutil.Scoring;

public class PtmMap {

	private int[][] ptmFreqTable;
	private ArrayList<PSM> history= new ArrayList<PSM>();
	
	private int lowestCol = (int)Math.round( Constants.minModifiedMass );
	private int highestCol = (int)Math.round( Constants.maxModifiedMass );
	private int columnSize= highestCol - lowestCol + 10;
	private int siteColumn= 23;
	
	public PtmMap(){			
		ptmFreqTable= new int[columnSize][siteColumn];		
	}

	public void hitPtmFreqTable( ScanCap spec, Peptide peptide, PTM ptm, PGraph pg ){
	
		double mass= ptm.getMassDifference(); 
		AminoAcid aa= ptm.getResidue();		
		String residue = ( ptm.getPTMPosition() == PTMPosition.ANY_N_TERM )? "nterm" : aa.toString();
		
		int delta = (int)Math.round(mass);
		if( delta < lowestCol || highestCol < delta ) return;
	
		ptmFreqTable[delta-lowestCol][aa.getIndex()+2]++;
		if( residue.compareToIgnoreCase("nterm") == 0 )
			ptmFreqTable[delta-lowestCol][0]++;
		
		double ptms[] = new double[peptide.size()];
		ptms[ptm.getID()] = delta;
		IonGraph iG = Scoring.PeptideSpectrumMatch(peptide.toString(), ptms, pg );
		
		history.add( new PSM(spec, peptide, ptm.getID(), delta, iG.getProb()) );
	}

	private class PSM implements Comparable<PSM> {
		ScanCap spec;
		Peptide pept;
		int site;
		int delta;
		double prob;
		public PSM( ScanCap s, Peptide p, int t, int dt, double score ){
			spec = s;
			pept = p;
			site = t;
			delta = dt;
			prob = score;
		}
		private int getModRes() { return pept.get(site).getIndex(); }

		public int compareTo(PSM p){
			
			if( this.getModRes() >  p.getModRes() ) return 1;
			else if( this.getModRes() <  p.getModRes() ) return -1;
			
			if( this.delta > p.delta ) return 1;
			else if( this.delta < p.delta ) return -1;
			else return 0;
		}
	}
	
}
