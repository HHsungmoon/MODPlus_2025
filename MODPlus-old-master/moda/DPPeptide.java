package moda;

import java.util.Comparator;

import msutil.IonGraph;
import msutil.PGraph;
import msutil.Scoring;

import modi.Constants;

public class DPPeptide implements Comparable<DPPeptide> {
	
	int score; 
	double prob;
	double peptMass;
	String peptide="X.XXX.X";
	double[] ptms;
	int stem;
	public int maxConsecutiveAALength;
	int protein;
	int NTT=0, noMods=0;
	
	public DPPeptide(){}
	public DPPeptide(String pept, int score, double[] ptms, int protein){
		this.peptide= pept;
		this.score= score;
		this.ptms= ptms;
		this.protein= protein;
	}
	
	public int getScore() 			{ return score; }
	public int getProtein() 		{ return protein; }
	public String getPeptide() 		{ return peptide; }
	public int getStem() 			{ return stem; }
	
	public boolean isConfident() {		
		if( (NTT-noMods) < 1 ) return false;
		if( prob < MODaConst.minCutOfProb ) return false;
		if( maxConsecutiveAALength < MODaConst.minConsecutiveAALength ) return false;
		return true;
	}

	public void setProteinAndNTT(int p, int n) { 
		protein += p;//put base in construct, and plus real pos
		NTT = n;
	}
	
	public void evaluatePSM(PGraph pg) {
		IonGraph iG = Scoring.PeptideSpectrumMatch( peptide, ptms, pg );
		prob = iG.getProb();
		score = iG.getRankScore();
		peptMass = iG.getCalculatedMW();
		maxConsecutiveAALength = iG.getMaxConsecutiveAALength();
		noMods = iG.getModifiedResd();
	}
	

	public String toString(){
		StringBuffer x= new StringBuffer();
		for( int i=0; i<ptms.length; i++){
			x.append(peptide.charAt(i));
			if( ptms[i] != 0 ){
				if( ptms[i] > 0 ) x.append("+");
				x.append( Constants.round(ptms[i]) );
			}
		}
		x.append("_"+score);
		return x.toString();
	}

	public boolean isSame( DPPeptide x ){
		if( this.peptide.compareTo(x.peptide) != 0 ) return false;		
		for(int i=0; i<ptms.length; i++){
			if( Math.abs(ptms[i]-x.ptms[i]) > 0.001 ) {
				return false;
			}
		}
		if( this.NTT != x.NTT ) return false;
		return true;
	}
		
	public int compareTo( DPPeptide x ) {			
		if( this.score < x.score ) return 1;
		else if( this.score > x.score ) return -1;
		
		if( this.score == 0 && x.score == 0 ) return -1;

		int thisLen = this.peptide.length(), xLen = x.peptide.length();	
		int thisMod = 0, xMod = 0;
		for(int i=0; i<thisLen; i++){
			if( this.ptms[i]!= 0 ) thisMod += Math.abs( this.ptms[i] );
		}
		for(int i=0; i<xLen; i++){
			if( x.ptms[i]!= 0 ) xMod += Math.abs( x.ptms[i] );
		}		
		if( thisMod > xMod ) return 1;
		else if( thisMod < xMod ) return -1;

		if( thisLen > xLen ) return 1;
		else if( thisLen < xLen ) return -1;
		else return 0;
	}	
}

class DPPeptideRefinedComparator implements Comparator<DPPeptide>{
	
	public int compare(DPPeptide AA, DPPeptide BB){
	
		if( AA.score < BB.score ) return 1;
		else if( AA.score > BB.score ) return -1;

		int aaMod = 0, bbMod = 0;
		int aaLen = AA.peptide.length(), bbLen = BB.peptide.length();
		
		for(int i=0; i<aaLen; i++){
			if( AA.ptms[i]!= 0 ) aaMod += Math.abs( AA.ptms[i] );
		}
		for(int i=0; i<bbLen; i++){
			if( BB.ptms[i]!= 0 ) bbMod += Math.abs( BB.ptms[i] );
		}
		
		if( aaMod > bbMod ) return 1;
		else if( aaMod < bbMod ) return -1;
		
		if( AA.prob < BB.prob ) return 1;
		else if( AA.prob > BB.prob ) return -1;
		
		if( aaLen > bbLen ) return 1;
		else if( aaLen < bbLen ) return -1;
		
		return 0;
	}	
	
	public boolean equals(DPPeptide AA, DPPeptide BB){
		return AA == BB;
	}
}























