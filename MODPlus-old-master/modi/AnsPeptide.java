package modi;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;

import msutil.*;
import processedDB.PeptideMatchToProtein;

public class AnsPeptide implements Comparable<AnsPeptide> {

	static final double  pointRounding = 1000;
	final Peptide 	peptide;
	final String 		ptmComb;
	final double[] 	ptms;
	final PTM[] 		ptmList;
	double 		pepMass;

	int 		score;
	double 		prob;
	double 		mzDev = 0;
	double 		roundprob=0;

	public AnsPeptide(Peptide p, String comb, double[] ptm, int s){
		peptide = p;
		ptmComb = comb.trim();
		ptms = ptm;
		ptmList = null;
		score = s;
	}

	public AnsPeptide(Peptide p, String comb, double[] ptm, PTM[] pList, int s){
		peptide = p;
		ptmComb = comb.trim();
		ptms = ptm;
		ptmList = pList;
		score =s;
	}

	public String getPeptideSequence(){ return peptide.toString(); }

	int round(double a) {
		if (a > 0) return (int) (a + 0.5);
		else return (int) (a - 0.5);
	}

	IonGraph PeptideSpectrumMatch(String peptide, double[] ptms, PTM[] ptmList, PGraph graph) {//for modeye final scoring
		IonGraph iGraph;
		if (Constants.INSTRUMENT_TYPE == Constants.msms_type.QTOF) iGraph = new TOFGraph(peptide, ptms, ptmList, graph);
		else iGraph = new TRAPGraph(peptide, ptms, ptmList, graph);

		iGraph.evaluateMatchQuality(graph);
		return iGraph;
	}

	public void evaluatePSM(PGraph pg){
		IonGraph iG = PeptideSpectrumMatch(peptide.toString(), ptms, ptmList, pg );
		score = iG.getRankScore();
		prob 	= iG.getProb();
		roundprob 	= round(iG.getProb()*pointRounding)/pointRounding;

		mzDev 	= iG.getMassError();
		pepMass = iG.getCalculatedMW();
		score -= (2-peptide.getNTT())*3;
	}

	public String toMODPlus( double obMW, ArrayList<PeptideMatchToProtein> protMatch ){
		//For MOD Plus, OrgDB & DynamicDB
		StringBuffer x= new StringBuffer();
		x.append( String.format("%.4f\t", pepMass) );
		x.append( String.format("%.4f\t", obMW-pepMass) );

		x.append( String.format("%d\t", score) );
		x.append( String.format("%.4f\t", prob) );

		String peptStr = peptide.toString();
		String wrapAA = protMatch.getFirst().getWrappingAA();
		x.append( wrapAA.charAt(0)+".");
		for( int i=0; i<ptms.length; i++){
			x.append(peptStr.charAt(i));
			if( ptms[i] != 0 ){
				if( ptms[i] > 0 ) x.append("+");
				x.append( String.format("%.3f", ptms[i]) );
			}
		}
		x.append( "."+wrapAA.charAt(1)+"\t" );
		x.append( protMatch.getFirst().toString() );

		for( int i=1; i<protMatch.size(); i++ ){
			x.append( ";" + protMatch.get(i).toString() );
		}
		x.append( "\t"+ptmComb );

		return x.toString();
	}

	public int compareTo( AnsPeptide x ){

		if( this.roundprob < x.roundprob ) return 1;
		else if( this.roundprob > x.roundprob ) return -1;

		if( this.score < x.score ) return 1;
		else if( this.score > x.score ) return -1;

		double mePtm = getNumOfComb(), youPtm = x.getNumOfComb();
		if( mePtm > youPtm ) return 1;
		else if( mePtm < youPtm ) return -1;

		if( mePtm > 0 ){
			int nt = isNtermModified(), xnt = x.isNtermModified();
			if( nt < xnt ) return 1;
			else if( nt > xnt ) return -1;
		}

		if( mzDev > x.mzDev ) return 1;
		else if( mzDev < x.mzDev ) return -1;

		if( peptide.getNTT() < x.peptide.getNTT() ) return 1;
		else if( peptide.getNTT() > x.peptide.getNTT() ) return -1;

		return 0;

	}

	public int hashCode(){
		return peptide.hashCode()+ptmComb.hashCode();
	}
	public boolean equals(Object o){
		if(!(o instanceof AnsPeptide p))
			return false;
		return this.peptide.equals(p.peptide) && this.ptmComb.equals(p.ptmComb);
	}

	public double getNumOfComb(){
		double numComb = 0;
		HashSet<String> t = new HashSet<>();

		for(int i=0; i<ptmList.length; i++){
			if( ptmList[i] != null ) {
				//	if( ptmList[i].getPenalty() != 0)
				{
					numComb += ptmList[i].getPenalty();
					t.add( ptmList[i].getName() );
				}
			}
		}

		numComb *= t.size();
		return numComb;
	}

	public int isCtermModified(){ return ( ptms[ptms.length-1] != 0 )? 1 : 0; }
	public int isNtermModified(){
		if( ptms[0]==0 ) return 0;
		else{
			if( ptmList != null ) return ptmList[0].getModCount();
			else return 1;
		}
	}

	public int getModifiedDelta(){
		double modifiedDelta = 0;
		for( int i =0; i<ptms.length; i++){
			if( ptms[i] != 0 )
				modifiedDelta += Math.abs(ptms[i]);
		}
		return (int)Math.round(modifiedDelta);
	}
}