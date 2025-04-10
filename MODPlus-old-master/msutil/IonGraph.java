package msutil;

import java.util.ArrayList;

import modi.Constants;
import modi.PTM;

public abstract class IonGraph extends ArrayList<IonNode> {
	
	static int NTermType = 0, CTermType = 1;
	static int SinglyCharged = 1, DoublyCharged = 2, TriplyCharged = 3;
	static double supportingFactor = 0.9;
	
	String 	sequence;
	double[] ptmMass;
	PTM[] 	 ptmList;
	int charge;
	
	double 	calculatedMW = 0;	
	double 	modifiedResd = 0;
	
	int 	rankScore = 0;
	int 	coverScore =0;
	double	prob = 0;
	double	massError = 0;
	
	double ionCoverage = 0;
	double seqCoverage = 0;
	
	int	   maxConsecutive = 0; 	
	int	   maxAAs = 0; //AA Length, not inos

	public double getCalculatedMW(){ return calculatedMW; }
	public double getMassError(){ return massError; }
	public int getRankScore(){ return rankScore; }
	public int getModifiedResd(){ return Constants.round(modifiedResd); }
	public double getProb(){ return prob; }

	public int getMaxConsecutiveAALength(){ return maxAAs; }
	
	public IonGraph( String peptide, double[] ptms, PGraph graph ){	
		sequence = peptide;
		ptmMass  = ptms;
		charge   = graph.getCharge();
		ptmList = null;
		construct();
	}
	
	public IonGraph( String peptide, double[] ptms, PTM[] pList, PGraph graph ){	
		sequence = peptide;
		ptmMass = ptms;
		charge = graph.getCharge();
		ptmList = pList;
		construct();
	}
	
	protected abstract void construct();
	public abstract void setScore( PGraph graph );
	protected abstract void setSecondaryScore( PGraph graph );
	public abstract void evaluateMatchQuality( PGraph graph );
	
	public void setMatchCoverage(){

		int len = sequence.length()-1;
		int[][] ionHits = new int[2][len];
		int[] seqHits = new int[len];
		
		for( IonNode node : this ){
			ionHits[node.type][node.index]++;
			seqHits[node.index]++;
		}
		
		int bbHit=0, yyHit=0, seqHit=0;
		for( int i=0 ; i<len; i++ ){
			if( ionHits[NTermType][i] != 0 ) bbHit++; 
			if( ionHits[CTermType][i] != 0 ) yyHit++;
			if( seqHits[i] != 0 ) seqHit++;
		}
		if( bbHit < yyHit ){
			bbHit = yyHit;
		}
		
		int maxConsec = 0, localsec = 0;
		for( int i=0 ; i<len; i++ ){// BION
			if( ionHits[NTermType][i] != 0 ) localsec++;
			else{
				if( maxConsec < localsec ) maxConsec = localsec;
				localsec = 0;		
			}
		}
		if( maxConsec < localsec ) maxConsec = localsec;
		
		localsec = 0;
		for( int i=0 ; i<len; i++ ){// YION
			if( ionHits[CTermType][i] != 0 ) localsec++;
			else{
				if( maxConsec < localsec ) maxConsec = localsec;
				localsec = 0;
			}
		}
		if( maxConsec < localsec ) maxConsec = localsec;
		
		maxConsecutive = maxConsec;
		maxAAs  = getMaxAAs(ionHits);
		
		double modifiedMass = 0;
		for(int i=0; i<ptmMass.length ; i++){
			if( ptmMass[i] != 0 ){
				modifiedMass += ptmMass[i];
			}
		}
		
		int modAcid = ( modifiedMass < Constants.precursorTolerance )?  0 : (int)Math.ceil( modifiedMass / 110 );	
		double penalty = modAcid + modifiedResd/2;
		seqCoverage = (double)seqHit/(len+penalty);		
		ionCoverage = (double)bbHit /(len+penalty);
	}
	
	protected int getMaxAAs(int[][] ionHits){

		int whichIonIntense = 2; //0:bion, 1:yion, 2:dont know	

		int len = ionHits[0].length;
		int maxConsec = 0, localsec = 0;		
		for( int i=0 ; i<len; i++ ){// YION
			if( ionHits[CTermType][i] != 0 ) localsec++;
			else{
				if( maxConsec < localsec ) maxConsec = localsec;
				localsec = -1;
			}
		}
		localsec++;
		if( maxConsec < localsec ) maxConsec = localsec;
		if( whichIonIntense == CTermType ) return maxConsec;

		if( whichIonIntense == NTermType ) maxConsec = 0;
		localsec = 0;	
		for( int i=0 ; i<len; i++ ){// BION
			if( ionHits[NTermType][i] != 0 ) localsec++;
			else{
				if( maxConsec < localsec ) maxConsec = localsec;
				localsec = -1;		
			}
		}
		localsec++;
		if( maxConsec < localsec ) maxConsec = localsec;		

		return maxConsec;
	}
	
}













