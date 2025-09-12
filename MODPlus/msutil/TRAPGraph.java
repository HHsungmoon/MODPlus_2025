package msutil;

import java.util.ArrayList;
import java.util.Collections;

import modi.Constants;
import modi.PTM;

public class TRAPGraph extends IonGraph {
		
	public TRAPGraph( String peptide, double[] ptms, PGraph graph ){	
		super(peptide, ptms, graph);
	}
	public TRAPGraph( String peptide, double[] ptms, PTM[] ptmList, PGraph graph ){	
		super(peptide, ptms, ptmList, graph);
	}
	
	protected void construct(){
		int len = sequence.length()-1;	
		int i;
		
		double ion_mz= Constants.B_ION_OFFSET + Constants.NTERM_FIX_MOD;					
		for(i=0; i<len; i++ ){
			ion_mz += MSMass.getAAMass( sequence.charAt(i) );
			if( ptmMass[i] != 0 ){
				ion_mz += ptmMass[i];
				if( ptmList == null ) modifiedResd += 1;
				else modifiedResd += ptmList[i].getPenalty();
			}
			this.add( new IonNode(i, NTermType, SinglyCharged, ion_mz, 3) );
			if( charge > 2 ) this.add( new IonNode(i, NTermType, DoublyCharged, (ion_mz+Constants.Proton)/2, 4) );
			if( charge > 3 ) this.add( new IonNode(i, NTermType, TriplyCharged, (ion_mz+Constants.Proton*2)/3, 6) );					
		}
		ion_mz += MSMass.getAAMass( sequence.charAt(i) );
		if( ptmMass[i] != 0 ){
			ion_mz += ptmMass[i];
			if( ptmList == null ) modifiedResd += 1;
			else modifiedResd += ptmList[i].getPenalty();
		}
		calculatedMW = ion_mz - Constants.B_ION_OFFSET + Constants.H2O + Constants.CTERM_FIX_MOD;
		
		ion_mz= Constants.Y_ION_OFFSET + Constants.CTERM_FIX_MOD;		
		for(i=len; i>0; i-- ){
			ion_mz += MSMass.getAAMass( sequence.charAt(i) );
			ion_mz += ptmMass[i];		
			this.add( new IonNode(i-1, CTermType, SinglyCharged, ion_mz, 2) );
			if( charge > 2 ) this.add( new IonNode(i-1, CTermType, DoublyCharged, (ion_mz+Constants.Proton)/2, 1) );
			if( charge > 3 ) this.add( new IonNode(i-1, CTermType, TriplyCharged, (ion_mz+Constants.Proton*2)/3, 5) );
		}
		Collections.sort(this);
	}
	
	public void setScore( PGraph graph ){ // iontype(n/c), charge(1/2)		
		graph.refresh();
		ArrayList<Integer> observedPaeks = new ArrayList<Integer>();		
		double score = 0;						
		for(int i=0; i<this.size(); i++){
			int index = graph.getIndexOfMatchedPeak( this.get(i).mz );
			if( index == -1 ){
				this.remove(i);
				i--;			
			}
			else{				
				score += graph.get(index).norm;
				this.get(i).observed = graph.get(index).mass;
				this.get(i).intensity = graph.get(index).norm;
				observedPaeks.add( index );
			}
		}

		int i = 0;
		for( IonNode node : this ){
			score += graph.getSupportingScore( observedPaeks.get(i++), node.type, node.charge, supportingFactor );			
		}
		
		rankScore = (int)(score - modifiedResd*Constants.rNorm[0]);
	}
	
	protected void setSecondaryScore( PGraph graph ){ // iontype(n/c), charge(1/2)		
		
		if( Constants.MSResolution == 0 ){		
			massError = Math.abs( graph.getCorrectedMW() - calculatedMW );
		}
		else
		{
			massError = Math.abs( graph.getObservedMW() - calculatedMW );
			int isoerr = (int)Math.round( massError / Constants.IsotopeSpace );
			massError = Math.abs( massError - isoerr*Constants.IsotopeSpace );
		}
		
		ArrayList<Integer> observedPaeks = new ArrayList<Integer>();
		graph.refresh();				
		double score = 0;						

		for(int i=0; i<this.size(); i++){
						
			int index = graph.getIndexOfMatchedPeak( this.get(i).mz );
			if( index == -1 ){
				this.remove(i);
				i--;			
			}
			else{				
				score += graph.get(index).norm;
				this.get(i).observed = graph.get(index).mass;
				this.get(i).intensity = graph.get(index).norm;				
				observedPaeks.add( index );
			}
		}

		int i = 0;
		for( IonNode node : this ){
			score += graph.getSupportingScore( observedPaeks.get(i++), node.type, node.charge, supportingFactor );	
		}	

		double EPS = score - modifiedResd*Constants.rNorm[0];
		double UnEPS = graph.getPenaltySCoreForUnassignedIons(this.size());

		rankScore = (int)Math.round(EPS);	
		coverScore = (int)Math.round(EPS-UnEPS);
	}
	
	public void evaluateMatchQuality(PGraph graph){ 
	
		setSecondaryScore( graph );
		setMatchCoverage();
	
		double discriminant = 0;		
		
		if( graph.getCharge() > 2 )//extend isotope
			discriminant= 0.1729*coverScore + 13.2683*ionCoverage + 1.8631*seqCoverage - 14.7483;		
		else
			discriminant= 0.1248*coverScore + 6.7925*ionCoverage + 5.5328*seqCoverage - 11.9468;//*/			
		
		prob = Scoring.getOddProbability( discriminant );
	}
	
}
