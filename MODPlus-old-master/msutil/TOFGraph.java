package msutil;

import java.util.ArrayList;
import java.util.Collections;

import modi.Constants;
import modi.PTM;

public class TOFGraph extends IonGraph {
	
	public TOFGraph( String peptide, double[] ptms, PGraph graph ){	
		super(peptide, ptms, graph);
	}
	
	public TOFGraph( String peptide, double[] ptms, PTM[] ptmList, PGraph graph ){	
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
			this.add( new IonNode(i, NTermType, SinglyCharged, ion_mz, 2) );
			if( charge > 2 && i > 2 ) this.add( new IonNode(i, NTermType, DoublyCharged, (ion_mz+Constants.Proton)/2, 4) );
		}
		
		ion_mz += MSMass.getAAMass( sequence.charAt(i) );
		if( ptmMass[i] != 0 ){
			ion_mz += ptmMass[i];
			if( ptmList == null ) modifiedResd += 1;
			else modifiedResd += ptmList[i].getPenalty();
		}
		calculatedMW = ion_mz - Constants.B_ION_OFFSET + Constants.H2O + Constants.CTERM_FIX_MOD;
		
		ion_mz= Constants.Y_ION_OFFSET + Constants.CTERM_FIX_MOD;/////////////////////////////				
		for(i=len; i>0; i-- ){
			ion_mz += MSMass.getAAMass( sequence.charAt(i) );
			ion_mz += ptmMass[i];		
			this.add( new IonNode(i-1, CTermType, SinglyCharged, ion_mz, 1) );
			if( charge > 2 && i < len-1 ) this.add( new IonNode(i-1, CTermType, DoublyCharged, (ion_mz+Constants.Proton)/2, 3) );
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
			if( node.type == NTermType && node.index < 2 ) score += graph.getNormOfMatchedPeak( node.mz-Constants.A_ION_OFFSET );	
			score += graph.getSupportingScore( observedPaeks.get(i++), node.type, node.charge, supportingFactor );
		}
		
		double immSupprot = 0;
		for(int imm=0; imm<sequence.length(); imm++ ){
			immSupprot += graph.getNormOfMatchedPeak( MSMass.getAAMass(sequence.charAt(imm))+ptmMass[imm] + Constants.IMM_OFFSET );
		}
		score += immSupprot*supportingFactor;
		
		rankScore = (int)(score - modifiedResd*Constants.rNorm[0]);
	}
	
	protected void setSecondaryScore( PGraph graph ){ // iontype(n/c), charge(1/2), Using PTMList,				
	
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
		}//*/
		
		if( ptmList != null ) {
			for(int k=0; k<ptmList.length; k++){
				if( ptmList[k] != null && ptmList[k].getNeutralLoss() > 0 ){					
					for( IonNode node : this ){	
						if( (node.type==0 && node.index>=k) || (node.type==1 && node.index<k) ){
							score += graph.getNormOfMatchedPeak( node.mz - ptmList[k].getNeutralLoss()/node.charge );
						}
					}
				}
			}
		}
		
		int i = 0;
		for( IonNode node : this ){	
			if( node.index == 1 && node.type == NTermType ) score += graph.getNormOfMatchedPeak( node.mz-Constants.A_ION_OFFSET );
			score += graph.getSupportingScore( observedPaeks.get(i++), node.type, node.charge, supportingFactor );	
		}

		double immSupprot = 0;
		for(int imm=0; imm<sequence.length(); imm++ ){
			immSupprot += graph.getNormOfMatchedPeak( MSMass.getAAMass(sequence.charAt(imm))+ptmMass[imm] + Constants.IMM_OFFSET );
			if( ptmList != null ) 
			if( ptmList[imm] != null && ptmList[imm].getDiagnosticIon() > 0 ) immSupprot += graph.getNormOfMatchedPeak( ptmList[imm].getDiagnosticIon() );
		}	
		score += immSupprot*0.5;
		
		double EPS = score - modifiedResd*Constants.rNorm[0];
		double UnEPS = graph.getPenaltySCoreForUnassignedIons(this.size());
		
		rankScore = (int)EPS;	
		coverScore = (int)(EPS-UnEPS);
	}
	
	public void evaluateMatchQuality( PGraph graph ){ 	
			
		setSecondaryScore( graph );
		setMatchCoverage();	
		
		double discriminant = 0;
		
		//SEQ
		if( charge < 3 ) discriminant= 0.2073*coverScore + 14.2075*seqCoverage + 0.9258*maxConsecutive - 17.2304;
		else discriminant= 0.1773*coverScore + 5.402*seqCoverage + 0.8193*maxConsecutive - 8.1191; //*/
		
		prob = Scoring.getOddProbability( discriminant );			
		
	}

}




