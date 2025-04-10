package moda;

import java.util.ArrayList;

import modi.ThreadPoolManager;
import processedDB.CandidateContainer;
import processedDB.ChainTagPeptide;
import processedDB.SequenceTag;
import processedDB.StemTagTrie;
import processedDB.TagTrie;

import msutil.MSMass;
import msutil.PGraph;
import msutil.PRM;
import msutil.PRMforHighCharge;

import modi.Constants;
import modi.Spectrum;
import modi.SpectrumAnalyzer;
import modi.TagPool;

public class MultiMOD {

	private static int bestOnlineScore = 2;
	
	public static DPHeap getHeatedPeptides( StemTagTrie stemDB, PGraph graph, TagPool tPool, boolean dynamicPMCorrection ){		
		bestOnlineScore = 2;
		DPHeap annotation = null;
		
		for( TagTrie stem : stemDB ){
			CandidateContainer cpool= DBSearch.construct_multimod_cpool(tPool, graph.getCorrectedMW(), stem);
						
			DPHeap sanno = null;
			if( dynamicPMCorrection ) 
				sanno= run_dynamic_mass_mode( cpool, graph, stem );
			else sanno= run_static_mass_mode( cpool, graph, stem );
		
			if( annotation == null ) annotation = sanno;
			else annotation.insertAll(sanno);		
		}
		if( annotation.reArr(graph) < 1 ) return null;

		return annotation;
	}

	public static DPHeap run_static_mass_mode( CandidateContainer cpool, PGraph graph, TagTrie ixPDB ){

		PRM prmTable;
		if( graph.getCharge() > 2 ) prmTable= new PRMforHighCharge(graph);
		else prmTable= new PRM(graph);	
		
		DPHeap topList= new DPHeap();
		ChainTagPeptide[] DBList = (ChainTagPeptide[])cpool.getList();

		for(int i=0; i<cpool.size(); i++){					
			static_multi_align(DBList[i], prmTable, ixPDB, topList);
		}			
		topList.setStemNo(ixPDB.getStemNo());
		
		return topList;
	}
	
	public static DPHeap run_dynamic_mass_mode( CandidateContainer cpool, PGraph graph, TagTrie ixPDB ){

		PRM prmTable;
		if( graph.getCharge() > 2 ) prmTable= new PRMforHighCharge(graph);
		else prmTable= new PRM(graph);
		
		DPHeap topList= new DPHeap();	
		
		ChainTagPeptide[] DBList = (ChainTagPeptide[])cpool.getList();
		for(int i=0; i<cpool.size(); i++){			
			dynamic_multi_align(DBList[i], prmTable, ixPDB, topList);
		}	
		topList.setStemNo(ixPDB.getStemNo());
		
		return topList;
	}
	
	private static void static_multi_align(ChainTagPeptide entry, PRM prmTable, TagTrie ixPDB, DPHeap topList){	
		double 				observedMass = prmTable.getPeptMass();
		String 				peptide = entry.getPeptide(ixPDB);
		double 				peptMass = MSMass.getPepMass(peptide);
		ArrayList<SequenceTag> 	mTagList = entry.getMatchedTags();
				
		int rowMax= mTagList.size()+2, colMax= peptide.length()+1;		
		MatCell[][] specMatrix = new MatCell [rowMax][colMax];		
		double[] deltas = new double[rowMax];
		deltas[0]= 0;	
		deltas[rowMax-1]= observedMass - peptMass - MODaConst.TERMINALMOD;
		
		for(int n=0; n<colMax ; n++){
			specMatrix[0][n]= new MatCell(2);
			specMatrix[rowMax-1][n]= new MatCell(0);
		}
		for(int m=1 ; m<rowMax-1 ; m++){
			SequenceTag st = mTagList.get(m-1);
			int n;
			for( n = 0 ; n <= st.getStart() ; n++ )
				specMatrix[m][n] = new MatCell(0);
			for( n = st.getStart()+1 ; n <= st.getEnd() ; n++ )
				specMatrix[m][n] = new MatCell(1);
			for( n = st.getEnd()+1 ; n < colMax ; n++ )
				specMatrix[m][n] = new MatCell(2);								
			deltas[m]= st.getNGap();
		}/// initialization done.*/
		
		DPPeptide temp= null;
		double nTermDeletion= 0., cTermDeletion= 0.;
		int npi = 0, cpi = colMax;
		
		for(int i=entry.getStart(); i<=entry.getLeft(); i++){
			cTermDeletion= 0.;
			cpi = colMax;
			for(int j=entry.getEnd(); j>=entry.getRight(); j--) {// MODPLUS
				
				int noOfET = ixPDB.getNTTOfPeptide(i, j, Constants.protease );
				if( noOfET >= Constants.numberOfEnzymaticTermini )
				{
					double massRange = deltas[rowMax-1] + nTermDeletion + cTermDeletion;
					int slot = ThreadPoolManager.getSlotIndex();
					if( (massRange < Constants.maxModifiedMass && massRange > Constants.minModifiedMass) || Math.abs(massRange) < Constants.precursorTolerance[slot] )
					{
						double[] ptms = new double[rowMax];
						int[] intptms = new int[rowMax];
						for(int k=1;k<rowMax-1; k++){
							ptms[k] = deltas[k] + nTermDeletion;
							intptms[k] = Constants.round(ptms[k]);													
						}
						ptms[rowMax-1] = massRange;
						intptms[rowMax-1] = Constants.round(ptms[rowMax-1]);
						MODaConst.ptmUnit.setPtmMasses(ptms, intptms);
						double pmzErr = massRange - ptms[rowMax-1];
						
						for(int m=0 ; m<rowMax ; m++){
							specMatrix[m][npi].setMass(0, ptms[m], intptms[m]);
							specMatrix[m][npi].refresh();
						}
						
						double cellMass = Constants.NTERM_FIX_MOD;
						for(int n=npi+1 ; n<cpi ; n++){
							cellMass += MSMass.getAAMass(peptide.charAt(n-1));
							for(int m=0 ; m<rowMax ; m++){
								specMatrix[m][n].setMass(cellMass, ptms[m], intptms[m]);
								specMatrix[m][n].refresh();
							}
						}
						for(int m=0 ; m<rowMax ; m++){						
							specMatrix[m][cpi-1].mass += Constants.CTERM_FIX_MOD;
						} 	

						temp = dynamicProgramming(peptide.substring(npi, cpi-1), observedMass-pmzErr, rowMax, npi, cpi, 
								specMatrix, prmTable, pmzErr);
						if( temp != null ) {
							temp.setProteinAndNTT(entry.getStart(), noOfET);
							topList.insert( temp );
						}
					}	
				}
				cpi--;	
				if( cpi - npi < 4 ) break;
				cTermDeletion += MSMass.getAAMass(peptide.charAt(cpi-1));
			}
			nTermDeletion += MSMass.getAAMass(peptide.charAt(npi++));
		}	
	}
	
	private static void dynamic_multi_align(ChainTagPeptide entry, PRM prmTable, TagTrie ixPDB, DPHeap topList){	
		double 				observedMass = prmTable.getPeptMass();
		String 				peptide = entry.getPeptide(ixPDB);
		double 				peptMass = MSMass.getPepMass(peptide);
		ArrayList<SequenceTag> 	mTagList = entry.getMatchedTags();
		
		int rowMax= mTagList.size()+2, colMax= peptide.length()+1;		
		MatCell[][] specMatrix = new MatCell [rowMax][colMax];		
		double[] deltas = new double[rowMax];
		int[] ionType = new int[rowMax];
		deltas[0]= 0;	
		ionType[0]= 0;
		deltas[rowMax-1]= observedMass - peptMass - MODaConst.TERMINALMOD;
		ionType[rowMax-1]= 1;
		
		for(int n=0; n<colMax ; n++){
			specMatrix[0][n]= new MatCell(2);
			specMatrix[rowMax-1][n]= new MatCell(0);
		}
		for(int m=1 ; m<rowMax-1 ; m++){
			SequenceTag st = mTagList.get(m-1);
			int n;
			for( n = 0 ; n <= st.getStart() ; n++ )
				specMatrix[m][n] = new MatCell(0);
			for( n = st.getStart()+1 ; n <= st.getEnd() ; n++ )
				specMatrix[m][n] = new MatCell(1);
			for( n = st.getEnd()+1 ; n < colMax ; n++ )
				specMatrix[m][n] = new MatCell(2);								
			deltas[m]= st.getNGap();
			ionType[m]= st.getType();
		}/// initialization done.*/
		
		DPPeptide temp= null;
		double nTermDeletion= 0., cTermDeletion= 0.;
		int npi = 0, cpi = colMax;
		
		for(int i=entry.getStart(); i<=entry.getLeft(); i++){
			cTermDeletion= 0.;
			cpi = colMax;
			for(int j=entry.getEnd(); j>=entry.getRight(); j--) {
				
				int noOfET = ixPDB.getNTTOfPeptide(i, j, Constants.protease );
				if( noOfET >= Constants.numberOfEnzymaticTermini )
				{
					double massRange = deltas[rowMax-1] + nTermDeletion + cTermDeletion;
					int slot = ThreadPoolManager.getSlotIndex();
					if( (massRange < Constants.maxModifiedMass && massRange > Constants.minModifiedMass) || Math.abs(massRange) < Constants.precursorTolerance[slot] )
					{
						double[] ptms = new double[rowMax];
						int[] intptms = new int[rowMax];
						for(int k=1;k<rowMax-1; k++){
							ptms[k] = deltas[k] + nTermDeletion;
							intptms[k] = Constants.round(ptms[k]);													
						}
						ptms[rowMax-1] = massRange;
						intptms[rowMax-1] = Constants.round(ptms[rowMax-1]);
						MODaConst.ptmUnit.setPtmMasses(ptms, intptms);
						double pmzErr = massRange - ptms[rowMax-1];

						for(int m=0 ; m<rowMax ; m++){
							if( ionType[m] == 1 ) {
								ptms[m]    += MODaConst.maxIsotopeError;
								intptms[m] += MODaConst.maxIntIsotopeError;
							}
							specMatrix[m][npi].setMass(0, ptms[m], intptms[m]);
							specMatrix[m][npi].refresh();
						}
						
						double cellMass = Constants.NTERM_FIX_MOD;
						for(int n=npi+1 ; n<cpi ; n++){
							cellMass += MSMass.getAAMass(peptide.charAt(n-1));
							for(int m=0 ; m<rowMax ; m++){
								specMatrix[m][n].setMass(cellMass, ptms[m], intptms[m]);	
								specMatrix[m][n].refresh();
							}
						}
						for(int m=0 ; m<rowMax ; m++){						
							specMatrix[m][cpi-1].mass += Constants.CTERM_FIX_MOD;
						}
						
						temp = DPwithMassCorrection(peptide.substring(npi, cpi-1), observedMass-pmzErr, rowMax, npi, cpi, 
								specMatrix, prmTable, pmzErr, ionType);

						temp.setProteinAndNTT(entry.getStart(), noOfET);
						topList.insert( temp );
					}	
				}
				cpi--;
				if( cpi - npi < 4 ) break;
				cTermDeletion += MSMass.getAAMass(peptide.charAt(cpi-1));
			}
			nTermDeletion += MSMass.getAAMass(peptide.charAt(npi++));
		}	
	}
	
	static DPPeptide DPwithMassCorrection(String peptide, double obsMass, int rowMax, int smStart, int smEnd, 
			MatCell[][] specMatrix, PRM prmTable, double pmzErr, int[] ionType){	
		
		DPPeptide best= new DPPeptide(), temp= null;
		
		double massCorrection = MODaConst.maxIsotopeError;
		for(int dc=0; dc<MODaConst.isotopePointsToBeCorrected ; dc++){
			temp = dynamicProgramming(peptide, obsMass+massCorrection, rowMax, smStart, smEnd, 
					specMatrix, prmTable, pmzErr-massCorrection);
			
			if( temp != null && temp.score > best.score ){
				best= temp;				
			}

			massCorrection += MODaConst.isotopeUnit;	
			for(int m=0 ; m<rowMax ; m++){// set mass for matrix	
				for(int n=smStart ; n<smEnd ; n++){					
					specMatrix[m][n].refresh();					
					if( ionType[m] == 1 ) specMatrix[m][n].correctMass(MODaConst.isotopeUnit);
				}	
			}	
		}	
		return best;
	}
	
	static DPPeptide dynamicProgramming(String peptide, double obsMass, int rowMax, int smStart, int smEnd, 
			MatCell[][] specMatrix, PRM prmTable, double pmzErr){	
		
		int colMax= smEnd-smStart;
		
		double upperLimit = obsMass+Constants.fragmentTolerance;
		specMatrix[0][smStart].isAAJump= 1;
		
		MatCell currNode, prevNode;
		for(int n=smStart+1 ; n<smEnd ; n++){
			for(int m=0; m<rowMax ; m++){

				currNode = specMatrix[m][n];
				if( currNode.mass > upperLimit ) continue;	
				
				double max = MODaConst.baseScore;			
				for(int d=0; d<rowMax ; d++){

					prevNode = specMatrix[d][n-1];
					if( prevNode.isAAJump == -1 ) continue;
					
					if( m == d ){// AA Jump
						if( max <= prevNode.score ){ 
							max= prevNode.score;			
							currNode.setParent(prevNode);
							if( currNode.isInsideTag > 0 ) currNode.isAAJump = 1;
							else currNode.isAAJump = 0;								
						}
					}
					else{ // Modification Jump
						if( currNode.canObliqueJumpFrom(prevNode) ){
							if( currNode.mass - prevNode.mass < MODaConst.minimumDistance ||
									currNode.nominalDelta - prevNode.nominalDelta > Constants.maxModifiedMass  ) continue;	
							
							double prev = ( currNode.nominalDelta == prevNode.nominalDelta ) ? prevNode.score : prevNode.score - Constants.rNorm[0];		
							if( max < prev ){
								max= prev;		
								currNode.setParent(prevNode);
								currNode.isAAJump = 0;								
							}
						}
					}
				}			
				if( max == MODaConst.baseScore ) continue;			
				currNode.score= max + prmTable.getScore(currNode.mass, pmzErr);
			}				
		}
		
		// back tracking
		MatCell initNode = specMatrix[0][smStart], tarNode = specMatrix[rowMax-1][smEnd-1];
		if( tarNode.score < bestOnlineScore/2 ) return null;
		double idScore= tarNode.score;
		
		int modifiedSite = 0;
		double[] ptms= new double[colMax-1];
		double[] matchedList = new double[colMax];	
		int ixx = colMax-1;	
		while( tarNode != initNode ){			
			matchedList[ixx--]= tarNode.mass;
			if( tarNode.nominalDelta != tarNode.parent.nominalDelta ){		
				ptms[ixx] = tarNode.delta - tarNode.parent.delta;
				modifiedSite++;
			}			
			tarNode= tarNode.parent;			
		}
		
		//Solve Symmetric path problem
		int forward=1, backward= colMax-2;	
		int symMatch = 0;
		double PMCorr = obsMass+Constants.H2O;
		while( forward < backward ){
			double symmetric = matchedList[forward] + matchedList[backward];
			if( Constants.fEqual(symmetric, PMCorr) ){		
				idScore -= prmTable.getScore( matchedList[forward], pmzErr );
				symMatch++;
				forward++;
				backward--;				
			}
			else if( symmetric > PMCorr )
				backward--;
			else forward++;		
		}
		
		if( idScore > bestOnlineScore ) bestOnlineScore = (int)idScore;		
		DPPeptide identification= new DPPeptide(peptide, (int)idScore, ptms, smStart );	
		if( modifiedSite > 1 && symMatch > 0 ){
			DPPeptide temp = dynamicProgrammingWithoutTags(peptide, obsMass, rowMax, smStart, smEnd, specMatrix, 
					prmTable, pmzErr);

			if( temp != null && identification.compareTo( temp ) == 1 )
				identification = temp;
		}
		return identification;		
	}
		
	static DPPeptide dynamicProgrammingWithoutTags(String peptide, double obsMass, int rowMax, int smStart, int smEnd, 
			MatCell[][] specMatrix, PRM prmTable, double pmzErr){	
		
		int colMax= smEnd-smStart, endingTag= rowMax-1;			
		if( specMatrix[endingTag][smStart].nominalDelta > Constants.maxModifiedMass ) return null;
		for(int n=smStart ; n<smEnd ; n++){
			specMatrix[endingTag][n].refresh();
		}
		
		double upperLimit = obsMass+Constants.fragmentTolerance;
		specMatrix[0][smStart].isAAJump= 1;
		
		MatCell currNode, prevNode;
		for(int n=smStart+1 ; n<smEnd ; n++){	
			for(int m=0; m<rowMax ; m+=endingTag ){
				
				currNode = specMatrix[m][n];
				if( currNode.mass > upperLimit ) continue;				
				
				double max = MODaConst.baseScore;	
				for(int d=0; d<rowMax ; d+=endingTag){
					
					prevNode = specMatrix[d][n-1];
					if( prevNode.isAAJump == -1 ) continue;
					
					if( m == d ){
						if( max <= prevNode.score ){ 
							max= prevNode.score;			
							currNode.setParent(prevNode);
							currNode.isAAJump = 1;	
						}
					}
					else if( m != 0 ) {					
						if( currNode.mass - prevNode.mass < MODaConst.minimumDistance ) continue;
						
						if( max < prevNode.score ){
							max= prevNode.score;		
							currNode.setParent(prevNode);
							currNode.isAAJump = 0;
						}
					}
				}			
				if( max < 0 ) continue;			
				currNode.score= max + prmTable.getScore(currNode.mass, pmzErr);
			}				
		}
		
		MatCell initNode = specMatrix[0][smStart], tarNode = specMatrix[endingTag][smEnd-1];
		double idScore= ( tarNode.nominalDelta == 0 )? tarNode.score : tarNode.score - Constants.rNorm[0];
		if( idScore < bestOnlineScore/2 ) return null;
		
		double[] ptms= new double[colMax-1];
		double[] matchedList = new double[colMax];	
		int ixx = colMax-1;				
		while( tarNode != initNode ){			
			matchedList[ixx--]= tarNode.mass;
			if( tarNode.nominalDelta != tarNode.parent.nominalDelta ){		
				ptms[ixx] = tarNode.delta - tarNode.parent.delta; 
			}			
			tarNode= tarNode.parent;			
		}
		
		//Solve Symmetric path problem
		int forward=1, backward= colMax-2;	
		double PMCorr = obsMass+Constants.H2O;
		while( forward < backward ){
			double symmetric = matchedList[forward] + matchedList[backward];
			if( Constants.fEqual(symmetric, PMCorr) ){
				idScore -= prmTable.getScore( matchedList[forward], pmzErr );		
				forward++;
				backward--;				
			}
			else if( symmetric > PMCorr  )
				backward--;
			else forward++;		
		}
		if( idScore > bestOnlineScore ) bestOnlineScore = (int)idScore;	
		return new DPPeptide(peptide, (int)idScore, ptms, smStart );	
	}
	
}






