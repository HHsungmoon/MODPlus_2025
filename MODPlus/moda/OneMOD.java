package moda;

import processedDB.CandidateContainer;
import processedDB.MODPeptide;
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

public class OneMOD {

	public static int[] bestOnlineScore = new int[ThreadPoolManager.numSlots];

	public static DPHeap getHeatedPeptides( StemTagTrie stemDB, PGraph graph, TagPool tPool, boolean dynamicPMCorrection ){
		int slot = ThreadPoolManager.getSlotIndex();
		bestOnlineScore[slot] = 2;
		DPHeap annotation = null;

		for( TagTrie stem : stemDB ){
			CandidateContainer cpool= DBSearch.construct_onemod_cpool(tPool, graph.getCorrectedMW(), stem);

			DPHeap sanno = null;
			if( dynamicPMCorrection )
				sanno= run_dynamic_mass_mode( cpool, graph, stem );
			else
				sanno= run_static_mass_mode( cpool, graph, stem );

			if( annotation == null )
				annotation = sanno;
			else
				annotation.insertAll(sanno);
		}

		if( annotation.evaluate(graph) < 1 )
			return null;

		return annotation;
	}
	
	public static DPHeap run_static_mass_mode( CandidateContainer cpool, PGraph graph, TagTrie ixPDB ){

		PRM prmTable;
		if( graph.getCharge() > 2 ) prmTable= new PRMforHighCharge(graph);
		else prmTable= new PRM(graph);
		
		DPHeap topList= new DPHeap();
	
		MODPeptide[] DBList = cpool.getList();
		
		int poolsize = cpool.size();
		for(int i=0; i<poolsize; i++){			
			static_single_align(DBList[i], prmTable, ixPDB, topList);
		}
		topList.setStemNo(ixPDB.getStemNo());
		
		return topList;
	}
	public static DPHeap run_dynamic_mass_mode( CandidateContainer cpool, PGraph graph, TagTrie ixPDB ){
		
		PRM prmTable;
		if( graph.getCharge() > 2 ) prmTable= new PRMforHighCharge(graph);
		else prmTable= new PRM(graph);//*/
		
		DPHeap topList= new DPHeap();	

		MODPeptide[] DBList = cpool.getList();
		int poolsize = cpool.size();
		for(int i=0; i<poolsize; i++){			
			dynamic_single_align(DBList[i], prmTable, ixPDB, topList);
		}
		topList.setStemNo(ixPDB.getStemNo());

		return topList;
	}
	
	private static void static_single_align(MODPeptide entry, PRM prmTable, TagTrie ixPDB, DPHeap topList){	
		double 				observedMass = prmTable.getPeptMass();
		String 				peptide = entry.getPeptide(ixPDB);
		double 				peptMass = MSMass.getPepMass(peptide);
		
		int rowMax= 2, colMax= peptide.length()+1;		
		MatCell[][] specMatrix = new MatCell [rowMax][colMax];		

		for(int n=0; n<colMax ; n++){
			specMatrix[0][n]= new MatCell(2);
			specMatrix[1][n]= new MatCell(0);
		}

		DPPeptide temp= null;
		double delta = observedMass - peptMass - MODaConst.TERMINALMOD;
		double nTermDeletion= 0., cTermDeletion= 0.;
		int npi = 0, cpi = colMax;
		for(int i=entry.getStart(); i<=entry.getLeft(); i++){
			cTermDeletion= 0.;
			cpi = colMax;
			for(int j=entry.getEnd(); j>=entry.getRight(); j--) {

				int noOfET = ixPDB.getNTTOfPeptide(i, j, Constants.protease );
				if( noOfET >= Constants.numberOfEnzymaticTermini )
				{
					double massRange = delta + nTermDeletion + cTermDeletion;
					int slotIdx = ThreadPoolManager.getSlotIndex();
					if( (massRange < Constants.maxModifiedMass && massRange > Constants.minModifiedMass) ||
							Math.abs(massRange) < Constants.precursorTolerance[slotIdx] )
					{						
						double ptm = MODaConst.ptmUnit.getPtmMass(massRange);
						int intptm = Constants.round(massRange);
						double pmzErr = massRange - ptm;
						
						specMatrix[0][npi].setMass(0, 0., 0);//init start-columns
						specMatrix[0][npi].refresh();
						specMatrix[1][npi].setMass(0, ptm, intptm);
						specMatrix[1][npi].refresh();
						
						double cellMass = Constants.NTERM_FIX_MOD;
						for(int n=npi+1 ; n<cpi ; n++){
							cellMass += MSMass.getAAMass(peptide.charAt(n-1));
							specMatrix[0][n].setMass(cellMass, 0., 0);
							specMatrix[0][n].refresh();
							specMatrix[1][n].setMass(cellMass, ptm, intptm);
							specMatrix[1][n].refresh();
						}
						
						specMatrix[0][cpi-1].mass += Constants.CTERM_FIX_MOD;
						specMatrix[1][cpi-1].mass += Constants.CTERM_FIX_MOD;
						
						temp = dynamic_programming(peptide.substring(npi, cpi-1), observedMass-pmzErr, rowMax, npi, cpi, 
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
	
	private static void dynamic_single_align(MODPeptide entry, PRM prmTable, TagTrie ixPDB, DPHeap topList){	
		double 				observedMass = prmTable.getPeptMass();
		String 				peptide = entry.getPeptide(ixPDB);
		double 				peptMass = MSMass.getPepMass(peptide);

		int rowMax= 2, colMax= peptide.length()+1;		
		MatCell[][] specMatrix = new MatCell [rowMax][colMax];		

		for(int n=0; n<colMax ; n++){
			specMatrix[0][n]= new MatCell(2);
			specMatrix[1][n]= new MatCell(0);
		}

		DPPeptide temp= null;
		double delta = observedMass - peptMass - MODaConst.TERMINALMOD;
		double nTermDeletion= 0., cTermDeletion= 0.;
		int npi = 0, cpi = colMax;
		
		for(int i=entry.getStart(); i<=entry.getLeft(); i++){
			cTermDeletion= 0.;
			cpi = colMax;
			for(int j=entry.getEnd(); j>=entry.getRight(); j--)	{
				
				int noOfET = ixPDB.getNTTOfPeptide(i, j, Constants.protease );
				if( noOfET >= Constants.numberOfEnzymaticTermini )
				{
					double massRange = delta + nTermDeletion + cTermDeletion;
					int slotIdx = ThreadPoolManager.getSlotIndex();
					if( (massRange < Constants.maxModifiedMass && massRange > Constants.minModifiedMass) ||
							Math.abs(massRange) < Constants.precursorTolerance[slotIdx] )
					{
						double ptm = MODaConst.ptmUnit.getPtmMass(massRange);
						int intptm = Constants.round(massRange);
						double pmzErr = massRange - ptm;
						ptm    += MODaConst.maxIsotopeError;	
						intptm += MODaConst.maxIntIsotopeError;
						
						specMatrix[0][npi].setMass(0, 0., 0); //init start-columns
						specMatrix[0][npi].refresh();
						specMatrix[1][npi].setMass(0, ptm, intptm);
						specMatrix[1][npi].refresh();
							
						double cellMass = Constants.NTERM_FIX_MOD;
						for(int n=npi+1 ; n<cpi ; n++){
							cellMass += MSMass.getAAMass(peptide.charAt(n-1));
							specMatrix[0][n].setMass(cellMass, 0., 0);
							specMatrix[0][n].refresh();
							specMatrix[1][n].setMass(cellMass, ptm, intptm);
							specMatrix[1][n].refresh();
						}
						specMatrix[0][cpi-1].mass += Constants.CTERM_FIX_MOD;
						specMatrix[1][cpi-1].mass += Constants.CTERM_FIX_MOD;
					
						temp = DPwithMassCorrection(peptide.substring(npi, cpi-1), observedMass-pmzErr, rowMax, npi, cpi, 
								specMatrix, prmTable, pmzErr);
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
			MatCell[][] specMatrix, PRM prmTable, double pmzErr){	
		
		DPPeptide best= new DPPeptide(), temp= null;
		
		double massCorrection = MODaConst.maxIsotopeError;
		for(int dc=0; dc<MODaConst.isotopePointsToBeCorrected ; dc++){
			temp = dynamic_programming(peptide, obsMass+massCorrection, rowMax, smStart, smEnd, 
					specMatrix, prmTable, pmzErr-massCorrection);
			if( temp != null && temp.score > best.score ){
				best= temp;				
			}
			
			massCorrection += MODaConst.isotopeUnit;
			for(int n=smStart ; n<smEnd ; n++){					
				specMatrix[0][n].refresh();		
				specMatrix[1][n].refresh();		
				specMatrix[1][n].correctMass(MODaConst.isotopeUnit);//.mass += 1;
			}	
		}	
		return best;
	}
		
	static DPPeptide dynamic_programming(String peptide, double obsMass, int rowMax, int smStart, int smEnd, 
			MatCell[][] specMatrix, PRM prmTable, double pmzErr){	


		int colMax= smEnd-smStart, endingTag= 1;			
		if( specMatrix[endingTag][smStart].nominalDelta > Constants.maxModifiedMass ) return null;
	
		double upperLimit = obsMass+Constants.fragmentTolerance;
		specMatrix[0][smStart].isAAJump= 1;

		for (int n = smStart + 1; n < smEnd; n++) {
			for (int m = 0; m < rowMax; m += endingTag) {

				final MatCell currNode = specMatrix[m][n];
				final double currMass = currNode.mass;   // 필드 캐시(여러 번 쓰임)
				if (currMass > upperLimit) continue;

				double max = MODaConst.baseScore;
				MatCell bestParent = null;
				int bestJump = -1;

				for (int d = 0; d < rowMax; d += endingTag) {
					final MatCell prevNode = specMatrix[d][n - 1];

					// prev 필드도 로컬로 캐시해서 중복 로드 방지
					final int prevAA = prevNode.isAAJump;
					if (prevAA == -1) continue;

					final double prevScore = prevNode.score;

					if (m == d) {
						// same-row: tie 허용(<=) — 원본과 동일
						if (max <= prevScore) {
							max = prevScore;
							bestParent = prevNode;
							bestJump = 1;
						}
					} else if (m != 0) {
						// cross-row: 거리 조건 + 엄격 비교(<) — 원본과 동일
						final double prevMass = prevNode.mass;
						if (currMass - prevMass < MODaConst.minimumDistance) continue;
						if (max < prevScore) {
							max = prevScore;
							bestParent = prevNode;
							bestJump = 0;
						}
					}
				}

				if (max < 0) continue;
				if (bestParent != null) currNode.setParent(bestParent);
				if (bestJump   != -1)   currNode.isAAJump = bestJump;

				// getScore는 선택 확정 후 1회만 (원본 의미 동일)
				currNode.score = max + prmTable.getScore(currMass, pmzErr);
			}
		}
		
		MatCell initNode = specMatrix[0][smStart], tarNode = specMatrix[endingTag][smEnd-1];
		double idScore= ( tarNode.nominalDelta == 0 )? tarNode.score : tarNode.score - Constants.rNorm[0];
		int slot = ThreadPoolManager.getSlotIndex();
		if( idScore < bestOnlineScore[slot]/2 ) return null;
		
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

		if( idScore > bestOnlineScore[slot] )
			bestOnlineScore[slot] = (int)idScore;
		return new DPPeptide(peptide, (int)idScore, ptms, smStart);	
	}
	
}

