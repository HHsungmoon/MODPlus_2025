package modi;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;

import msutil.PGraph;

public class Gap implements SpecInterpretation {
	Spectrum			sourceSpectrum;
	Peptide				matchedPeptide;
	int					start, end;		// both are inclusive
	double				offset;
	double				bStartMass;
	double				yStartMass;		// since ver 1.0
	PTMPosition			position;		// Nterm, C
	GapInterpretation	gapInterpretation;		
	// since ver 1.0
	boolean				interpreted;
	int					interpretID;
	
	double				mostAbundantBYPeakIntensity;
	

	public				Gap(Peptide matchedPeptide, int start, int end, 
			double bStartMass, double yStartMass, double offset, 
			PTMPosition position, Spectrum sourceSpectrum, double mostAbundantBYPeakIntensity)
	{
		this.matchedPeptide = matchedPeptide;
		this.start = start;
		this.end = end;
		this.bStartMass = bStartMass;
		this.yStartMass = yStartMass;
		this.offset = offset;
		this.position = position;
		this.sourceSpectrum = sourceSpectrum;
		this.mostAbundantBYPeakIntensity = mostAbundantBYPeakIntensity;
	}

	
	public	void		setInterpretation( PTMSearchResult result, PGraph graph )	
	{
		this.gapInterpretation = new GapInterpretation( this, result.getPTMRun(), graph );
		this.interpretID = result.getID();
	}
	
	public	void		setInterpreted(boolean interpreted) { this.interpreted = interpreted; }
	public	double		getOffset()			{ return offset; }

	public	Peptide		getMatchedPeptide() { return matchedPeptide; }
	public	Sequence	sequence() 			{ return matchedPeptide.subSequence(start, end+1); }
	public	int			getStart()			{ return start; }
	public	int			getEnd()			{ return end; }
	public	int			getLength()			{ return end - start + 1; }
	public	PTMPosition	getPosition()		{ return position; }
	public	GapInterpretation getGapInterpretation()	{ return gapInterpretation; }
	public	int			compareTo(SpecInterpretation t)	{ return this.start - t.getStart(); }

	public	String		getSequenceStr()	
	{ 
		StringBuffer output = new StringBuffer();
		for(int i=start; i<=end; i++)
			output.append('_');
		
		return output.toString();
	}
	// should be optimized later
	public	double		getScore()			
	{ 
		if(gapInterpretation == null || gapInterpretation.size() == 0)
			return 0.;
		else
			return gapInterpretation.getGapScore(0);
	}

	
	public ArrayList<Peak> getTheoreticalPeaks(PTMRun run)
	{
		double [] ptmMass = new double[end-start+1];
		for(PTMOccurrence occr : run)
			ptmMass[occr.getPosition()] += occr.getPTM().getMassDifference();

		ArrayList<Peak> allPeaks = new ArrayList<Peak>();
		
		double mass = bStartMass, compMass= sourceSpectrum.getCorrectedMW() - bStartMass + 2*Constants.Proton;
		
		for(int i=start; i<end; i++) 
		{
			double addMass = matchedPeptide.get(i).getMass() + ptmMass[i-start];
			mass += addMass;
			compMass -= addMass;
			
			AminoAcid NTermResidue = matchedPeptide.get(i);
			AminoAcid CTermResidue = matchedPeptide.get(i+1);
			if(ptmMass[i-start] != 0)
				NTermResidue = AminoAcid.getAminoAcid(23);	// modified
			if(ptmMass[i-start+1] != 0)
				CTermResidue = AminoAcid.getAminoAcid(23);  // modified//*/

		//	if( type != GapType.btnCC )
			{
				allPeaks.add(new MatchedPeak(-1, mass, 0, 1, PeakProperty.B_ION, NTermResidue, CTermResidue));
				if( sourceSpectrum.getCharge() > 2 ){
					allPeaks.add(new MatchedPeak(-1, (mass+1)/2, 0, 2, PeakProperty.B_ION, NTermResidue, CTermResidue));
				}
			}

		//	if( type != GapType.btnNN )
			{
				allPeaks.add(new MatchedPeak(-1, compMass, 0, 1, PeakProperty.Y_ION, NTermResidue, CTermResidue));
				if( sourceSpectrum.getCharge() > 2 ){
					allPeaks.add(new MatchedPeak(-1, (compMass+1)/2, 0, 1, PeakProperty.Y_ION, NTermResidue, CTermResidue));
				}
			}
		}
		
		// add diagnostic peak
		HashSet<PTM> checkAdded = new HashSet<PTM>();
		DiagnosticPTMIonDB diagPTMDB = new DiagnosticPTMIonDB();
		for(PTMOccurrence occr : run)
		{
			PTM ptm = occr.getPTM();
			if(!checkAdded.contains(ptm))
			{
				checkAdded.add(ptm);
				for(DiagnosticPTMIon ion : diagPTMDB.getDiagIon(ptm)){
					allPeaks.add( new DiagnosticPTMPeak(ion.getMass(), ion.getProb()) );
				}
			}
		}
		
		Collections.sort(allPeaks);
		
		return allPeaks;
	}

	public ArrayList<Peak> getTheoreticalPeaks()
	{
		if(gapInterpretation == null || gapInterpretation.size() == 0)
			return getTheoreticalPeaks(new PTMRun());
		else
			return getTheoreticalPeaks(gapInterpretation.get(0));
	}
	
	public	static String getTheoreticalMassStr(ArrayList<Peak> peakList)
	{
		StringBuffer output = new StringBuffer();
		boolean start = true;
		for(Peak p : peakList)
		{
			if(start)
				start = false;
			else
				output.append(",");
			output.append(Constants.getString(p.getMass()));
		}
		return output.toString();
	}
	
	public double getScore(PTMRun run)
	{
		ArrayList<Peak> theoPeaks = this.getTheoreticalPeaks(run);
		
		double score = ( run.getError() > Constants.fragmentTolerance )? -0.5 : 0;

		ArrayList<PeakPair> sharedPeaks = new PeakListComparator(sourceSpectrum, theoPeaks).getSharedPeaks();
		for(PeakPair p : sharedPeaks)
		{
			PeakProperty property = p.getSecond().getPeakProperty();
			double weight;
			if(property == PeakProperty.DIAGNOSTIC_PTM_ION)
				weight = ((DiagnosticPTMPeak)p.getSecond()).prob;
			else
				weight= 1;
			score += weight*p.getFirst().getNormIntensity();
		}
		
		score -= run.size()*0.5;//Constants.PTM_ADD_PENALTY;
	//	score -= run.getPenalty()*0.5;//Constants.PTM_ADD_PENALTY;
		
		return score;
	}
	
	public boolean isReasonable(){
		if( !Constants.isInModifiedRange( this.offset ) ) return false;	
		
	//	if ( this.offset > 0 ) return true;
	/*	// prevent negative AA
		double coveredMass= matchedPeptide.getMonoMass(start, end+1);
		if( start == 0 ) coveredMass += Constants.NTERM_FIX_MOD;
		if( end+1 == matchedPeptide.size() ) coveredMass += Constants.CTERM_FIX_MOD;
		if( coveredMass + this.offset <= (end-start+1)*18-Constants.gapTolerance ) return false;//*/
		
		int pos = 0; //because no. of allowed mods is max two
		if( start == 0 ) pos = 1;
		if( end+1 == matchedPeptide.size() ) pos = 2;
		Sequence seq = matchedPeptide.subSequence(start, end+1);
		double minMod = Constants.variableModifications.minimumModifiedMass(seq, pos);		
		double maxMod = Constants.variableModifications.maximumModifiedMass(seq, pos);		
	//	System.out.println( this.offset + " " + maxMod + " " +minMod + " " + matchedPeptide.subSequence(start, end+1));
		if( this.offset < (minMod - Constants.gapTolerance) || this.offset > (maxMod + Constants.gapTolerance) ) return false;//*/
		
		return true;
	}
	
	public int predictModifiedSite()
	{// for novel ptms
		if( start == end ){	// gap length 1	
			if ( start==0 && Constants.fEqual(matchedPeptide.get(0).getMass()+this.offset, 0) )
				return -1;
			return start;
		}
		
		int i;
		double mass=0, compMass=0;
		double localScore=0, maxScore=0;
		int maxSite= -1;
		double [] scoreList = new double[end-start+2]; //for one, random score
	
		for(i=start; i<=end; i++) {		
			double [] ptmMass = new double[end-start+1];
			ptmMass[i-start] = this.offset;
			mass = bStartMass;
			localScore= 0;
			
			for(int k=start; k<end; k++) {
				mass = mass + matchedPeptide.get(k).getMass() + ptmMass[k-start];				
				if( Constants.fEqual(mass, Constants.B_ION_OFFSET) )
					return -1;
				compMass = sourceSpectrum.getCorrectedMW() - mass + 2*Constants.UNIT_MASS;
				if( Constants.fEqual(compMass, Constants.Y_ION_OFFSET) )
					return -2;
				localScore += sourceSpectrum.getMatchedPeak(mass);
				localScore += sourceSpectrum.getMatchedPeak(compMass);
			}
			mass = mass + matchedPeptide.get(end).getMass() + this.offset;				
			if( Constants.fEqual(mass, Constants.B_ION_OFFSET) )
				return -1;
		
			
			scoreList[i-start] = localScore;
			if( localScore > maxScore ){
				maxScore= localScore;
				maxSite= i;
			}
		}
		
		if ( start==0 && maxSite < 2 ) // nterm gap
			return 0;	
		
		scoreList[i-start] = sourceSpectrum.getRandomPeakIntensity() * (end-start);
		Arrays.sort(scoreList);
		double delt = scoreList[scoreList.length-2];
		
		delt= 1 - (delt/maxScore);
		if( delt > 1.0/(end-start+1) )
		 return maxSite;
		else
			return -10; //*/
	}
	
	public double getScore(PTMRun run, PGraph graph)
	{
		double score = - run.getError();
		double [] ptmMass = new double[end-start+1];
		for(PTMOccurrence occr : run)
			ptmMass[occr.getPosition()] += occr.getPTM().getMassDifference();
		
		graph.refresh();
		ArrayList<Integer> singly = new ArrayList<Integer>();
		ArrayList<Integer> doubly = new ArrayList<Integer>();
		
		double mass = bStartMass, compMass= sourceSpectrum.getCorrectedMW() - bStartMass + 2*Constants.Proton;
		int count = 0;
		for(int i=start; i<end; i++)
		{
			double addMass = matchedPeptide.get(i).getMass() + ptmMass[i-start];
			mass += addMass;
			compMass -= addMass;
			count++;
			
			int index = graph.getIndexOfMatchedPeak( mass );
			if( index != -1 ){
				singly.add(index);
				score += graph.get(index).getNorm();
			}
			index = graph.getIndexOfMatchedPeak( compMass );
			if( index != -1 ){
				singly.add(index);
				score += graph.get(index).getNorm();
			}
			if( sourceSpectrum.getCharge() > 2 ){
				index = graph.getIndexOfMatchedPeak( (mass+Constants.Proton)/2 );
				if( index != -1 ){
					doubly.add(index);
					score += graph.get(index).getNorm();
				}
				index = graph.getIndexOfMatchedPeak( (compMass+Constants.Proton)/2 );
				if( index != -1 ){
					doubly.add(index);
					score += graph.get(index).getNorm();
				}
			}
		}
		for(int i=0; i<singly.size();i++){
			score += graph.getSupportingScore( singly.get(i++), 1, 1, 0.5 );
		}
		for(int i=0; i<doubly.size();i++){
			score += graph.getSupportingScore( doubly.get(i++), 1, 2, 0.5 );
		}
	//	score -= run.size()*2;
		score -= run.getPenalty()*2;
		return score;
	}
	
	
}