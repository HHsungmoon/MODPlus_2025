package modi;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.TreeSet;

import msutil.PGraph;
import msutil.PNode;

public class Spectrum extends ArrayList<Peak> implements Comparable<Spectrum> {

	private String		name = "";	// id
	private double 		precursor;
	private int 		charge;
	private double 		observedMW;	// not used now
	private double 		correctedMW;
	private double		TIC= 0;
	private double		RTTime= 0;
	private double		bpIntensity= 0;
	private double		randomPeakIntensity= 0;
	
	private	ArrayList<Peak> 	selectedPeak = null;
	
	public	Spectrum(double precursor, int charge) 
	{
		correctedMW = observedMW = (precursor - Constants.Proton)*charge;
		this.precursor = precursor;
		this.charge = charge;
	}
	public 	Spectrum(double precursor, int charge, String name)
	{
		this(precursor, charge);
		this.name = name;
	}
	public	void 	setExtraInformation(double bpi, double t){
		bpIntensity= bpi;
		TIC= t;
	}

	public	double	getObservedMW()	{ return observedMW; }
	public	double	getCorrectedMW()	{ return correctedMW; }
	public	void	setCorrectedParentMW(double mw)	{ 
		correctedMW = mw; 
	}
	
	public	double	getPrecursor()			{ return precursor; }
	public	int		getCharge() 			{ return charge; }
	public	String	getName()				{ return name; }
	public	double	getRandomPeakIntensity() { return randomPeakIntensity; }

	public void normalizeIntensityLocally(){
		for(Peak p : this){			
			p.setNormIntensity(p.getIntensity()/getLocalMaxPeak(p.getMass()));	
		}
	}
	
	public double getLocalMaxPeak(double center){
		double lMax =0;
		int bin = 50;
		int index = binarySearchforPeaks( center-bin );
		
		while( this.get(index).getMass() <= center+bin ){
			if( this.get(index).getIntensity() > lMax )
				lMax= this.get(index).getIntensity();
			index++;
			if( index == this.size() )
				break;
		}	
		return lMax;
	}

	public PGraph getPeakGraph(){ //MOD SERIES
		
		int binSize= 100, considered= 10;
		double minimumCut= bpIntensity*0.001;
		int slot = ThreadPoolManager.getSlotIndex();
		double precursoeRange= Constants.precursorTolerance[slot]/this.charge + Constants.fragmentTolerance;
		
		PGraph graph = new PGraph(this.observedMW, charge);
		
		ArrayList<Peak> locBin = new ArrayList<Peak>();
		
		int lowerBound = (int)this.get(0).getMass()+100, upperBound = (int)this.get(this.size()-1).getMass()-100;
		
		int startLocal = (int)this.get(0).getMass()/binSize;		
		int local= startLocal;
		
		for( Peak p : this ){
			if( local == (int)p.getMass()/binSize )
				locBin.add(p);
			else{
				Collections.sort(locBin, Collections.reverseOrder( new IntensityComparator()) );
				
				int to= 1;
				for( Peak lp : locBin ){			
					if( startLocal == local && lp.intensity < minimumCut ) break;
					if( lp.mass < lowerBound || lp.mass > upperBound ) {
						if( lp.intensity < minimumCut ) {
							to++;
							if( considered < to ) break;
							continue;
						}//*/
					}
					
					if( this.precursor - lp.mass < precursoeRange && lp.mass - this.precursor < 2. ) {
						graph.add( new PNode(lp.mass, minimumCut, 10, Constants.rNorm[10]) );		
					}
					else //*/
					{
						graph.add( new PNode(lp.mass, lp.intensity, to, Constants.rNorm[to]) );					
						to++;
					}
					if( considered < to ) break;
				}
				
				local = (int)p.getMass()/binSize;
				locBin.clear();
				locBin.add(p);
			}
		}
		
		Collections.sort(locBin, Collections.reverseOrder( new IntensityComparator()) );
		int to= 1;
		for( Peak lp : locBin ){		
			if( lp.intensity < minimumCut ) break;				
			graph.add( new PNode(lp.mass, lp.intensity, to, Constants.rNorm[to]) );
			to++;
			if( considered < to ) break;
		}//Last Region	

		graph.ready(TIC);
		
		return graph;
	}

	public 	String	toString() 
	{
		StringBuffer buffer = new StringBuffer();
		buffer.append("Spectrum("+precursor+","+charge+")\n");
		Iterator e = this.iterator();
		while(e.hasNext())
		{
			buffer.append(((Peak)e.next()).toString());
		}
		return buffer.toString();
	}
	
	public int peakSelection (double selectionWindowSize, int numOfPeaksInWindow ){
		int globalSelectedSize = (int) (this.observedMW / selectionWindowSize * numOfPeaksInWindow);
		return peakSelection(globalSelectedSize, selectionWindowSize, numOfPeaksInWindow);
	}
	public int peakSelection (int globalSelectedSize, double selectionWindowSize, int numOfPeaksInWindow ){
		assert(this.checkSorted());
		if(globalSelectedSize > this.size())
			globalSelectedSize = this.size();
		if(this.size() == 0)
			return -1;	// later : error code
		// global selection
		Iterator<Peak> it = this.iterator();
		TreeSet<Peak> globalSelected = new TreeSet<Peak>(new IntensityComparator());
		for(int i=0; it.hasNext() && i<globalSelectedSize; i++)
			globalSelected.add(it.next());
		TreeSet<Peak> notSelected = new TreeSet<Peak>(new MassComparator());
		Peak curPeak;
		while(it.hasNext()) 
		{
			curPeak = it.next();
			if(curPeak.getIntensity() > globalSelected.first().getIntensity())
			{
				notSelected.add(globalSelected.first());
				globalSelected.remove(globalSelected.first());	// delete smallest peak
				globalSelected.add(curPeak);
			}
			else
				notSelected.add(curPeak);
		}
		
		// local selection
		TreeSet<Peak> selected = new TreeSet<Peak> (new MassComparator());
		selected.addAll(globalSelected);
		int currentBinPeakSize;
		double currentBinStart, currentBinEnd;
		// do not select additional peaks in (0, selectionWindowSize)
		for(int i=2; i<(int)(correctedMW/selectionWindowSize+1)*2; i++)
		{
			currentBinStart = i*selectionWindowSize/2;
			currentBinEnd = currentBinStart + selectionWindowSize;
			currentBinPeakSize = selected.subSet(new Peak(0, currentBinStart, 0), new Peak(0, currentBinEnd, 0)).size();
			if( currentBinPeakSize < numOfPeaksInWindow ){
				Peak [] newPeaks = notSelected.subSet(new Peak(0, currentBinStart, 0.), new Peak(0, currentBinEnd, 0.)).toArray(new Peak[0]);
				Arrays.sort(newPeaks, Collections.reverseOrder(new IntensityComparator()));
				for(int j=0; j<newPeaks.length && j < numOfPeaksInWindow - currentBinPeakSize; j++)
				{
					if( newPeaks[j].getNormIntensity() > Constants.minNormIntensity ){
						notSelected.remove(newPeaks[j]);
						selected.add(newPeaks[j]);
					}
					else break;
				} 
			}
		}
		
		selected.add(new Peak(-1, Constants.B_ION_OFFSET+Constants.NTERM_FIX_MOD, 0, 1, PeakProperty.N_TERM_B_ION_ONLY));
		selected.add(new Peak(-1, Constants.Y_ION_OFFSET+Constants.CTERM_FIX_MOD, 0, 1, PeakProperty.C_TERM_Y_ION_ONLY));
		selected.add(new Peak(-1, correctedMW-Constants.H2O+Constants.Proton-Constants.CTERM_FIX_MOD, 0, 1, PeakProperty.C_TERM_B_ION_ONLY));
		selected.add(new Peak(-1, correctedMW+Constants.Proton-Constants.NTERM_FIX_MOD, 0, 1, PeakProperty.N_TERM_Y_ION_ONLY));

		selectedPeak = new ArrayList<Peak>(selected);		
		setScoreOfSelectedPeaks(selectedPeak, 1, Constants.massToleranceForDenovo);
	/*	for(int si=0; si<selectedPeak.size(); si++)
			System.out.println("P: " + selectedPeak.get(si).mass + " / " + selectedPeak.get(si).intensity);
		for(int si=0; si<this.size(); si++)
			System.out.println("x: " + this.get(si).mass + " / " + this.get(si).intensity);//*/
		
		return selectedPeak.size();
	}	

	public	TagPool	generateTags( int minTagLength, int minTagLengthPeptideShouldContain, double massTolerance ) {// by NA
		// selectedPeak array should be sorted before generating Tags
		assert(selectedPeak != null);
	//	double maxGap = AminoAcid.getAminoAcid('W').getMonoMass() + massTolerance; // Mass of 'W' is biggest
		double maxGap = AminoAcid.getMaxMonoMass() + massTolerance; // Mass of 'W' is biggest
		
		TagPool tags = new TagPool();		
		for(int i=0; i<selectedPeak.size()-1; i++)
		{		
			for(int j=i+1; j<selectedPeak.size(); j++){
				double diff = Peak.getMassDifference(selectedPeak.get(i), selectedPeak.get(j));
				if( diff > maxGap )
					break;
				
				// one length tag
				ArrayList<AminoAcid> codeList = AminoAcid.getCode( diff, massTolerance );//massTolerance);
				for(AminoAcid code : codeList) {
					tags.add(new Tag(selectedPeak.get(i), selectedPeak.get(j), code, this));
				}							
			}			
		}
		
		TagPool oneLengthTagPool 	= (TagPool)tags.clone();
		TagPool tempTags			= (TagPool)tags.clone();
		
		for(int i=1; i<minTagLengthPeptideShouldContain; i++)
		{
			tempTags = TagPool.extendTags(tempTags, oneLengthTagPool);		
			tags.addAll(tempTags);
		}//*/		
		
		if( this.charge > 2 ){
			tags.addAll( this.generateDoublyTags(minTagLength, minTagLengthPeptideShouldContain, massTolerance*2) );			
		}//*/
		
		tags.setTagScores();
		tags = tags.extractQualifiedTag(minTagLength, Constants.MAX_TAG_SIZE);//*/	
	
	/*	for(int j=0; j<selectedPeak.size(); j++){
			if( selectedPeak.get(j).getPeakProperty() == PeakProperty.N_TERM_B_ION_ONLY ) selectedPeak.get(j).shiftMass(-Constants.NTERM_FIX_MOD);
			if( selectedPeak.get(j).getPeakProperty() == PeakProperty.C_TERM_Y_ION_ONLY ) selectedPeak.get(j).shiftMass(-Constants.CTERM_FIX_MOD);
			if( selectedPeak.get(j).getPeakProperty() == PeakProperty.C_TERM_B_ION_ONLY ) selectedPeak.get(j).shiftMass(Constants.CTERM_FIX_MOD);
			if( selectedPeak.get(j).getPeakProperty() == PeakProperty.N_TERM_Y_ION_ONLY ) selectedPeak.get(j).shiftMass(Constants.NTERM_FIX_MOD);
		}//*/
		
		return tags;
	}
	
	private	TagPool	generateDoublyTags(int minTagLength, int minTagLengthPeptideShouldContain, double massTolerance)
	{// by NA
		// selectedPeak array should be sorted before generating Tags
		assert(selectedPeak != null);
		
		ArrayList<Peak> virtualDoublyPeak = new ArrayList<Peak>();
		for(int i=0; i<selectedPeak.size(); i++){
			
			if( selectedPeak.get(i).property == PeakProperty.N_TERM_B_ION_ONLY || 
					selectedPeak.get(i).property == PeakProperty.N_TERM_Y_ION_ONLY ){
				virtualDoublyPeak.add( selectedPeak.get(i) );
				continue;
			}
			
			if( selectedPeak.get(i).property == PeakProperty.C_TERM_B_ION_ONLY || 
					selectedPeak.get(i).property == PeakProperty.C_TERM_Y_ION_ONLY ){
				virtualDoublyPeak.add( selectedPeak.get(i) );
				continue;
			}
			
			if( selectedPeak.get(i).getMass()*2 < correctedMW )
			{
				Peak virtualPeak = new Peak(-1, selectedPeak.get(i).getMass()*2-Constants.Proton, selectedPeak.get(i).getIntensity(), 1, PeakProperty.VIRTUAL_PEAK);
				virtualPeak.setNormIntensity( selectedPeak.get(i).getNormIntensity() );				
				virtualDoublyPeak.add( virtualPeak );
			}
		}
		Collections.sort( virtualDoublyPeak );
		setScoreOfSelectedPeaks(virtualDoublyPeak, 2, massTolerance);
		
	//	double maxGap = AminoAcid.getAminoAcid('W').getMass() + massTolerance; // Mass of 'W' is biggest
		double maxGap = AminoAcid.getMaxMonoMass() + massTolerance; // Mass of 'W' is biggest
		TagPool tags = new TagPool();		
		for(int i=0; i<virtualDoublyPeak.size()-1; i++)
		{
			for(int j=i+1; j<virtualDoublyPeak.size(); j++){
				double diff = Peak.getMassDifference(virtualDoublyPeak.get(i), virtualDoublyPeak.get(j));
				if( diff > maxGap )
					break;
				
				// one length tag
				ArrayList<AminoAcid> codeList = AminoAcid.getCode(diff, massTolerance);
				for(AminoAcid code : codeList){//public Peak(int index, double mass, double intensity, int charge, PeakProperty property)				
					tags.add(new Tag(virtualDoublyPeak.get(i), virtualDoublyPeak.get(j), code, this));
				}							
			}			
		}
		
		TagPool oneLengthTagPool 	= (TagPool)tags.clone();
		TagPool tempTags			= (TagPool)tags.clone();
		
		for(int i=1; i<minTagLengthPeptideShouldContain; i++) {
			tempTags = TagPool.extendTags(tempTags, oneLengthTagPool);		
			tags.addAll(tempTags);
		}	

		return tags;
	}
	
	public Tag getB2Tag(Sequence seq)
	{
		assert(selectedPeak != null);
		
		int virB0 = -1;
		for(int j=0; j<selectedPeak.size(); j++){
			if( selectedPeak.get(j).getPeakProperty() == PeakProperty.N_TERM_B_ION_ONLY ){
				virB0 = j;
				break;
			}
		}
		
		double b2mz = selectedPeak.get(virB0).getMass() + seq.getMonoMass();
		
		for(int j=virB0; j<selectedPeak.size(); j++){
			
			double diff = b2mz - selectedPeak.get(j).getMass();			
			if( diff < -Constants.fragmentTolerance ) break;
			
			if( Math.abs(diff) <= Constants.fragmentTolerance ){
			//	System.out.println(name + " | " + b2mz);
				Peak virtualPeak = new Peak(-1, b2mz - seq.getMonoMass(1, 2), 0, 1, PeakProperty.VIRTUAL_PEAK);
				virtualPeak.setNormIntensity(0);
				virtualPeak.setProbability(0);
				Tag t1 = new Tag(selectedPeak.get(virB0), virtualPeak, seq.get(0), this);
				Tag t2 = new Tag(virtualPeak, selectedPeak.get(j), seq.get(1), this);				
				return Tag.merge(t1, t2);	
			}
		}//*/
		return null;
	}
	
	public int binarySearchforPeaks( double left )
	{
		int index;	
		if( left <= this.get(0).getMass() )
			index= 0;
		else if( left > this.get(this.size()-1).getMass() )
			index= this.size()-1;
		else
		{
			int M, L= 0, R= this.size()-1;
			while( R - L > 1 )
			{
				M= ( L + R ) /2;

				if( left <= this.get(M).getMass() )
					R= M;
				else
					L= M;
			}
			index= R;
		}	
		return index;
	}
	
	public double getMatchedPeak( double mz ){
		double it=0;

		int id= binarySearchforPeaks( mz- Constants.fragmentTolerance );
		if( this.get(id).getMass() < mz - Constants.fragmentTolerance )
			return it;
		
		while( this.get(id).getMass() <= mz + Constants.fragmentTolerance )
		{
			if( this.get(id).getNormIntensity() > it )
				it = this.get(id).getNormIntensity();		
			id++;
			if( id == this.size() )
				break;
		}
		return it;
	}

	public boolean isConfidentPeak( Peak tar, double factor ){
	//	System.out.println(tar.mass +" "+ tar.getIndex() + " " + this.size());
		int me = tar.getIndex();
		if( me < 0 ) return true;
		if( me > 0  && tar.getMass() - this.get(me-1).getMass() < 1.5 ){
			if( tar.getIntensity()/5 < this.get(me-1).getIntensity() ) return false;
		}

		if( me < this.size()-1 && this.get(me+1).getMass() - tar.getMass() < 1.5 ){
			if( tar.getIntensity() < this.get(me+1).getIntensity() ) return false;
		}			
			
		return true;		
	}

	
	public boolean checkSorted()	// identical to Tag's method
	{
		Peak tmp = null;
		for (Peak p : this)
		{
			if (tmp != null)
			{
				if (p.compareTo(tmp)<0) return false;
			}
			tmp = p;
		}
		return true;
	}
	
	public int compareTo(Spectrum s) 
	{
		if( observedMW > s.observedMW ) return 1;
		else if( observedMW == s.observedMW ) return 0;
		else return -1;
	}
	
	private void setScoreOfSelectedPeaks( ArrayList<Peak> selected, int assumedCS, double tolerance )
	{
		double isotopeDelta = Constants.IsotopeSpace, NH3Delta = Constants.NH3, H2ODelta = Constants.H2O;
		for( int node=0; node<selected.size(); node++){
			
			if( selected.get(node).property == PeakProperty.N_TERM_B_ION_ONLY || 
					selected.get(node).property == PeakProperty.N_TERM_Y_ION_ONLY ){
				selected.get(node).setProbability(1.);
				continue;
			}
			
			if( selected.get(node).property == PeakProperty.C_TERM_B_ION_ONLY || 
					selected.get(node).property == PeakProperty.C_TERM_Y_ION_ONLY ){
				selected.get(node).setProbability(1.);
				continue;
			}
			
			double targetmz=0;		
			double score = selected.get(node).getNormIntensity(); 
			double ionMZ= selected.get(node).getMass(); 
			double ionIT= selected.get(node).getIntensity();
			
			// check this is isotope???
			int OK = -1;    	
			double okmax = 0;
			targetmz= ionMZ-isotopeDelta;
			for(int i=node-1; i>-1; i--){
				if( selected.get(i).getMass() < targetmz-tolerance ) break;			
				else if( Math.abs(selected.get(i).getMass()-targetmz) < tolerance ){
					if( selected.get(i).getIntensity() > okmax ) {
						okmax = selected.get(i).getIntensity();
						OK = i;
					}
				}
			}
			if( OK != -1 ) {
				if( okmax > ionIT ) selected.get(node).setProbability(0.);
				else selected.get(node).setProbability(score);
				continue;
			}//*/
			

			int ISO = node + 1; // plus isotope peak    
			double prevISO = ionIT;
			boolean isoDecent = false;
			targetmz= ionMZ + isotopeDelta;
			for( int nth_iso=1 ;  ; nth_iso++  ){
				int H = -1;
				double imax = 0;				
				for(int i=ISO; i<selected.size(); i++){
					if( selected.get(i).getMass() > targetmz+Constants.massToleranceForDenovo ) break;			
					else if( Math.abs(selected.get(i).getMass()-targetmz) < Constants.massToleranceForDenovo ){
						if( selected.get(i).getIntensity() > imax ) {
							imax = selected.get(i).getIntensity();
							H = i;
						}
					}
				}

				if( H == -1 || prevISO < imax/5 ) break;
				if( isoDecent && prevISO < imax ) break;
				if( prevISO > imax ) isoDecent=true;
				
				score += selected.get(H).getNormIntensity()/nth_iso;	
				ISO = H+1;
				targetmz= selected.get(H).getMass() + isotopeDelta;		
				prevISO = imax;
			}//*/

			int NLOSS = -1;    	
			double lossmax=0, lossmz=0;
			targetmz= ionMZ-H2ODelta;
			for(int i=node-1; i>-1; i--){
				if( selected.get(i).getMass() < targetmz-tolerance ) break;			
				else if( Math.abs(selected.get(i).getMass()-(ionMZ-NH3Delta)) < tolerance ){
					if( selected.get(i).getIntensity() > lossmax ) {
						lossmax = selected.get(i).getIntensity();
						lossmz = selected.get(i).getMass();
						NLOSS = i;
					}
				}
				else if( Math.abs(selected.get(i).getMass()-targetmz) < tolerance ){
					if( selected.get(i).getIntensity() > lossmax ) {
						lossmax = selected.get(i).getIntensity();
						lossmz = selected.get(i).getMass();
						NLOSS = i;
					}
				}
			}
			if( NLOSS != -1 ){
				double lossScore = selected.get(NLOSS).getNormIntensity();				
				int NL_H = -1; // plus isotope peak    	
				double nsmax = 0;
				targetmz= lossmz+isotopeDelta;
				for(int i=NLOSS+1; i<selected.size(); i++){
					if( selected.get(i).getMass() > targetmz+tolerance ) break;			
					else if( Math.abs(selected.get(i).getMass()-targetmz) < tolerance ){
						if( selected.get(i).getIntensity() > nsmax ) {
							nsmax = selected.get(i).getIntensity();
							NL_H = i;
						}
					}
				}
				if( NL_H != -1 ) lossScore += selected.get(NL_H).getNormIntensity();				
				score += lossScore*0.5;
				
			}//*/
		
			double complement = correctedMW - ionMZ + Constants.Proton*2;	
			score += getMatchedPeak( complement ); 
			if( assumedCS == 1 && charge > 2 ) {
				score += getMatchedPeak( (complement+Constants.Proton)/2 ); 
			}			
			selected.get(node).setProbability(score);
		}
	}
}

