package scaniter;

import java.util.ArrayList;
import java.util.Collections;

import moda.ThreadPoolManager;
import modi.Constants;
import modi.Peak;
import modi.Spectrum;
import msutil.MSMass;

public class MSMScan {
	
	static int minPeaksCount = 4;
	static double minMW = 8 * MSMass.getMinAAMass() + 18;
	
	private String 		title;
	private int 		specIndex;
	private int 		scanNo;
	private double 		pmz;
	private double 		neutralMW;
	private int 		charge;
	private long 		offset;
	private Spectrum 	peaklist;
	private static double tolerance= Constants.massToleranceForDenovo;

	private double 	precursorTolerance = 0;
	private double 	precursorAccuracy= 0;
	private double 	gapTolerance = 0;
	private double 	nonModifiedDelta = 0;
	private int		maxNoOfC13 = 0;
	
	public MSMScan(int index, double pmz, int charge){
		this.title 		= "";
		this.specIndex	= index;
		this.scanNo		= 0;	
		this.pmz 		= pmz;
		this.charge 	= charge;
		this.neutralMW 	= (pmz - Constants.Proton)*charge;
	}
	
	public MSMScan(String title, int index, int sn, double pmz, int charge){		
		this.title 		= title;
		this.specIndex  = index;
		this.scanNo		= sn;
		this.pmz 		= pmz;
		this.charge 	= charge;
		this.neutralMW 	= (pmz - Constants.Proton)*charge;
	}

	public double 	getObservedMW(){ return neutralMW; }
	public int 		getCharge(){ return charge; }
	public long 	getOffset(){ return offset; }	
	public String 	getHeader(){ return String.format("%d\t%.4f\t%d\t%d\t%s",
													specIndex, neutralMW, charge, scanNo, title); }
	
	public Spectrum getSpectrum() {
		int slotIdx = ThreadPoolManager.getSlotIndex();
		Constants.precursorTolerance[slotIdx]= precursorTolerance;
		Constants.precursorAccuracy[slotIdx]	= precursorAccuracy;
		Constants.gapTolerance[slotIdx] 		= gapTolerance;
		Constants.gapAccuracy[slotIdx] 		= precursorAccuracy + 2*Constants.fragmentTolerance;
		Constants.nonModifiedDelta[slotIdx] 	= nonModifiedDelta;
		Constants.maxNoOfC13[slotIdx] 		= maxNoOfC13;
		return peaklist; 
	}

	public boolean setSpectrum(ArrayList<RawPeak> rawPL) {
		int slotIdx = ThreadPoolManager.getSlotIndex();
		if( neutralMW < minMW || Constants.maxPeptideMass < neutralMW ) return false; 
				
		if( Constants.reporterMassOfIsobaricTag != null ) removeReporterIons(rawPL, Constants.reporterMassOfIsobaricTag);
		
		if( Constants.rangeForIsotopeIncrement != 0 ) maxNoOfC13 = (int)Math.ceil( neutralMW / Constants.rangeForIsotopeIncrement );
		else maxNoOfC13 = Constants.maxNoOfC13[slotIdx];
		
		if( Constants.PPMTolerance != 0 ) precursorAccuracy = Constants.PPMtoDalton( neutralMW, Constants.PPMTolerance );
		else
			precursorAccuracy = Constants.precursorAccuracy[slotIdx];
		
		precursorTolerance = precursorAccuracy + maxNoOfC13*Constants.IsotopeSpace;
		
		int index = 0;
		Spectrum spectrum = new Spectrum( this.pmz, this.charge, this.title );
		
		double basePeakIntensity=0, TIC=0;
		double tarMass=0, tarInten=0;
		for( RawPeak rp : rawPL ) {		
			double mass = rp.mz;
			double intensity = rp.it;
			if( intensity <= 0 || mass <= 0 ) continue;
			if( mass > neutralMW ) continue;
			
			if( ( mass - tarMass ) < tolerance ){
				double sum = tarInten + intensity;
				tarMass = tarMass*(tarInten/sum)+ mass*(intensity/sum);
				tarInten += intensity;
				spectrum.get(index-1).set(tarMass, tarInten);
			}
			else{
				spectrum.add( new Peak(index++, mass, intensity) );
				tarMass = mass;
				tarInten = intensity; 
			}
			TIC += intensity;
			if( tarInten > basePeakIntensity )
				basePeakIntensity= tarInten;
		}	
		spectrum.setExtraInformation( basePeakIntensity, TIC );
		
		gapTolerance = Constants.fragmentTolerance*2;
		nonModifiedDelta = (precursorTolerance < Constants.massToleranceForDenovo)? precursorTolerance : Constants.massToleranceForDenovo;
				
		if( precursorTolerance > gapTolerance ) gapTolerance += precursorTolerance;
		
		if( spectrum.size() < minPeaksCount ) peaklist = null; 
		else peaklist = spectrum;
		
		return (peaklist!=null);
	}
	
	private void removeReporterIons( ArrayList<RawPeak> rawPL, double[] removedMasses ){
	
		ArrayList<RawPeak> reporters = new ArrayList<RawPeak>();
		for(int i=1; i<removedMasses.length; i++)
			reporters.add( new RawPeak(removedMasses[i], Constants.fragmentTolerance) );
		
		reporters.add( new RawPeak(removedMasses[0] + Constants.Proton, Constants.fragmentTolerance) );

		int fragCS = 1;
		while( true ){
			double compItraqTag = (this.neutralMW - removedMasses[0] + Constants.Proton*fragCS)/fragCS;
	
			double secondIso = compItraqTag + Constants.IsotopeSpace/fragCS;
			double thirdIso  = secondIso + Constants.IsotopeSpace/fragCS;
			double forthIso  = thirdIso + Constants.IsotopeSpace/fragCS;
			
			reporters.add( new RawPeak(compItraqTag, Constants.fragmentTolerance) );
			reporters.add( new RawPeak(secondIso, Constants.fragmentTolerance) );
			reporters.add( new RawPeak(thirdIso, Constants.fragmentTolerance) );
			reporters.add( new RawPeak(forthIso, Constants.fragmentTolerance) );

			for(int i=1; i<=maxNoOfC13; i++){
				reporters.add( new RawPeak(compItraqTag-i*Constants.IsotopeSpace/fragCS, Constants.fragmentTolerance) );
			}
			if( ++fragCS >= this.charge ) break;
		}
		
		Collections.sort(reporters);
		
		int start = 0;
		for( RawPeak rp : reporters ){
			for (int i=start; i<rawPL.size(); i++){
				if( rawPL.get(i).mz < rp.mz-rp.it ) continue;
				else if( rawPL.get(i).mz > rp.mz+rp.it ) {
					start = i;
					break;
				}			
				rawPL.remove(i);
				i--;
			}
		}	
	}
}

























