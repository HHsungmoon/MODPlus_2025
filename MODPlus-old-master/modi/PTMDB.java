package modi;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

import msutil.MSMass;

class GapKey extends Object {
	private Sequence seq;
	private int massDiff;
	private static int currentID = 0;
	private static final int unitSize = 10;
	private PTMPosition pos;
	
	public GapKey( Sequence seq, double massDiff, PTMPosition pos )
	{
		this.seq = seq;
		this.massDiff = (int)(massDiff*unitSize); // should be considered later. 4 just means a big enough tolerance for 0.5
		this.pos = pos;
	}
	
	public String	getString()		{ return seq.toString(); }
	public Sequence	getSeq()		{ return seq; }
	public double	getMassDiff()	{ return massDiff/(double)unitSize; }
	public int hashCode()
	{
		return seq.hashCode() + massDiff + pos.hashCode() ;
	}
	
	public boolean equals(Object key)
	{
		if (key instanceof GapKey)
			return this.seq.equals( ((GapKey)key).seq ) && this.massDiff==((GapKey)key).massDiff 
				&& this.pos==((GapKey)key).pos;
		else return false;
	}
}


public class PTMDB extends ArrayList<PTM> {
	private ArrayList<PTM>[][]		PTMTable;
	private TreeMap<Integer, String> classifications = new TreeMap<Integer, String>();

	public boolean isNovelPtm(int site, int res, double mass){
		int i;
		for(i=0; i<PTMTable[res][PTMPosition.ANYWHERE.ordinal()].size(); i++ ){
			int delta= (int)Math.round(PTMTable[res][PTMPosition.ANYWHERE.ordinal()].get(i).getMassDifference());
			if( delta == (int)Math.round(mass) )
				return false;		
		}		
		if( site == 0 ){
			for(i=0; i<PTMTable[res][PTMPosition.ANY_N_TERM.ordinal()].size(); i++ ){
				int delta= (int)Math.round(PTMTable[res][PTMPosition.ANY_N_TERM.ordinal()].get(i).getMassDifference());
				if( delta == (int)Math.round(mass) )
					return false;	
			}
		}
		return true;	
	}
	
	public void constructPTMLookupTable(){
		constructPTMTable();
	}

	private int constructPTMTable() {	
		this.sortByPTMPosition();
				
		PTMTable = new ArrayList[AminoAcid.getIndexSize()][PTMPosition.PTMPOSITION_COUNT.ordinal()];
		for (int i=0; i<AminoAcid.getIndexSize(); i++) {
			for (int j=0; j<PTMPosition.PTMPOSITION_COUNT.ordinal(); j++)
				PTMTable[i][j] = new ArrayList<PTM>();		
		}
		
		
		for ( PTM ptm : this ) {
			if( ptm.getMassDifference() < Constants.minModifiedMass || ptm.getMassDifference() > Constants.maxModifiedMass ) continue;
			if ( ptm.getResidue() != null ) {
				int aa = ptm.getResidue().getIndex();
				PTMTable[ptm.getResidue().getIndex()][ptm.getPTMPosition().ordinal()].add(ptm);
				
				if( ptm.getPTMPosition() == PTMPosition.ANYWHERE ){
					if( modRange[0][aa] > ptm.getMassDifference() ) modRange[0][aa] = ptm.getMassDifference(); // Min
					else if( modRange[1][aa] < ptm.getMassDifference() ) modRange[1][aa] = ptm.getMassDifference(); // Max
				}
				else if( ptm.getPTMPosition() == PTMPosition.ANY_N_TERM || ptm.getPTMPosition() == PTMPosition.PROTEIN_N_TERM ){
					if( ntermModRange[0][aa] > ptm.getMassDifference() ) ntermModRange[0][aa] = ptm.getMassDifference(); // Min
					else if( ntermModRange[1][aa] < ptm.getMassDifference() ) ntermModRange[1][aa] = ptm.getMassDifference(); // Max
				}
				else if( ptm.getPTMPosition() == PTMPosition.ANY_C_TERM || ptm.getPTMPosition() == PTMPosition.PROTEIN_C_TERM ){
					if( ctermModRange[0][aa] > ptm.getMassDifference() ) ctermModRange[0][aa] = ptm.getMassDifference(); // Min
					else if( ctermModRange[1][aa] < ptm.getMassDifference() ) ctermModRange[1][aa] = ptm.getMassDifference(); // Max
				}
			}
			else {
				if ( !ptm.isResidueCTerm() && !ptm.isResidueNTerm() ) return -1;
				
				for (int j=0; j<AminoAcid.getIndexSize(); j++){				
					
					if( PTMTable[j][PTMPosition.ANYWHERE.ordinal()].contains(ptm) ) continue;
					
					if ( ptm.getPTMPosition() == PTMPosition.PROTEIN_N_TERM ) {
						if( PTMTable[j][PTMPosition.ANY_N_TERM.ordinal()].contains(ptm) )
							continue;
					}					
					if ( ptm.getPTMPosition() == PTMPosition.PROTEIN_C_TERM ) {
						if( PTMTable[j][PTMPosition.ANY_C_TERM.ordinal()].contains(ptm) )
							continue;
					}//*/

					PTMTable[j][ptm.getPTMPosition().ordinal()].add(ptm);
					
					if( ptm.getPTMPosition() == PTMPosition.ANY_N_TERM || ptm.getPTMPosition() == PTMPosition.PROTEIN_N_TERM ){
						if( ntermModRange[0][j] > ptm.getMassDifference() ) ntermModRange[0][j] = ptm.getMassDifference(); // Min
						else if( ntermModRange[1][j] < ptm.getMassDifference() ) ntermModRange[1][j] = ptm.getMassDifference(); // Max
					}
					else if( ptm.getPTMPosition() == PTMPosition.ANY_C_TERM || ptm.getPTMPosition() == PTMPosition.PROTEIN_C_TERM ){
						if( ctermModRange[0][j] > ptm.getMassDifference() ) ctermModRange[0][j] = ptm.getMassDifference(); // Min
						else if( ctermModRange[1][j] < ptm.getMassDifference() ) ctermModRange[1][j] = ptm.getMassDifference(); // Max
					}
				}
			}
		}
		
		Collections.sort(this);	
		return 0;
	}
	
	// Temporary valuables for recursive search
	private Sequence			seq;
	private	double				massDiff;
	private	PTMPosition			position;
	private ArrayList<PTMRun>	result;
	private PTM[]				occur;
	private int[]				numNextFixSite;
	private int					numMaxMods;

	// Finding PTMs by Exhaustive Depth-first-search
	private void findPTM_DFS( double mass, int pos, int cnt, int extra )
	{	
		// when search has reached the leaf node(end of sequence), check result
		if ( pos == seq.size() ) {			
			double ierror = Math.abs( mass - massDiff );
		//	System.out.println(Constants.gapTolerance);
		//	if ( ierror <= Constants.gapTolerance ) 
			if ( ierror <= Constants.gapTolerance && Constants.isWithinAccuracy(ierror) ) 
			{
				PTMRun run = new PTMRun();
				for (int i=0; i<seq.size(); i++)
					if ( occur[i] != null ){
						run.add(new PTMOccurrence(i, occur[i]));
					}
				
				if ( run.size() > 0 ) {
					run.setError( ierror );
					this.result.add(run);
				}
			}
			return;
		}
		
		// recursion call without checking PTM on this pos
		findPTM_DFS( mass, pos+1, cnt, extra );		
		
		if( cnt >= Constants.maxPTMPerGap ) return;
		if( extra == 1 && numMaxMods == 1 ) {
			if ( numNextFixSite[pos] == 0 ) return;
		}

		int residueIndex = seq.get(pos).getIndex();
		for (int i=1; i<PTMPosition.PTMPOSITION_COUNT.ordinal(); i++) {
			if (pos>0 && pos<seq.size()-1 && i!=1) continue;
			if ((i==2 || i==4) && pos!=0) continue;
			if ((i==3 || i==5) && pos!=seq.size()-1) continue;
			
			if ((i==2) && this.position!=PTMPosition.ANY_N_TERM 
					&& this.position!=PTMPosition.PROTEIN_N_TERM ) continue;
			if ((i==3) && this.position!=PTMPosition.ANY_C_TERM 
					&& this.position!=PTMPosition.PROTEIN_C_TERM ) continue;
			
			if ((i==4) && this.position!=PTMPosition.PROTEIN_N_TERM ) continue;
			if ((i==5) && this.position!=PTMPosition.PROTEIN_C_TERM ) continue;
			
			for (PTM ptm : PTMTable[residueIndex][i]) {				
				occur[pos] = ptm;				
				if( extra + ptm.getModCount() > numMaxMods ) continue;
				findPTM_DFS( mass + ptm.getMassDifference(), pos+1, cnt+1, extra+ptm.getModCount() );
						
			}
		}
		occur[pos] = null;
	}

	// Hash table for reusing the result for Gaps already searched once
	private HashMap< GapKey, PTMSearchResult >	hashTable = new HashMap< GapKey, PTMSearchResult >();
	
	public void clearHashTable()		{ hashTable.clear(); }
	
	// before ver 0.8, deprecated, Not Cashing
	public PTMSearchResult searchPTM( Sequence seq, double massDiff, PTMPosition position ) {		
		PTMSearchResult searchResult;
		// Hash miss
		ArrayList<PTMRun> newGapInterpret = new ArrayList<PTMRun>();
	
		double ierror = Math.abs( massDiff );
		if( ierror < Constants.nonModifiedDelta ){
			searchResult = new PTMSearchResult( newGapInterpret, true );			
			return searchResult;
		}
		
		if( ierror < Constants.gapTolerance ){
			PTMRun run = new PTMRun();
			run.setError(ierror);
			newGapInterpret.add( run );			
		}
		
		// PTM search
		this.seq 		= seq;
		this.massDiff	= massDiff;
		this.position	= position;
		this.result		= newGapInterpret;
		occur = new PTM[seq.size()];
		numNextFixSite = new int[seq.size()];				
		numMaxMods = Constants.getMaxPTMOccurrence(seq.size());
				
		int cum = ( Constants.CTERM_FIX_MOD !=0 && (this.position == PTMPosition.ANY_C_TERM || this.position == PTMPosition.PROTEIN_C_TERM) )? 1 : 0;
		for(int i=seq.size()-1; i>-1; i--){
			if( seq.get(i).isLabelled() ) {
				numNextFixSite[i] = ++cum;
			}
			else {				
				numNextFixSite[i] = cum;
			}
		}
		
		findPTM_DFS( 0.0, 0, 0, 0 );
		this.result		= null;		// release reference from PTMDB class
		
		searchResult = new PTMSearchResult(newGapInterpret, newGapInterpret != null && newGapInterpret.size() > 0);
		
	//	hashTable.put( new GapKey(seq, massDiff, position), searchResult );
		
		return searchResult;
	}
	
	public org.jdom.Element getGapInterpretListElement()
	{
		if(Constants.ANALYSIS_VERSION < 1.0)
			return null;
		
		org.jdom.Element gapInterpretListElement = new org.jdom.Element("gapInterpretList");
		// HashMap< GapKey, ArrayList<PTMRun> > hashTable
		Iterator<Map.Entry<GapKey, PTMSearchResult>> it = hashTable.entrySet().iterator();
		
		// to print actualInterpret entry by order of id
		
		HashMap<Integer, GapKey> gapKeyOrder = new HashMap<Integer, GapKey>();
		while(it.hasNext())
		{
			Map.Entry<GapKey, PTMSearchResult> entry = it.next();
			gapKeyOrder.put(entry.getValue().getID(), entry.getKey());
		}

		ArrayList<Integer> orderedIDSet = new ArrayList<Integer> (gapKeyOrder.keySet());
		Collections.sort(orderedIDSet);
			
		for(Integer keyID : orderedIDSet)
		{
			GapKey key = gapKeyOrder.get(keyID);
			PTMSearchResult result = hashTable.get(key);
			
			org.jdom.Element elemGapInterpret = new org.jdom.Element("gapInterpret");
			elemGapInterpret.setAttribute(new org.jdom.Attribute("id", String.valueOf(result.getID())));
			elemGapInterpret.setAttribute(new org.jdom.Attribute("sequence", key.getString()));
			elemGapInterpret.setAttribute(new org.jdom.Attribute("offset", String.valueOf(key.getMassDiff())));
			String hasPTMStr;
			ArrayList<PTMRun> ptmRunList = result.getPTMRun();
			if(ptmRunList == null || ptmRunList.isEmpty())
				hasPTMStr = "no";
			else
				hasPTMStr = "yes";
			elemGapInterpret.setAttribute(new org.jdom.Attribute("hasPTM", hasPTMStr));
			int id = 0;
			for(PTMRun ptmRun : ptmRunList)
			{
				org.jdom.Element elemActualInterpret = new org.jdom.Element("actualInterpret");
				elemActualInterpret.setAttribute(new org.jdom.Attribute("id", String.valueOf(id++)));
				
				org.jdom.Element elemPTMRun = new org.jdom.Element("PTMRun");
				for(PTMOccurrence occr : ptmRun)
				{
					org.jdom.Element elemPTMOccr = new org.jdom.Element("PTMOccr");
					elemPTMOccr.setAttribute(new org.jdom.Attribute("id", String.valueOf(occr.getPTM().getID())));
					elemPTMOccr.setAttribute(new org.jdom.Attribute("position", String.valueOf(occr.getPosition())));
					elemPTMRun.addContent(elemPTMOccr);
				}
				elemActualInterpret.addContent(elemPTMRun);
				
				org.jdom.Element elemBTheoPeaks = new org.jdom.Element("bTheoPeaks");
				elemBTheoPeaks.setText(Gap.getTheoreticalMassStr(getBTheoreticalPeaks(key.getSeq(), ptmRun)));
				elemActualInterpret.addContent(elemBTheoPeaks);

				org.jdom.Element elemYTheoPeaks = new org.jdom.Element("yTheoPeaks");
				elemYTheoPeaks.setText(Gap.getTheoreticalMassStr(getYTheoreticalPeaks(key.getSeq(), ptmRun)));
				elemActualInterpret.addContent(elemYTheoPeaks);
				
				elemGapInterpret.addContent(elemActualInterpret);
			}
			gapInterpretListElement.addContent(elemGapInterpret);
		}
		return gapInterpretListElement;
	}

	public	ArrayList<Peak> getBTheoreticalPeaks(Sequence seq, PTMRun run)
	{
		double [] ptmMass = new double[seq.size()];
		for(PTMOccurrence occr : run)
			ptmMass[occr.getPosition()] += occr.getPTM().getMassDifference();

		ArrayList<Peak> bPeaks = new ArrayList<Peak>();
		double mass = 0;
		for(int i=0; i<seq.size()-1; i++)
		{
			mass = mass + seq.get(i).getMass() + ptmMass[i];
			bPeaks.add(new Peak(-1, mass, 0.));
		}
		return bPeaks;
	}
	
	public	ArrayList<Peak> getYTheoreticalPeaks(Sequence seq, PTMRun run)
	{
		double [] ptmMass = new double[seq.size()];
		for(PTMOccurrence occr : run)
			ptmMass[occr.getPosition()] += occr.getPTM().getMassDifference();

		ArrayList<Peak> yPeaks = new ArrayList<Peak>();
		double mass = 0;
		for(int i=0; i<seq.size()-1; i++)
		{
			mass = mass + seq.get(seq.size()-1-i).getMass() + ptmMass[i];
			yPeaks.add(new Peak(-1, mass, 0.));
		}
		return yPeaks;
	}

	public PTM	getPTM(String name, char residue){
		int AA = AminoAcid.getAminoAcid(residue).getIndex();
		for(int j=0; j<PTMPosition.PTMPOSITION_COUNT.ordinal(); j++){
			for ( PTM ptm : PTMTable[AA][j] ){
				if( ptm.getName().equals(name) ){
					return ptm;
				}
			}
		}
		return null;
	}
	
	public void sortByPTMPosition(){
		Collections.sort( this, new PTMPosComparator() );
	}

	
	public int setVariableModificatinos( String fileName, double[] fixedAA, boolean canBeModifiedOnFixedAA ) throws Exception
	{
		int id = this.size();
		org.jdom.Document document = null;
		try {
			 document = new org.jdom.input.SAXBuilder().build(new File(fileName));
		} catch(org.jdom.JDOMException e) {
			System.out.println(e);
		}
		
		org.jdom.Element rootElement = document.getRootElement();
		org.jdom.Element classificationsElement = rootElement.getChild("classifications");
		if(classificationsElement != null)	// for backward compatibility
		{
			classifications.clear();
			for(Object obj : classificationsElement.getChildren())
			{
				org.jdom.Element classificationRowElement = (org.jdom.Element)obj;
				classifications.put(Integer.parseInt(classificationRowElement.getAttributeValue("recordID")),
						classificationRowElement.getAttributeValue("classification"));
			}
		}		

		for(Object obj : rootElement.getChildren("PTM"))
		{
			org.jdom.Element elemPTM = (org.jdom.Element)obj;

			String pname = elemPTM.getChildText("name");
			
			String category = elemPTM.getChildText("classification");
			
			String residueStr = elemPTM.getChildText("residue");
			AminoAcid residue = null;
			
			double massdelta = 0;
			if(residueStr != null && residueStr.compareToIgnoreCase("N-term") != 0 && residueStr.compareToIgnoreCase("C-term") != 0) {
				assert( residueStr.length() == 1 );
				residue = AminoAcid.getAminoAcid(residueStr.charAt(0));
				if( fixedAA[residueStr.charAt(0)-'A'] != 0 ) {
					if( canBeModifiedOnFixedAA ) massdelta = fixedAA[residueStr.charAt(0)-'A'];
					else continue;
				}
				if( "AA substitution".compareTo(category) == 0 ){
					char tarAA = AminoAcid.getAminoAcidBy3Letter(pname.substring(pname.length()-3));
					massdelta -= fixedAA[tarAA-'A'];
					if( tarAA == 'I' && fixedAA['I'-'A'] != fixedAA['L'-'A'] ){
						System.out.println("PROGRAM SHOULD BE FIXED");
						System.exit(1);
					}
				}
			}		
			
			PTMPosition position = null;
			String positionStr = elemPTM.getChildText("position");
			if(positionStr.equalsIgnoreCase("ANYWHERE"))
				position = PTMPosition.ANYWHERE;
			else if(positionStr.equalsIgnoreCase("ANY_N_TERM"))
				position = PTMPosition.ANY_N_TERM;
			else if(positionStr.equalsIgnoreCase("ANY_C_TERM"))
				position = PTMPosition.ANY_C_TERM;
			else if(positionStr.equalsIgnoreCase("PROTEIN_N_TERM"))
				position = PTMPosition.PROTEIN_N_TERM;
			else if(positionStr.equalsIgnoreCase("PROTEIN_C_TERM"))
				position = PTMPosition.PROTEIN_C_TERM;
			else
				assert(false);

		//	if( Constants.NTERM_FIX_MOD != 0 && (position == PTMPosition.ANY_N_TERM || position == PTMPosition.PROTEIN_N_TERM) ) continue;
		//	if( Constants.CTERM_FIX_MOD != 0 && (position == PTMPosition.ANY_C_TERM || position == PTMPosition.PROTEIN_C_TERM) ) continue;

			if( Constants.NTERM_FIX_MOD != 0 && (position == PTMPosition.ANY_N_TERM || position == PTMPosition.PROTEIN_N_TERM) ) {
				if( canBeModifiedOnFixedAA ) {
					pname += "/Nterm";
					massdelta += Constants.NTERM_FIX_MOD;
				}
				else continue;
			}
			if( Constants.CTERM_FIX_MOD != 0 && (position == PTMPosition.ANY_C_TERM || position == PTMPosition.PROTEIN_C_TERM) ) {
				if( canBeModifiedOnFixedAA ) {
					massdelta += Constants.CTERM_FIX_MOD;
					pname += "/Cterm";
				}
				else continue;
			}
			
			double ac_delta = Double.parseDouble(elemPTM.getChildText("massDifference"));
			if( ac_delta < Constants.minModifiedMass || ac_delta > Constants.maxModifiedMass ) continue;
			ac_delta -= massdelta;
			
			this.add(new PTM(
						id++,
						pname, 
						"", 
						ac_delta,
						0,
						residue,
						position,
						category
						)
					);
		}
		return 1;
	}
	
	public int setVariableModificatinos( org.jdom.Element modifications, double[] fixedAA, boolean canBeModifiedOnFixedAA ) throws Exception
	{
		int id = this.size();
		for(Object obj : modifications.getChildren("mod"))
		{
			org.jdom.Element elemPTM = (org.jdom.Element)obj;
			
			String residueStr = elemPTM.getAttributeValue("site");
			AminoAcid residue = null;
			double massdelta = 0;
			if(residueStr != null && residueStr.compareToIgnoreCase("N-term") != 0 && residueStr.compareToIgnoreCase("C-term") != 0) {
				assert(residueStr.length() == 1);
				residue = AminoAcid.getAminoAcid(residueStr.charAt(0));
				if( fixedAA[residueStr.charAt(0)-'A'] != 0 ) {
					if( canBeModifiedOnFixedAA ) massdelta = fixedAA[residueStr.charAt(0)-'A'];
					else continue;
				}
			}			
			
			PTMPosition position = null;
			String positionStr = elemPTM.getAttributeValue("position");
			if(positionStr.equalsIgnoreCase("ANYWHERE"))
				position = PTMPosition.ANYWHERE;
			else if(positionStr.equalsIgnoreCase("ANY_N_TERM"))
				position = PTMPosition.ANY_N_TERM;
			else if(positionStr.equalsIgnoreCase("ANY_C_TERM"))
				position = PTMPosition.ANY_C_TERM;
			else if(positionStr.equalsIgnoreCase("PROTEIN_N_TERM"))
				position = PTMPosition.PROTEIN_N_TERM;
			else if(positionStr.equalsIgnoreCase("PROTEIN_C_TERM"))
				position = PTMPosition.PROTEIN_C_TERM;
			else
				assert(false);

			String pname = elemPTM.getAttributeValue("name");
			
			if( Constants.NTERM_FIX_MOD != 0 && (position == PTMPosition.ANY_N_TERM || position == PTMPosition.PROTEIN_N_TERM) ) {
				if( canBeModifiedOnFixedAA ) {
					pname += "/Nterm";
					massdelta += Constants.NTERM_FIX_MOD;
				}
				else continue;
			}
			if( Constants.CTERM_FIX_MOD != 0 && (position == PTMPosition.ANY_C_TERM || position == PTMPosition.PROTEIN_C_TERM) ) {
				if( canBeModifiedOnFixedAA ) {
					pname += "/Cterm";
					massdelta += Constants.CTERM_FIX_MOD;
				}
				else continue;
			}
			
			double ac_delta = Double.parseDouble(elemPTM.getAttributeValue("massdiff"));
			if( ac_delta < Constants.minModifiedMass || ac_delta > Constants.maxModifiedMass ) continue;
			ac_delta -= massdelta;
			
			this.add(new PTM(
						id++,
						pname, 
						"", 
						ac_delta,
						0,
						residue,
						position
						)
					);
		}
		return 1;//constructPTMTable();
	}

	public int setFixedModificatinos(org.jdom.Element modifications, double[] fixedAA) throws Exception
	{	
		int id=0;
		for(Object obj : modifications.getChildren("mod"))
		{
			org.jdom.Element elemPTM = (org.jdom.Element)obj;
			
			String residueStr = elemPTM.getAttributeValue("site");
			AminoAcid residue = null;
			
			if( residueStr == null ) continue;
			else if( residueStr.compareToIgnoreCase("N-term") == 0 ){
				Constants.NTERM_FIX_MOD += Double.parseDouble(elemPTM.getAttributeValue("massdiff"));
			}
			else if( residueStr.compareToIgnoreCase("C-term") == 0 ){
				Constants.CTERM_FIX_MOD += Double.parseDouble(elemPTM.getAttributeValue("massdiff"));
			}
			else residue = AminoAcid.getAminoAcid(residueStr.charAt(0));
			
			double diff = Double.parseDouble(elemPTM.getAttributeValue("massdiff"));		
			
			PTMPosition position = null;
			String positionStr = elemPTM.getAttributeValue("position");
			if(positionStr.equalsIgnoreCase("ANYWHERE"))
				position = PTMPosition.ANYWHERE;
			else if(positionStr.equalsIgnoreCase("ANY_N_TERM"))
				position = PTMPosition.ANY_N_TERM;
			else if(positionStr.equalsIgnoreCase("ANY_C_TERM"))
				position = PTMPosition.ANY_C_TERM;
			else if(positionStr.equalsIgnoreCase("PROTEIN_N_TERM"))
				position = PTMPosition.PROTEIN_N_TERM;
			else if(positionStr.equalsIgnoreCase("PROTEIN_C_TERM"))
				position = PTMPosition.PROTEIN_C_TERM;
			else
				assert(false);

			this.add(new PTM(
						id++,
						elemPTM.getAttributeValue("name"), 
						"", 
						diff,
						diff,
						residue,
						position
						)
					);
			if( residue != null ){
				if( fixedAA[residueStr.charAt(0)-'A'] != 0 ) return 0;
				AminoAcid.modifiedAminoAcidMass( residueStr.charAt(0), diff );
				MSMass.modifiedAminoAcidMass( residueStr.charAt(0), diff );
				fixedAA[residueStr.charAt(0)-'A']= diff;
			}
		}
		return 1;
	}

	private double[][] modRange = new double[2][26];
	private double[][] ntermModRange = new double[2][26];
	private double[][] ctermModRange = new double[2][26];
	
	public double minimumModifiedMass( Sequence seq, int pos ){
		double min1=0, min2=0, extra = 0;
		for( AminoAcid aa : seq ){
			int aix= aa.getIndex();
			if( min1 > modRange[0][aix] ) {
				min2 = min1;
				min1 = modRange[0][aix];
			}
			else if( min2 > modRange[0][aix] ) min2 = modRange[0][aix];
		}
		
		if( pos == 1 ){
			int aix= seq.get(0).getIndex();
			if( min1 > ntermModRange[0][aix] ) {
				min2 = min1;
				min1 = ntermModRange[0][aix];
			}
			else if( min2 > ntermModRange[0][aix] ) min2 = ntermModRange[0][aix];
		}
		else if( pos == 2 ){
			int aix= seq.get(seq.size()-1).getIndex();
			if( min1 > ctermModRange[0][aix] ) {
				min2 = min1;
				min1 = ctermModRange[0][aix];
			}
			else if( min2 > ctermModRange[0][aix] ) min2 = ctermModRange[0][aix];
		}
		if( seq.size() == 1 ) return min1;
		return min1 + min2;
	}
	
	public double maximumModifiedMass(Sequence seq, int pos){
		double max1=0, max2=0;
		for( AminoAcid aa : seq ){
			int aix= aa.getIndex();
			if( max1 < modRange[1][aix] ) {
				max2 = max1;
				max1 = modRange[1][aix];
			}
			else if( max2 < modRange[1][aix] ) max2 = modRange[1][aix];
		}
		if( pos == 1 ){
			int aix= seq.get(0).getIndex();
			if( max1 < ntermModRange[1][aix] ) {
				max2 = max1;
				max1 = ntermModRange[1][aix];
			}
			else if( max2 < ntermModRange[1][aix] ) max2 = ntermModRange[1][aix];
		}
		else if( pos == 2 ){
			int aix= seq.get(seq.size()-1).getIndex();
			if( max1 < ctermModRange[1][aix] ) {
				max2 = max1;
				max1 = ctermModRange[1][aix];
			}
			else if( max2 < ctermModRange[1][aix] ) max2 = ctermModRange[1][aix];
		}
		if( seq.size() == 1 ) return max1;
		return max1 + max2;
	}
	
	public void setPTMDiagnosticIon() {		
		for (PTM ptm : this) {
			
			//setting common modifications
			if( Math.abs(15.994915-ptm.getMassDifference()) < 0.01 ){ 
				if( ptm.getAbbAA() =='M' ) {
					ptm.setPenalty(0.5);
					ptm.setModCount(0);
					AminoAcid.canBeEasilyModified('M');
					ptm.setNeutralLoss( ptm.getMassDifference()+48.003371 );
				}
				if( ptm.getAbbAA() =='W' ) ptm.setPenalty(0.8);	
			}
			
			if( ptm.getAbbAA() =='Q' && ptm.getPTMPosition() == PTMPosition.ANY_N_TERM && Constants.round(ptm.getMassDifference()+Constants.NTERM_FIX_MOD) == -17 ){
				ptm.setPenalty(0.5);
				ptm.setModCount(0);
			}

			//setting isobaric labeling 
			if( Constants.reporterMassOfIsobaricTag != null && Math.abs(Constants.reporterMassOfIsobaricTag[0]-ptm.getMassDifference()) < 0.01 ){
				if( ptm.getAbbAA() == 'S' || ptm.getAbbAA() == 'T' ) {
					if( "itraq4plex".compareTo(Constants.isobaricTag) == 0 ) ptm.setDiagnosticIon(163.1199);
					ptm.setPenalty(0.5);
					ptm.setModCount(0);
					AminoAcid.canBeEasilyModified(ptm.getAbbAA());
				}
			}
			
			//setting specific modifications
			if( ptm.getName().startsWith("Acetyl") ) {
				if( ptm.getAbbAA() =='K' ) ptm.setDiagnosticIon(126.0913);
				if( ptm.getPTMPosition() == PTMPosition.PROTEIN_N_TERM ) ptm.setPenalty(0.5);				
				
				if( "Acetyl".compareToIgnoreCase(Constants.enrichedModification) == 0 ) {
					ptm.setPenalty(0.5);					
					if( ptm.getAbbAA() =='K' ) {
						ptm.setModCount(0);
						AminoAcid.canBeEasilyModified('K');
					}
					if( ptm.getPTMPosition() == PTMPosition.PROTEIN_N_TERM ) {
						ptm.setModCount(0);
					}
				}
			}
			
			if( "Phospho".compareTo(ptm.getName()) == 0 ){
				if( ptm.getAbbAA() =='S' || ptm.getAbbAA() =='T' ) {
					ptm.setNeutralLoss(ptm.getMassDifference()+Constants.H2O);
					if( "Phospho".compareToIgnoreCase(Constants.enrichedModification) == 0 ) {
						ptm.setPenalty(0.5);
						ptm.setModCount(0);
						AminoAcid.canBeEasilyModified(ptm.getAbbAA());
					}
				}
			}
		}
	}
	
}


