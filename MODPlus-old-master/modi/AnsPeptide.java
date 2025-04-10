package modi;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.StringTokenizer;

import org.jdom.Attribute;
import org.jdom.Element;

import msutil.IonGraph;
import msutil.PGraph;
import msutil.Scoring;
import processedDB.PeptideMatchToProtein;
import processedDB.ProtDatabase;

public class AnsPeptide implements Comparable<AnsPeptide> {

	static double  pointRounding = 1000;
	Peptide 	peptide;
	String 		ptmComb="";
	double[] 	ptms;
	PTM[] 		ptmList;
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
	
	public Peptide getPeptide(){ return peptide; }
	public String getPeptideSequence(){ return peptide.toString(); }
	public String getPTM(){ return ptmComb; }
	public double getScore(){ return score; }

	public void evaluatePSM(PGraph pg){//using ptm information
		IonGraph iG = Scoring.PeptideSpectrumMatch(peptide.toString(), ptms, ptmList, pg );
		score = iG.getRankScore();
		prob 	= iG.getProb();
		roundprob 	= Constants.round(iG.getProb()*pointRounding)/pointRounding;
		
		mzDev 	= iG.getMassError();		
		pepMass = iG.getCalculatedMW();		
		score -= (2-peptide.getNTT())*3;
	}

	public Element get_peptide_hit_element( int rank, ProtDatabase proDB, double obmw, PTMDB ptmDB ){ //for modeye output
		Element peptide_hit = new Element( "peptide_hit" );
		peptide_hit.setAttribute(new Attribute( "rank", String.valueOf(rank) ));
		peptide_hit.setAttribute(new Attribute( "sequence", peptide.toString() ));
		peptide_hit.setAttribute(new Attribute( "prev_next_res", 
				String.valueOf( proDB.getTeminalResiduesOfPeptide( peptide ) ) ));
		peptide_hit.setAttribute(new Attribute( "calculated_MW", String.format("%.4f", pepMass) ));
		
		Element delta_mass = new Element( "delta_mass" );
		delta_mass.setAttribute(new Attribute( "value", String.format("%.4f", obmw-pepMass) ));
		peptide_hit.addContent( delta_mass );
		
		Element score = new Element( "score" );
		score.setAttribute(new Attribute( "value", String.format("%.4f", this.prob) ));
		peptide_hit.addContent( score );
		
		StringTokenizer token = new StringTokenizer(ptmComb);
		if( token.countTokens() > 0 ){
			Element pep_modifications = new Element( "pep_modifications" );
			while( token.hasMoreTokens() ){
				Element modified = new Element( "modified" );
				String tp = token.nextToken();
				int delimiter = tp.lastIndexOf('(');
				modified.setAttribute(new Attribute( "site", tp.substring(delimiter+2, tp.length()-1) ));
				modified.setAttribute(new Attribute( "mod", 
						String.valueOf( ptmDB.getPTM(tp.substring(0, delimiter), tp.charAt(delimiter+1)).getID()+1) ));
				pep_modifications.addContent( modified );
			}
			peptide_hit.addContent( pep_modifications );
		}
			
		SrcProteinInfo prot = peptide.getSrcProteinInfo().get(0);
		Element protein = new Element( "protein" );
		protein.setAttribute(new Attribute( "name", proDB.getProteinIdentity(prot.getStartPos()) ));
		protein.setAttribute(new Attribute( "index", String.valueOf(proDB.getProteinIndex(prot.getStartPos())+1) ));
		protein.setAttribute(new Attribute( "pep_start", String.valueOf(proDB.getRealPositionInProtein( prot.getStartPos() )) ));
		peptide_hit.addContent( protein );
		
		if( peptide.getSrcProteinInfo().size() > 1 ){
			Element alternative_proteins = new Element( "alternative_proteins" );
			alternative_proteins.setAttribute(new Attribute( "num", String.valueOf(peptide.getSrcProteinInfo().size()-1) ));
			for( int pr=1; pr<peptide.getSrcProteinInfo().size(); pr++){
				Element entry = new Element( "entry" );
				entry.setAttribute(new Attribute( "name", proDB.getProteinIdentity(peptide.getSrcProteinInfo().get(pr).getStartPos()) ));
				entry.setAttribute(new Attribute( "index", String.valueOf(proDB.getProteinIndex(peptide.getSrcProteinInfo().get(pr).getStartPos())+1) ));
				entry.setAttribute(new Attribute( "pep_start", String.valueOf(proDB.getRealPositionInProtein( peptide.getSrcProteinInfo().get(pr).getStartPos() )) ));
				alternative_proteins.addContent( entry );
			}		
			peptide_hit.addContent( alternative_proteins );
		}
		return peptide_hit;
	}
	
	public String toMODPlus( double obMW, ArrayList<PeptideMatchToProtein> protMatch ){
		//For MOD Plus, OrgDB & DynamicDB
		StringBuffer x= new StringBuffer();	
		x.append( String.format("%.4f\t", pepMass) );
		x.append( String.format("%.4f\t", obMW-pepMass) );		

		x.append( String.format("%d\t", score) );
		x.append( String.format("%.4f\t", prob) );

		String peptStr = peptide.toString();
		String wrapAA = protMatch.get(0).getWrappingAA();
		x.append( wrapAA.charAt(0)+".");
		for( int i=0; i<ptms.length; i++){
			x.append(peptStr.charAt(i));
			if( ptms[i] != 0 ){
				if( ptms[i] > 0 ) x.append("+");
				x.append( String.format("%.3f", ptms[i]) );
			}
		}
		x.append( "."+wrapAA.charAt(1)+"\t" );		
		x.append( protMatch.get(0).toString() );
		
		for( int i=1; i<protMatch.size(); i++ ){
			x.append( ";" + protMatch.get(i).toString() );
		}
		x.append( "\t"+ptmComb );
		
		return x.toString();
	}
	
	public Element get_peptide_hit_element( int rank, double obmw, PTMDB ptmDB, ArrayList<PeptideMatchToProtein> protMatch ){
		//For MOD Plus, OrgDB & DynamicDB
		Element peptide_hit = new Element( "peptide_hit" );
		peptide_hit.setAttribute(new Attribute( "rank", String.valueOf(rank) ));
		peptide_hit.setAttribute(new Attribute( "sequence", peptide.toString() ));
		peptide_hit.setAttribute(new Attribute( "calculated_MW", String.format("%.4f", pepMass) ));
		
		Element delta_mass = new Element( "delta_mass" );
		delta_mass.setAttribute(new Attribute( "value", String.format("%.4f", obmw-pepMass) ));
		peptide_hit.addContent( delta_mass );
		
		Element score = new Element( "score" );
		score.setAttribute(new Attribute( "value", String.format("%.4f", this.prob) ));
		peptide_hit.addContent( score );

		StringTokenizer token = new StringTokenizer( ptmComb );
		if( token.countTokens() > 0 ){
			Element pep_modifications = new Element( "pep_modifications" );
			while( token.hasMoreTokens() ){
				Element modified = new Element( "modified" );
				String tp = token.nextToken();
				int delimiter = tp.lastIndexOf('(');
				modified.setAttribute(new Attribute( "site", tp.substring(delimiter+2, tp.length()-1) ));
				modified.setAttribute(new Attribute( "mod", 
						String.valueOf( ptmDB.getPTM(tp.substring(0, delimiter), tp.charAt(delimiter+1)).getID()+1) ));
				pep_modifications.addContent( modified );
			}
			peptide_hit.addContent( pep_modifications );
		}
		
		PeptideMatchToProtein pm = protMatch.get(0);
		Element protein = new Element( "protein" );
		protein.setAttribute(new Attribute( "name", pm.getProtName() ));
		protein.setAttribute(new Attribute( "index", String.valueOf(pm.getProtIndex()) ));
		protein.setAttribute(new Attribute( "pep_start", String.valueOf(pm.getStartSite()) ));
		peptide_hit.addContent( protein );	
		peptide_hit.setAttribute(new Attribute( "prev_next_res", pm.getWrappingAA() ) );
		
		if( protMatch.size() > 1 ){
			Element alternative_proteins = new Element( "alternative_proteins" );
			alternative_proteins.setAttribute(new Attribute( "num", String.valueOf(protMatch.size()-1) ));
			for( int i=1; i<protMatch.size(); i++ ){
				Element entry = new Element( "entry" );
				entry.setAttribute(new Attribute( "name", protMatch.get(i).getProtName() ));
				entry.setAttribute(new Attribute( "index", String.valueOf(protMatch.get(i).getProtIndex()) ));
				entry.setAttribute(new Attribute( "pep_start", String.valueOf(protMatch.get(i).getStartSite()) ));
				alternative_proteins.addContent( entry );
			}		
			peptide_hit.addContent( alternative_proteins );
		}
		return peptide_hit;
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
		if(!(o instanceof AnsPeptide))
			   return false;		
		AnsPeptide p= (AnsPeptide)o;
		if( this.peptide.equals(p.peptide) && this.ptmComb.equals(p.ptmComb) )
			return true;
		return false;
	}
	
	public double getNumOfComb(){
		double numComb = 0;
		HashSet<String> t = new HashSet<String>();		
		
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

