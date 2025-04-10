package msutil;

import modi.Constants;
import modi.PTM;

public class Scoring {
	
	//MODa Scoring
	public static int getModARankScore( String peptide, double[] ptms, PGraph graph ){ //for modA preliminary ranking	 
		IonGraph iGraph;
		if( Constants.INSTRUMENT_TYPE == Constants.msms_type.QTOF ) iGraph = new TOFGraph(peptide, ptms, graph);
		else iGraph= new TRAPGraph(peptide, ptms, graph);
		
		iGraph.setScore(graph);
		return iGraph.getRankScore();
	}
	
	public static IonGraph PeptideSpectrumMatch( String peptide, double[] ptms, PGraph graph ){//for moda final scoring
		IonGraph iGraph;
		if( Constants.INSTRUMENT_TYPE == Constants.msms_type.QTOF ) iGraph = new TOFGraph(peptide, ptms, graph);
		else iGraph= new TRAPGraph(peptide, ptms, graph);
		
		iGraph.evaluateMatchQuality(graph);
		return iGraph;
	}
	
    //MODeye Scoring
	public static int getModEyeRankScore( String peptide, double[] ptms, PGraph graph ){	//for modeye preliminary ranking
		IonGraph iGraph;
		if( Constants.INSTRUMENT_TYPE == Constants.msms_type.QTOF ) iGraph = new TOFGraph(peptide, ptms, graph);
		else iGraph= new TRAPGraph(peptide, ptms, graph);
		
		if( !Constants.isWithinTolerance(iGraph.getCalculatedMW(), graph.getObservedMW(), Constants.precursorTolerance) ) 
			return -1;
		
		iGraph.setScore(graph);
		return iGraph.getRankScore();
	}
	
	public static IonGraph PeptideSpectrumMatch( String peptide, double[] ptms, PTM[] ptmList, PGraph graph ){//for modeye final scoring
		IonGraph iGraph;
		if( Constants.INSTRUMENT_TYPE == Constants.msms_type.QTOF ) iGraph = new TOFGraph(peptide, ptms, ptmList, graph);
		else iGraph= new TRAPGraph(peptide, ptms, ptmList, graph);
		
		iGraph.evaluateMatchQuality(graph);
		return iGraph;
	}	
	
	public static double getOddProbability( double tod ){
		return Math.exp(tod)/(1+Math.exp(tod));
	}
	
}






















