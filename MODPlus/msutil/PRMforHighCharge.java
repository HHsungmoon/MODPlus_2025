package msutil;

import msutil.PGraph;
import modi.Constants;

public class PRMforHighCharge extends PRM {
	
	static int extolerance= tolerance;
	
	public PRMforHighCharge(PGraph graph){
		super(graph);			
		graph.setPRMScores(2);
		
		double[] b2Table= new double[bTable.length];
		double[] y2Table= new double[yTable.length];
		
		for( int pn=0; pn<graph.size(); pn++ ){
			
			double singlyChargedM = graph.get(pn).getMass()*2-Constants.Proton;
			
			if( singlyChargedM < 56 || singlyChargedM > (MW-74) ) continue;						
			
			int index = (int)Math.round( (singlyChargedM-Constants.Proton)*accuracy );

			for( int i=index-extolerance; i<=index+extolerance; i++ ){
				if( b2Table[i] < graph.get(pn).getBPRMScore() );
					b2Table[i] = graph.get(pn).getBPRMScore();
			}
			
			index= (int)Math.round( (MW-singlyChargedM+Constants.Proton)*accuracy );
			for( int i=index-extolerance; i<=index+extolerance; i++ ){
				if( y2Table[i] < graph.get(pn).getYPRMScore() );
					y2Table[i] = graph.get(pn).getYPRMScore();
			}
		}
		
		for(int i=0; i<bTable.length; i++){
			bTable[i] += b2Table[i];
			yTable[i] += y2Table[i];
		}		
	}
}

