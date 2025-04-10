package msutil;

import modi.Constants;
import modi.ThreadPoolManager;


public class PRM {
	double MW, prmMW;
	
	static int accuracy = ( Constants.MSMSResolution==1 )? 100 : 10;
	static int tolerance = (Constants.fragmentTolerance < 1.0/accuracy)? 1 : (int)(Constants.fragmentTolerance*accuracy);

	double[] bTable, yTable;
	
	public PRM(PGraph graph){
		
		MW = graph.getCorrectedMW();
		prmMW = MW - Constants.H2O;

		int matSize= (int)Math.round(MW*accuracy);
		
		graph.setPRMScores(1);
		
		bTable= new double[matSize];
		yTable= new double[matSize];
		
		for( int pn=0; pn<graph.size(); pn++ ){
			if( graph.get(pn).getMass() < 56 || graph.get(pn).getMass() > (MW-74) ) continue;
									
			int index = (int)Math.round( (graph.get(pn).getMass()-Constants.Proton)*accuracy );
			for( int i=index-tolerance; i<=index+tolerance; i++ ){
				if( bTable[i] < graph.get(pn).getBPRMScore() )
					bTable[i] = graph.get(pn).getBPRMScore();			
			}
			
			index= (int)Math.round( (MW-graph.get(pn).getMass()+Constants.Proton)*accuracy );
			for( int i=index-tolerance; i<=index+tolerance; i++ ){
				if( yTable[i] < graph.get(pn).getYPRMScore() )
					yTable[i] = graph.get(pn).getYPRMScore();				
			}
		}
	}
	
	public double getPeptMass(){ return prmMW; }
	
	public double getScore(double mass){
		return bTable[(int)(mass*accuracy)] + yTable[(int)(mass*accuracy)];
	}
	
	public double getScore(double mass, double delta){
		return bTable[(int)(mass*accuracy)] + yTable[(int)((mass+delta)*accuracy)];
	}

	public double massCorrection( boolean dynamicCorrection )
	{
		int slack = 56*accuracy;
		int size = bTable.length - slack;
	
		int startpoint;
		int slot = ThreadPoolManager.getSlotIndex();
		int endpoint = (int)(Constants.precursorTolerance[slot] * accuracy);
		if( Constants.maxNoOfC13[slot] == 0 ){
			startpoint	= 	-endpoint;			
		}//*/
		else{
			startpoint = (int)((Constants.minNoOfC13 - Constants.precursorAccuracy[slot]) * accuracy);
		}
		
		if( slack < Math.abs(startpoint) || slack < Math.abs(endpoint) ) return 0;
		
		double maxConv = 0, secondConv = 0;
		int maxShift = 0, secondShift=0;
		
		for( int pos=startpoint; pos<=endpoint; pos++){
			double conv = 0;			
			for( int i=slack; i<size ; i++){
				conv += bTable[i]*yTable[i+pos];
			}
			
			if( conv > maxConv ){
				secondConv = maxConv;
				secondShift = maxShift;
				maxConv = conv;
				maxShift = pos;
			}
			else if( conv > secondConv ){
				secondConv = conv;
				secondShift = pos;
			}
		}
		if( maxConv == 0 ) maxShift = secondShift = 0;
		else if( dynamicCorrection && maxShift > endpoint-tolerance ){
			int newendpoint = endpoint*2;			
			for( int pos=endpoint; pos<=newendpoint; pos++){
				double conv = 0;				
				for( int i=slack; i<size ; i++){
					conv += bTable[i]*yTable[i+pos];
				}
				
				if( conv > maxConv ){
					secondConv = maxConv;
					secondShift = maxShift;
					maxConv = conv;
					maxShift = pos;
				}
				else if( conv > secondConv ){
					secondConv = conv;
					secondShift = pos;
				}
			}
			endpoint = newendpoint;
		}		
		
		return ((double)(maxShift+secondShift)/2)/accuracy;
	}
	
	
}











