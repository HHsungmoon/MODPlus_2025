package modi;

import java.util.ArrayList;

public class PTMRun extends ArrayList<PTMOccurrence> {

	double score = 0.;
	double error = 0 ;
	
	public void setError(double err){ error = err; }
	public double getError(){ return error; }
	
	public String toString(int basePos) {
		StringBuffer output = new StringBuffer();
		boolean start = true;

		for (PTMOccurrence ptmOcc : this)
		{
			if(start)
				start = false;
			else
				output.append("\\");
				
			output.append( ptmOcc.toString(basePos) );			
		}
		return output.toString();
	}
	
	public String toString(Peptide pept, int basePos){
		StringBuffer output = new StringBuffer();
		for ( PTMOccurrence ptmOcc : this ){
			int pos = basePos + ptmOcc.position;
			output.append( ptmOcc.ptm.getName() );
			output.append( "("+pept.get(pos).getResidue()+(pos+1)+")" );
			output.append( " " );
		}
		
		return output.toString();		
	}
	public double getPenalty(){
		double pen = 0;
		for (PTMOccurrence ptmOcc : this){			
			pen += ptmOcc.getPTM().getPenalty();
		}
		return pen;
	}
	public int setPTMMass(double[] ptms, int basePos){
		for (PTMOccurrence ptmOcc : this){			
			ptms[basePos + ptmOcc.position] = ptmOcc.ptm.getMassDifference();
		}
		return this.size();
	}
	public int setPTMs(PTM[] ptms, int basePos){
		for (PTMOccurrence ptmOcc : this){			
			ptms[basePos + ptmOcc.position] = ptmOcc.ptm;
		}	
		return this.size();
	}
}
