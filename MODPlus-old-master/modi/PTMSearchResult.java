package modi;

import java.util.ArrayList;

public class PTMSearchResult 
{
	static int currentID = 0;
	ArrayList<PTMRun> ptmRun;
	int id;
	boolean interpreted;
	
	public PTMSearchResult(ArrayList<PTMRun> ptmRun, boolean interpreted) {
		this.id = currentID++;
		this.ptmRun = ptmRun;
		this.interpreted = interpreted;
	}
	public	ArrayList<PTMRun>	getPTMRun()		{ return ptmRun; }
	public	int					getID()			{ return id; }
	public	boolean				isInterpreted()	{ return interpreted; }

}
