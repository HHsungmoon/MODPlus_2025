package processedDB;

import java.util.ArrayList;

public class ProxDB extends ArrayList<Prox> {
	
	int sizeOfEntries;
	int sizeOfResidues;	
	String name, orgName;
	
	public int getSizeOfResidues() { return sizeOfResidues; }
	public String getName() { return orgName; }
	
	public void setSizeOfEntries(int x) { sizeOfEntries = x; }
	public void setSizeOfResidues(int x) { sizeOfResidues = x; }
}























