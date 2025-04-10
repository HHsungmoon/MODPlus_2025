package modi;

import java.util.ArrayList;
import java.util.Collections;

import msutil.PGraph;

public class AnsHeap {
	static int Capacity = 100; // max heap capacity
	int Size; // current heap size
	int nth = 10; //print
	AnsPeptide[] Elements;
	
	public AnsHeap() {
		Size = 0;
		Elements = new AnsPeptide[Capacity+1];
		Elements[0] = new AnsPeptide(null, "", null, -1000);
	}
	
	public void add( AnsPeptide X ) {
		int i=0;
		if( Size == Capacity ){
			if( Elements[1].compareTo(X) > 0 ) deleteMin();
			else return;
		}
		
		for( i = ++Size ; Elements[i/2].compareTo(X) < 0 ; i /= 2 )
			Elements[i] = Elements[i/2];
		
		Elements[i] = X;
	}
	
	public void deleteMin() {
		if( Size == 0 ) return;	
		int i, Child;
		AnsPeptide LastElement = Elements[Size--];		
		
		for( i = 1; i*2 <= Size; i = Child )
		{
			Child = i * 2;
			if( Child != Size && Elements[Child+1].compareTo(Elements[Child]) > 0 )
				Child++;
			if( LastElement.compareTo(Elements[Child]) < 0 )
				Elements[i] = Elements[Child];
			else
				break;
		}
		Elements[ i ] = LastElement;
	}

	public ArrayList<AnsPeptide> getFinalList(PGraph graph) {
		ArrayList<AnsPeptide> tempList = new ArrayList<AnsPeptide>();
		for( int i=1; i<= Size; i++){
			Elements[i].evaluatePSM(graph);	
			tempList.add(Elements[i]);
		}			
		Collections.sort( tempList );
		
		int toIndex = ( Size > nth )? nth : Size;
		ArrayList<AnsPeptide> fList = new ArrayList<AnsPeptide>();
		fList.addAll( tempList.subList(0, toIndex) );

		return fList;		
	}
	
}