package moda;

import java.io.PrintWriter;
import java.util.Collections;
import java.util.LinkedList;

import processedDB.StemTagTrie;

import msutil.PGraph;
import msutil.Scoring;
import scaniter.MSMScan;

public class DPHeap extends LinkedList<DPPeptide> {
	
	static int Capacity = 100; // max heap capacity //MODplus
		
	public DPHeap(){
		for(int i=0; i<Capacity; i++ )
			this.add( new DPPeptide() );
	}
	
	public void setStemNo(int stem){
		for(int i=0; i<this.size(); i++) this.get(i).stem = stem;
	}
	
	public boolean isConfident() {
		for( int i=0; i<this.size()&&i<10; i++ ){
			if( this.get(i).isConfident() ) return true;
		}
		return false;
	}
	
	public boolean insert(DPPeptide dp){	
		
		if( dp.score < 1 || dp.compareTo( this.getLast() ) == 1 ) return false;
		
		int i=0;
		for( DPPeptide x : this ) {
			if( dp.isSame(x) ) break;			
			if( x.compareTo(dp) == 1 ){
				this.add(i, dp);
				this.removeLast();
				break;
			}
			i++;
		}
		return true;
	}
	
	public void insertAll(DPHeap heap){		
		for( DPPeptide x : heap ){
			this.insert(x);
		}
	}

	public int evaluate( PGraph graph ) { //for moda
		int i = 0;
		int maxScore = 0;
		for( i=0; i<this.size(); i++ ){
			if( this.get(i).score < 1 ) {
				this.remove(i);
				i--;
				continue;
			}
			this.get(i).evaluatePSM(graph);
			if( maxScore < this.get(i).score ) maxScore = this.get(i).score;
		}
		
		for( i=0; i<this.size(); i++ ){
			if( this.get(i).score == maxScore ) {
				this.get(i).score -= 1;
			}
		}		
		Collections.sort( this, new DPPeptideRefinedComparator() );		
		if( this.size() == 0 || this.get(0).score < 1 ) return 0;	
		return i;
	}

	public int reArr( PGraph graph ) { //for modplus
		int i = 0;
		int maxScore = 0;
		for( i=0; i<this.size(); i++ ) {
			if( this.get(i).score < 1 ) {
				this.remove(i);
				i--;
				continue;
			}
			this.get(i).score = Scoring.getModARankScore(this.get(i).peptide, this.get(i).ptms, graph);
			if( maxScore < this.get(i).score ) maxScore = this.get(i).score;
		}
		
		Collections.sort( this, new DPPeptideRefinedComparator() );		
		int cut = maxScore/2;
		for( i=0; i<this.size(); i++ ){
			if( this.get(i).score < cut ) this.removeRange(i, this.size());
		}
		
		if( this.size() == 0 || this.get(0).score < 1 ) return 0;	
		return i;
	}
	
}

