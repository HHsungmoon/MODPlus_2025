package processedDB;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedList;

import modi.Constants;
import modi.IonDirection;
import modi.Tag;

import msutil.MSMass;
import msutil.ProtCutter;

public class TagTrie extends ProtDatabase {	
	
	private static int 		AASIZE= 26;
	private static double 	isoShift= 2;
	private int[][][]		tagHit;
	private int[][][][] 	tagTrie;
	private int stem = 0;

	public TagTrie( String fileName ) throws Exception{		
		super( fileName );
		constructTagTrie();
	}
	
	public TagTrie( ProxDB pdb ) throws Exception {		
		super( pdb );
		constructTagTrie();
	}
	public void setStemNo(int s) { stem = s; }
	public int getStemNo() { return stem; }
	
	private void constructTagTrie() {	
		
		tagHit  = new int[AASIZE][AASIZE][AASIZE];
		for( int k=0; k<at_end_of_line; k++ ){
			if( sequence[k] == delimeter ){				
				int i = k + 3;
				for( ; sequence[i] != delimeter; i++){					
					tagHit[sequence[i-2]-'A'][sequence[i-1]-'A'][sequence[i]-'A']++;
				}
				k= i-1;			
			}
		}
		
		tagTrie = new int[AASIZE][AASIZE][AASIZE][];
		for (int i=0; i<AASIZE; i++){
			for (int j=0; j<AASIZE; j++){
				for (int k=0; k<AASIZE; k++){					
					tagTrie[i][j][k] = new int[tagHit[i][j][k]];
				}
			}
		}	
		
		int[][][] tIndex = new int[AASIZE][AASIZE][AASIZE];		
		for( int k=0; k<at_end_of_line; k++ ){
			if( sequence[k] == delimeter ){				
				int i = k + 3;
				for( ; sequence[i] != delimeter; i++){									
					tagTrie[sequence[i-2]-'A'][sequence[i-1]-'A'][sequence[i]-'A'][tIndex[sequence[i-2]-'A'][sequence[i-1]-'A'][sequence[i]-'A']++] = i-2;
				}
				k= i-1;			
			}
		}
	}
	
	public LinkedList<MODPeptide> getOneModPeptides(double nGap, String tag, double cGap, int tagType,
			double minModified, double maxModified, double gapTolerance) {	
		double tolerance = gapTolerance + isoShift;
		LinkedList<MODPeptide> matLocal = new LinkedList<MODPeptide>();
		HashSet<String> peptPool = new HashSet<String>();
		
		int	listSize  = tagHit[tag.charAt(0)-'A'][tag.charAt(1)-'A'][tag.charAt(2)-'A'];
		int[] tagList = tagTrie[tag.charAt(0)-'A'][tag.charAt(1)-'A'][tag.charAt(2)-'A'];	
					
		for( int thx=0; thx<listSize; thx++ ){	
			
			int tagIndex = tagList[thx];
			
			boolean isGapMatched = false;
			int left = 0, right = 0;
			
			int preLen= -1, sufLen= -1;
			double preMass= 0, sufMass= 0;
			
			int start = tagIndex;
			double delta = nGap-preMass;
			while( delta > minModified ){	
				if( delta < maxModified ) preLen++;	
				if( Math.abs(delta) <= tolerance ) {
					left = start;
					isGapMatched = true;
					start--;
					break;
				}
					
				start--;
				if( sequence[start]==delimeter || MSMass.getAAMass( sequence[start] ) == 0 ) break;
				preMass += MSMass.getAAMass( sequence[start] );
				delta = nGap-preMass;
			}
			
			int end = tagIndex+2;
			delta = cGap-sufMass;
			while( delta > minModified ){	
				if( delta < maxModified ) sufLen++;			
				if( Math.abs(delta) <= tolerance ) {
					right = end+1;
					isGapMatched = true;
					end++;
					break;
				}
				
				end++;
				if( sequence[end]==delimeter || MSMass.getAAMass( sequence[end] ) == 0 ) break;
				sufMass += MSMass.getAAMass( sequence[end] );
				delta = cGap-sufMass;
			}
			
			if( !isGapMatched ) continue;
			
			if( preLen > -1 && sufLen > -1 ){
				start++;
				if( !peptPool.add(getPeptide(start, end)) ) continue;
				
				MODPeptide mPept = new MODPeptide(start, end);
				if( left == 0 ) left = start + preLen;
				if( right == 0 ) right = end - sufLen;
				mPept.setConservedRegion( left, right );
				matLocal.add( mPept );
			}
		}
		return matLocal;
	}
	
	public LinkedList<TagPeptide> getMultiModPeptides(double nGap, String tag, double cGap, int tagType,
			double minModified, double maxModified, double gapTolerance) {	
		
		LinkedList<TagPeptide> matLocal = new LinkedList<TagPeptide>();
		HashSet<String> peptPool = new HashSet<String>();
		
		int	listSize  = tagHit[tag.charAt(0)-'A'][tag.charAt(1)-'A'][tag.charAt(2)-'A'];
		int[] tagList = tagTrie[tag.charAt(0)-'A'][tag.charAt(1)-'A'][tag.charAt(2)-'A'];	

		for( int thx=0; thx<listSize; thx++ ){	
			
			int tagIndex = tagList[thx];
			
			int preLen= -1, sufLen= -1;
			double preMass= 0, sufMass= 0;
			
			int start = tagIndex;
			double delta = nGap-preMass;
			while( delta > minModified ){	
				if( delta < maxModified ) preLen++;				
				start--;
				if( sequence[start] == delimeter || MSMass.getAAMass( sequence[start] ) == 0 ) break;
				preMass += MSMass.getAAMass( sequence[start] );
				delta = nGap-preMass;
			}
			
			int end = tagIndex+2;
			delta = cGap-sufMass;
			while( delta > minModified ){	
				if( delta < maxModified ) sufLen++;				
				end++;
				if( sequence[end] == delimeter || MSMass.getAAMass( sequence[end] ) == 0 ) break;
				sufMass += MSMass.getAAMass( sequence[end] );
				delta = cGap-sufMass;
			}
			
			if( preLen > -1 && sufLen > -1 ){
				start++;
				if( !peptPool.add( getPeptide(start, end)+(tagIndex-start) ) ) continue;
				SequenceTag seqTag = new SequenceTag(nGap-MSMass.getPepMass(getPeptide(start, tagIndex)), 
						tagIndex-start, tagIndex-start+2, tagType);	
				TagPeptide xPept = new TagPeptide(start, end, seqTag);
				xPept.setConservedRegion( start+preLen, end-sufLen );
				matLocal.add( xPept );
			}
		}
		return matLocal;
	}

	public RetrivedPeptideMap getRetrivedPeptides(double orbMass, ProtCutter cutter, int NTT, double nGap, Tag queryTag, double cGap, 
			IonDirection ir, double minModified, double maxModified, double gapTolerance) {	
		//for modeye identification in second search				
		double tolerance = gapTolerance + isoShift;	
		RetrivedPeptideMap rpMap = new RetrivedPeptideMap();

		String tag= queryTag.sequence().toString();
		double tagMass= MSMass.getPepMass(tag);
		int	listSize  = tagHit[tag.charAt(0)-'A'][tag.charAt(1)-'A'][tag.charAt(2)-'A'];
		int[] tagList = tagTrie[tag.charAt(0)-'A'][tag.charAt(1)-'A'][tag.charAt(2)-'A'];
		for( int thx=0; thx<listSize; thx++ ){
			
			int tagIndex = tagList[thx];
			
			int start, end;
			int preLen=0, sufLen=0;
			double preMass= 0, sufMass= 0;
		
			if( Math.abs(nGap) <= tolerance ){ 
				start = tagIndex-1;
				preLen++;
			}
			else{
				for(start=tagIndex-1; sequence[start]!=delimeter ; start--){		
					if( MSMass.getAAMass( sequence[start] ) == 0 ) break;
					preMass += MSMass.getAAMass( sequence[start] );
					double offset = nGap-preMass;
					if( offset < maxModified &&  offset > minModified )
						preLen++;
					else if( offset < minModified )
						break;					
				}
			}
			
			if( Math.abs(cGap) <= tolerance ){ 
				end=tagIndex+3;
				sufLen++;
			}
			else{
				for(end=tagIndex+3; sequence[end]!=delimeter ; end++){
					if( MSMass.getAAMass( sequence[end] ) == 0 ) break;
					sufMass += MSMass.getAAMass( sequence[end] );
					double offset = cGap-sufMass;
					if( offset < maxModified && offset > minModified )
						sufLen++;
					else if( offset < minModified )
						break;
				}
			}
			
			double ntermMass, ctermMass;
			if( start+1+preLen < tagList[thx] )
				ntermMass = MSMass.getPepMass(getPeptide(start+1+preLen, tagIndex));
			else
				ntermMass = - MSMass.getPepMass(getPeptide(tagIndex, start+1+preLen));
			if( tagList[thx]+3 < end-sufLen )
				ctermMass = MSMass.getPepMass(getPeptide(tagIndex+3, end-sufLen));
			else
				ctermMass = - MSMass.getPepMass(getPeptide(end-sufLen, tagIndex+3));
		
			double sharedMass= ntermMass + tagMass + ctermMass;
			if( sharedMass < 0 ) continue;
			
			preMass= 0; 
			sufMass= 0;
			int nterm=0, cterm=0;
			
			for(int i=preLen; i>0 ; i--){
				int n_leftAAPos = start+i;
				preMass += MSMass.getAAMass( sequence[n_leftAAPos] );
				nterm= 0;
				if( cutter.isCleavage(sequence[n_leftAAPos-1], sequence[n_leftAAPos]) || isPotentialProteinNterm(n_leftAAPos) )
					nterm++;
				sufMass=0;
				for(int j=sufLen; j>0; j--){
					int c_rightAAPos = end-j;
					sufMass += MSMass.getAAMass( sequence[c_rightAAPos] );
					cterm= 0;
					if( cutter.isCleavage(sequence[c_rightAAPos], sequence[c_rightAAPos+1]) )
						cterm++;
					
					if( nterm + cterm < NTT ) continue;
					if( NTT != 0 ) {
						int missCleavage = 0;
						for(int k=n_leftAAPos; k<c_rightAAPos; k++){
							if( cutter.isCleavage(sequence[k], sequence[k+1]) ) missCleavage++;
						}
						if( missCleavage > Constants.missCleavages ) continue;
					}
					
					double deltaM= orbMass - (preMass+sharedMass+sufMass+Constants.NTERM_FIX_MOD+Constants.CTERM_FIX_MOD);	
					if( ( deltaM < maxModified && deltaM > minModified ) || Math.abs(deltaM) <= tolerance  ){						
						int peptPos = 1;					
						if( sequence[n_leftAAPos-1] == delimeter || isPotentialProteinNterm(n_leftAAPos) ) peptPos = 0;//protein nterm
						else if( sequence[c_rightAAPos+1] == delimeter ) peptPos = 2;//protein cterm
						RetrivedPeptide rp= new RetrivedPeptide( getPeptide(n_leftAAPos, c_rightAAPos+1), peptPos, n_leftAAPos, c_rightAAPos, nterm+cterm );		
						rp.addTag( queryTag, nGap-(preMass+ntermMass), cGap-(sufMass+ctermMass), tagIndex-n_leftAAPos, tagIndex-n_leftAAPos+2, ir );
						rpMap.add(rp);	
					}
				}
			}	
		}									
		return rpMap;
	}
	
	public ArrayList<PeptideMatchToProtein> getMatchProteins(String peptide, int prior){
		
		ArrayList<PeptideMatchToProtein> protMatch = new ArrayList<PeptideMatchToProtein>();
		
		String tag= peptide.substring(0, 3);
		int posShift = peptide.length()-3;
		String end= peptide.substring(posShift, peptide.length());

		int	ASize  = tagHit[tag.charAt(0)-'A'][tag.charAt(1)-'A'][tag.charAt(2)-'A'];
		int[] AList = tagTrie[tag.charAt(0)-'A'][tag.charAt(1)-'A'][tag.charAt(2)-'A'];
		
		int	BSize  = tagHit[end.charAt(0)-'A'][end.charAt(1)-'A'][end.charAt(2)-'A'];
		int[] BList = tagTrie[end.charAt(0)-'A'][end.charAt(1)-'A'][end.charAt(2)-'A'];
		
		ArrayList<Integer> matPos = new ArrayList<Integer>();
		int a = 0, b=0;
		while( a < ASize && b < BSize ){
			if( AList[a] > BList[b]-posShift ) b++;
			else if( AList[a] < BList[b]-posShift ) a++;
			else { // matched
				matPos.add(AList[a]);
				a++;
				b++;
			}
		}		
		
		for( int thx=0; thx<matPos.size(); thx++ ){
			int start = matPos.get(thx);
			
			int comp_pos = start+3;			
			boolean matched = true;
			for(int i=3; i<posShift; i++){
				if( sequence[comp_pos] != peptide.charAt(i) ){
					matched = false;
					break;
				}
				comp_pos++;
			}
			
			if( matched ){				
				int endPos = start+peptide.length()-1;
				char prevAA = (char)sequence[start-1];
				char nextAA = (char)sequence[endPos+1];

				int ntt = 0;
				if( Constants.protease.isCleavage(sequence[endPos], nextAA) ) ntt++; // checking Cterm			
				if( Constants.protease.isCleavage(prevAA, peptide.charAt(0)) || isPotentialProteinNterm(start) ) ntt++;  // checking Nterm
				if( ntt < Constants.numberOfEnzymaticTermini ) continue;
				
				int protein = searchProtein(start);
				String pname =  proteins[protein]; 
				int realPos = start;
				if( protein != 0 ) realPos = start - protPos[protein-1];
				protMatch.add( new PeptideMatchToProtein(pname, (protein+prior+1), prevAA, realPos, realPos+peptide.length()-1, nextAA, ntt)  );			
			}
		}
		
		return protMatch;
	}
	
}





































