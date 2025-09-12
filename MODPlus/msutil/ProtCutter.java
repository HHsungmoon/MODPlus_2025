package msutil;

import modi.Constants;

public class ProtCutter {
	
	private static int 	N_TERM = 1, C_TERM = 2, BOTH_TERM = 3;
	private static char TERMINAL_RES = '-';
	
	private String 		name="";
	private String		cleave="";
	private String		term="";
	private boolean 	specified;
	
	private String		nTermCleave="", cTermCleave="";
	private int[]		cleaveMap = new int[26]; 
	
	public ProtCutter( String n, String c, String t ){
		name 	= n;
		cleave 	= c;
		term 	= t;
		specified = true;

		for (int i=0; i<term.length(); i++){
			if( term.charAt(i) == 'C' ) cTermCleave = c;
			if( term.charAt(i) == 'N' ) nTermCleave = c;
		}
		setCleaveMap( nTermCleave, cTermCleave, true );
	}
	
	public ProtCutter( String na, String nc, String cc, boolean isEnzyme ){
		name = na;
		nTermCleave 	= nc;
		cTermCleave 	= cc;		
		specified		= isEnzyme;
		setCleaveMap( nTermCleave, cTermCleave, isEnzyme );
		
	}
	
	public void setCleaveMap( String nc, String cc, boolean isEnzyme ){

		if( !isEnzyme ){
			for (int i=0; i<cleaveMap.length; i++){
				cleaveMap[i] = C_TERM;
			}
			Constants.numberOfEnzymaticTermini = 0;
		}
		else{
			for (int i=0; i<nTermCleave.length(); i++){
				if( 'A' <= nTermCleave.charAt(i) && nTermCleave.charAt(i) <= 'Z' );
					cleaveMap[nTermCleave.charAt(i)-'A'] = N_TERM;
			}
			for (int i=0; i<cTermCleave.length(); i++){
				if( 'A' <= cTermCleave.charAt(i) && cTermCleave.charAt(i) <= 'Z' );{
					if( cleaveMap[cTermCleave.charAt(i)-'A'] == N_TERM ) cleaveMap[cTermCleave.charAt(i)-'A'] = BOTH_TERM;
					else cleaveMap[cTermCleave.charAt(i)-'A'] = C_TERM;
				}
			}
		}
	}
	public String getName()	{ return name; }
	public boolean isSpecified() { return specified; }


	public boolean isCleavage( int prev_aa, int next_aa ){
		if( prev_aa == TERMINAL_RES || next_aa == TERMINAL_RES ) return true;

		int tmp = cleaveMap[prev_aa-'A'];
		if( tmp == C_TERM || tmp == BOTH_TERM )
			return true;

		tmp = cleaveMap[next_aa-'A'];
		if( tmp == N_TERM || tmp == BOTH_TERM )
			return true;

		return false;
	}

	
	public static ProtCutter getCutter( String name ){
		for(ProtCutter e : enzymeTable)
			if( e.getName().equalsIgnoreCase(name) )
				return e;
		return null;
	}
	
	private static final ProtCutter[] enzymeTable = 
	{
		new ProtCutter("Trypsin",		"KR",		"C" ),
		new ProtCutter("Arg-C",			"R",		"C" ),
		new ProtCutter("Glu-C",			"E",		"C" ),
		new ProtCutter("Asp-N",			"BD",		"N" ),
		new ProtCutter("Asp-N_ambic",	"DE",		"N" ),
		new ProtCutter("Chymotrypsin",	"FYWL",		"C" ),
		new ProtCutter("CNBr",			"M",		"C" ),
		new ProtCutter("Formic_acid",	"D",		"C" ),
		new ProtCutter("Lys-C",			"K",		"C" ),
		new ProtCutter("Lys-C/P",		"K",		"C" ),
		new ProtCutter("PepsinA",		"FL",		"C" ),
		new ProtCutter("Tryp-CNBr",		"KRM",		"C" ),
		new ProtCutter("TrypChymo",		"FYWLKR",	"C" ),
		new ProtCutter("Trypsin/P",		"KR",		"C" ),
		new ProtCutter("V8-DE",			"BDEZ",		"C" ),
		new ProtCutter("V8-E",			"EZ",		"C" ),
		new ProtCutter("CNBr+Trypsin",	"KRM",		"C" ),
	};//*/

}
