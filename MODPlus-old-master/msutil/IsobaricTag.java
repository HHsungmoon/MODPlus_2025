package msutil;

public class IsobaricTag {

	public static double[] getReporterMasses(String name){
		
		     if( "iTRAQ4plex".compareToIgnoreCase(name) == 0 ) return massTable[0];
		else if( "iTRAQ8plex".compareToIgnoreCase(name) == 0 ) return massTable[1];
		     
		else if( "TMT2plex".compareToIgnoreCase(name) == 0 ) return massTable[2];	
		else if( "TMT6plex".compareToIgnoreCase(name) == 0 ) return massTable[3];	
		     
		return null;
	}
	
	private static final double[][] massTable =  {
		//index 0: reagent mass, reporters are sorted
		{144.102063, 114.11123, 115.10826, 116.11162, 117.11497}, //iTRAQ 4plex
		{304.205363, 113.107873, 114.111228, 115.108263, 116.111618, 117.114973, 118.112008, 119.115363, 121.122072}, //iTRAQ 8plex
		{225.155833, 126.1283, 127.1316}, // TMT 2plex
		{229.162932, 126.1283, 127.1316, 128.1350, 129.1383, 130.1417, 131.1387 }, // TMT 6plex		
	};

}
