package modi;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.TreeMap;

import moda.ThreadPoolManager;
import msutil.MSMass;

public class PTMDB extends ArrayList<PTM> {
	private ArrayList<PTM>[][]		PTMTable;
	private TreeMap<Integer, String> classifications = new TreeMap<Integer, String>();
	private double[][] modRange = new double[2][26];
	private double[][] ntermModRange = new double[2][26];
	private double[][] ctermModRange = new double[2][26];

	public PTMDB deepCopy() {
		PTMDB copy = new PTMDB();

		// PTM 리스트 복제
		for (PTM ptm : this) {
			copy.add(ptm.clone()); // PTM 클래스에 clone() 구현 필요
		}

		// PTMTable 복제
		if (this.PTMTable != null) {
			copy.PTMTable = new ArrayList[this.PTMTable.length][];
			for (int i = 0; i < this.PTMTable.length; i++) {
				copy.PTMTable[i] = new ArrayList[this.PTMTable[i].length];
				for (int j = 0; j < this.PTMTable[i].length; j++) {
					copy.PTMTable[i][j] = new ArrayList<>(this.PTMTable[i][j]);
				}
			}
		}

		// classifications 복제
		copy.classifications = new TreeMap<>(this.classifications);

		// 배열 복제
		copy.modRange = copyMatrix(this.modRange);
		copy.ntermModRange = copyMatrix(this.ntermModRange);
		copy.ctermModRange = copyMatrix(this.ctermModRange);

		return copy;
	}

	private double[][] copyMatrix(double[][] original) {
		double[][] newMatrix = new double[original.length][original[0].length];
		for (int i = 0; i < original.length; i++) {
			System.arraycopy(original[i], 0, newMatrix[i], 0, original[i].length);
		}
		return newMatrix;
	}

	public void constructPTMLookupTable(){
		constructPTMTable();
	}

	private int constructPTMTable() {
		this.sortByPTMPosition();

		PTMTable = new ArrayList[AminoAcid.getIndexSize()][PTMPosition.PTMPOSITION_COUNT.ordinal()];
		for (int i=0; i<AminoAcid.getIndexSize(); i++) {
			for (int j=0; j<PTMPosition.PTMPOSITION_COUNT.ordinal(); j++)
				PTMTable[i][j] = new ArrayList<PTM>();
		}

		for ( PTM ptm : this ) {
			if( ptm.getMassDifference() < Constants.minModifiedMass || ptm.getMassDifference() > Constants.maxModifiedMass ) continue;
			if ( ptm.getResidue() != null ) {
				int aa = ptm.getResidue().getIndex();
				PTMTable[ptm.getResidue().getIndex()][ptm.getPTMPosition().ordinal()].add(ptm);

				if( ptm.getPTMPosition() == PTMPosition.ANYWHERE ){
					if( modRange[0][aa] > ptm.getMassDifference() ) modRange[0][aa] = ptm.getMassDifference(); // Min
					else if( modRange[1][aa] < ptm.getMassDifference() ) modRange[1][aa] = ptm.getMassDifference(); // Max
				}
				else if( ptm.getPTMPosition() == PTMPosition.ANY_N_TERM || ptm.getPTMPosition() == PTMPosition.PROTEIN_N_TERM ){
					if( ntermModRange[0][aa] > ptm.getMassDifference() ) ntermModRange[0][aa] = ptm.getMassDifference(); // Min
					else if( ntermModRange[1][aa] < ptm.getMassDifference() ) ntermModRange[1][aa] = ptm.getMassDifference(); // Max
				}
				else if( ptm.getPTMPosition() == PTMPosition.ANY_C_TERM || ptm.getPTMPosition() == PTMPosition.PROTEIN_C_TERM ){
					if( ctermModRange[0][aa] > ptm.getMassDifference() ) ctermModRange[0][aa] = ptm.getMassDifference(); // Min
					else if( ctermModRange[1][aa] < ptm.getMassDifference() ) ctermModRange[1][aa] = ptm.getMassDifference(); // Max
				}
			}
			else {
				if ( !ptm.isResidueCTerm() && !ptm.isResidueNTerm() ) return -1;

				for (int j=0; j<AminoAcid.getIndexSize(); j++){

					if( PTMTable[j][PTMPosition.ANYWHERE.ordinal()].contains(ptm) ) continue;

					if ( ptm.getPTMPosition() == PTMPosition.PROTEIN_N_TERM ) {
						if( PTMTable[j][PTMPosition.ANY_N_TERM.ordinal()].contains(ptm) )
							continue;
					}
					if ( ptm.getPTMPosition() == PTMPosition.PROTEIN_C_TERM ) {
						if( PTMTable[j][PTMPosition.ANY_C_TERM.ordinal()].contains(ptm) )
							continue;
					}//*/

					PTMTable[j][ptm.getPTMPosition().ordinal()].add(ptm);

					if( ptm.getPTMPosition() == PTMPosition.ANY_N_TERM || ptm.getPTMPosition() == PTMPosition.PROTEIN_N_TERM ){
						if( ntermModRange[0][j] > ptm.getMassDifference() ) ntermModRange[0][j] = ptm.getMassDifference(); // Min
						else if( ntermModRange[1][j] < ptm.getMassDifference() ) ntermModRange[1][j] = ptm.getMassDifference(); // Max
					}
					else if( ptm.getPTMPosition() == PTMPosition.ANY_C_TERM || ptm.getPTMPosition() == PTMPosition.PROTEIN_C_TERM ){
						if( ctermModRange[0][j] > ptm.getMassDifference() ) ctermModRange[0][j] = ptm.getMassDifference(); // Min
						else if( ctermModRange[1][j] < ptm.getMassDifference() ) ctermModRange[1][j] = ptm.getMassDifference(); // Max
					}
				}
			}
		}

		Collections.sort(this);
		return 0;
	}

	public ArrayList<PTM>[][] getPTMTable() {
		return this.PTMTable;
	}


	public PTMSearchResult searchPTM(Sequence seq, double massDiff, PTMPosition position) {
		// 슬롯별 상수 사용 (기존 구조 유지)
		final int slotIdx = ThreadPoolManager.getSlotIndex();

		ArrayList<PTMRun> out = new ArrayList<>();
		double err = Math.abs(massDiff);

		// (1) 무변형(비수정) 허용 오차 내부 → 바로 성공 표시
		if (err < Constants.nonModifiedDelta[slotIdx]) {
			return new PTMSearchResult(out, true);
		}

		// (2) 갭 허용오차 내부 → 빈 PTMRun(에러만 세팅) 하나 넣어둔다
		if (err < Constants.gapTolerance[slotIdx]) {
			PTMRun run = new PTMRun();
			run.setError(err);
			out.add(run);
		}

		// (3) DFS 준비: 전부 지역 상태로 관리 (스레드 안전)
		final PTM[] occur = new PTM[seq.size()];
		final int[] numNextFixSite = new int[seq.size()];
		final int numMaxMods = Constants.getMaxPTMOccurrence(seq.size());

		// 다음 위치들 중 fix site 개수 누적 (C-말단 고정 변형 고려 규칙 동일)
		int cum = (Constants.CTERM_FIX_MOD != 0 &&
			(position == PTMPosition.ANY_C_TERM || position == PTMPosition.PROTEIN_C_TERM)) ? 1 : 0;

		for (int i = seq.size() - 1; i >= 0; i--) {
			numNextFixSite[i] = seq.get(i).isLabelled() ? ++cum : cum;
		}

		// (4) 재귀 DFS 실행 — 원본 findPTM_DFS와 동치
		dfs(
			out,
			seq, massDiff, position,
			occur, numNextFixSite, numMaxMods,
			/*mass*/ 0.0, /*pos*/ 0, /*cnt*/ 0, /*extra*/ 0
		);

		return new PTMSearchResult(out, !out.isEmpty());
	}

	private void dfs(
		ArrayList<PTMRun> result,
		Sequence seq, double massDiff, PTMPosition position,
		PTM[] occur, int[] numNextFixSite, int numMaxMods,
		double mass, int pos, int cnt, int extra
	) {
		final int slotIdx = ThreadPoolManager.getSlotIndex();

		// 리프 도달: 질량 차이 검사
		if (pos == seq.size()) {
			double ierror = Math.abs(mass - massDiff);
			if (ierror <= Constants.gapTolerance[slotIdx] && Constants.isWithinAccuracy(ierror)) {
				PTMRun run = new PTMRun();
				for (int i = 0; i < seq.size(); i++) {
					if (occur[i] != null) {
						run.add(new PTMOccurrence(i, occur[i]));
					}
				}
				// 비수정 케이스는 searchPTM 상단에서 이미 처리했으므로 size>0만 추가
				if (run.size() > 0) {
					run.setError(ierror);
					result.add(run);
				}
			}
			return;
		}

		// (A) 현재 위치 PTM 적용 안 함
		dfs(result, seq, massDiff, position, occur, numNextFixSite, numMaxMods, mass, pos + 1, cnt, extra);

		// (B) 최대 갯수 제한
		if (cnt >= Constants.maxPTMPerGap) return;

		// (C) 고정 변형 1개만 허용되는 상황에서 뒤에 fix site 없으면 더 볼 필요 없음
		if (extra == 1 && numMaxMods == 1) {
			if (numNextFixSite[pos] == 0) return;
		}

		// (D) 현재 위치에 올 수 있는 PTM들 시도
		int residueIndex = seq.get(pos).getIndex();

		// 위치 종류 제약(원본과 동일)
		for (int i = 1; i < PTMPosition.PTMPOSITION_COUNT.ordinal(); i++) {
			// 내부 위치는 i==1만 허용
			if (pos > 0 && pos < seq.size() - 1 && i != 1) continue;
			// N-term 전용
			if ((i == 2 || i == 4) && pos != 0) continue;
			// C-term 전용
			if ((i == 3 || i == 5) && pos != seq.size() - 1) continue;

			// 사용자가 요청한 position 제약
			if ((i == 2) && position != PTMPosition.ANY_N_TERM && position != PTMPosition.PROTEIN_N_TERM) continue;
			if ((i == 3) && position != PTMPosition.ANY_C_TERM && position != PTMPosition.PROTEIN_C_TERM) continue;
			if ((i == 4) && position != PTMPosition.PROTEIN_N_TERM) continue;
			if ((i == 5) && position != PTMPosition.PROTEIN_C_TERM) continue;

			ArrayList<PTM> list = PTMTable[residueIndex][i];
			if (list == null || list.isEmpty()) continue;

			for (PTM ptm : list) {
				// 변형 개수 초과 시 skip
				int nextExtra = extra + ptm.getModCount();
				if (nextExtra > numMaxMods) continue;

				occur[pos] = ptm;
				dfs(result, seq, massDiff, position, occur, numNextFixSite, numMaxMods,
					mass + ptm.getMassDifference(), pos + 1, cnt + 1, nextExtra);
				occur[pos] = null; // 백트래킹
			}
		}
	}


	public void sortByPTMPosition(){
		Collections.sort( this, new PTMPosComparator() );
	}

	public int setVariableModificatinos( String fileName, double[] fixedAA, boolean canBeModifiedOnFixedAA ) throws Exception
	{
		int id = this.size();
		org.jdom.Document document = null;
		try {
			 document = new org.jdom.input.SAXBuilder().build(new File(fileName));
		} catch(org.jdom.JDOMException e) {
			System.out.println(e);
		}

		org.jdom.Element rootElement = document.getRootElement();
		org.jdom.Element classificationsElement = rootElement.getChild("classifications");
		if(classificationsElement != null)	// for backward compatibility
		{
			classifications.clear();
			for(Object obj : classificationsElement.getChildren())
			{
				org.jdom.Element classificationRowElement = (org.jdom.Element)obj;
				classifications.put(Integer.parseInt(classificationRowElement.getAttributeValue("recordID")),
						classificationRowElement.getAttributeValue("classification"));
			}
		}

		for(Object obj : rootElement.getChildren("PTM"))
		{
			org.jdom.Element elemPTM = (org.jdom.Element)obj;

			String pname = elemPTM.getChildText("name");

			String category = elemPTM.getChildText("classification");

			String residueStr = elemPTM.getChildText("residue");
			AminoAcid residue = null;

			double massdelta = 0;
			if(residueStr != null && residueStr.compareToIgnoreCase("N-term") != 0 && residueStr.compareToIgnoreCase("C-term") != 0) {
				assert( residueStr.length() == 1 );
				residue = AminoAcid.getAminoAcid(residueStr.charAt(0));
				if( fixedAA[residueStr.charAt(0)-'A'] != 0 ) {
					if( canBeModifiedOnFixedAA ) massdelta = fixedAA[residueStr.charAt(0)-'A'];
					else continue;
				}
				if( "AA substitution".compareTo(category) == 0 ){
					char tarAA = AminoAcid.getAminoAcidBy3Letter(pname.substring(pname.length()-3));
					massdelta -= fixedAA[tarAA-'A'];
					if( tarAA == 'I' && fixedAA['I'-'A'] != fixedAA['L'-'A'] ){
						System.out.println("PROGRAM SHOULD BE FIXED");
						System.exit(1);
					}
				}
			}

			PTMPosition position = null;
			String positionStr = elemPTM.getChildText("position");
			if(positionStr.equalsIgnoreCase("ANYWHERE"))
				position = PTMPosition.ANYWHERE;
			else if(positionStr.equalsIgnoreCase("ANY_N_TERM"))
				position = PTMPosition.ANY_N_TERM;
			else if(positionStr.equalsIgnoreCase("ANY_C_TERM"))
				position = PTMPosition.ANY_C_TERM;
			else if(positionStr.equalsIgnoreCase("PROTEIN_N_TERM"))
				position = PTMPosition.PROTEIN_N_TERM;
			else if(positionStr.equalsIgnoreCase("PROTEIN_C_TERM"))
				position = PTMPosition.PROTEIN_C_TERM;
			else
				assert(false);

		//	if( Constants.NTERM_FIX_MOD != 0 && (position == PTMPosition.ANY_N_TERM || position == PTMPosition.PROTEIN_N_TERM) ) continue;
		//	if( Constants.CTERM_FIX_MOD != 0 && (position == PTMPosition.ANY_C_TERM || position == PTMPosition.PROTEIN_C_TERM) ) continue;

			if( Constants.NTERM_FIX_MOD != 0 && (position == PTMPosition.ANY_N_TERM || position == PTMPosition.PROTEIN_N_TERM) ) {
				if( canBeModifiedOnFixedAA ) {
					pname += "/Nterm";
					massdelta += Constants.NTERM_FIX_MOD;
				}
				else continue;
			}
			if( Constants.CTERM_FIX_MOD != 0 && (position == PTMPosition.ANY_C_TERM || position == PTMPosition.PROTEIN_C_TERM) ) {
				if( canBeModifiedOnFixedAA ) {
					massdelta += Constants.CTERM_FIX_MOD;
					pname += "/Cterm";
				}
				else continue;
			}

			double ac_delta = Double.parseDouble(elemPTM.getChildText("massDifference"));
			if( ac_delta < Constants.minModifiedMass || ac_delta > Constants.maxModifiedMass ) continue;
			ac_delta -= massdelta;

			this.add(new PTM(
						id++,
						pname,
						"",
						ac_delta,
						0,
						residue,
						position,
						category
						)
					);
		}
		return 1;
	}

	public int setVariableModificatinos( org.jdom.Element modifications, double[] fixedAA, boolean canBeModifiedOnFixedAA ) throws Exception
	{
		int id = this.size();
		for(Object obj : modifications.getChildren("mod"))
		{
			org.jdom.Element elemPTM = (org.jdom.Element)obj;

			String residueStr = elemPTM.getAttributeValue("site");
			AminoAcid residue = null;
			double massdelta = 0;
			if(residueStr != null && residueStr.compareToIgnoreCase("N-term") != 0 && residueStr.compareToIgnoreCase("C-term") != 0) {
				assert(residueStr.length() == 1);
				residue = AminoAcid.getAminoAcid(residueStr.charAt(0));
				if( fixedAA[residueStr.charAt(0)-'A'] != 0 ) {
					if( canBeModifiedOnFixedAA ) massdelta = fixedAA[residueStr.charAt(0)-'A'];
					else continue;
				}
			}

			PTMPosition position = null;
			String positionStr = elemPTM.getAttributeValue("position");
			if(positionStr.equalsIgnoreCase("ANYWHERE"))
				position = PTMPosition.ANYWHERE;
			else if(positionStr.equalsIgnoreCase("ANY_N_TERM"))
				position = PTMPosition.ANY_N_TERM;
			else if(positionStr.equalsIgnoreCase("ANY_C_TERM"))
				position = PTMPosition.ANY_C_TERM;
			else if(positionStr.equalsIgnoreCase("PROTEIN_N_TERM"))
				position = PTMPosition.PROTEIN_N_TERM;
			else if(positionStr.equalsIgnoreCase("PROTEIN_C_TERM"))
				position = PTMPosition.PROTEIN_C_TERM;
			else
				assert(false);

			String pname = elemPTM.getAttributeValue("name");

			if( Constants.NTERM_FIX_MOD != 0 && (position == PTMPosition.ANY_N_TERM || position == PTMPosition.PROTEIN_N_TERM) ) {
				if( canBeModifiedOnFixedAA ) {
					pname += "/Nterm";
					massdelta += Constants.NTERM_FIX_MOD;
				}
				else continue;
			}
			if( Constants.CTERM_FIX_MOD != 0 && (position == PTMPosition.ANY_C_TERM || position == PTMPosition.PROTEIN_C_TERM) ) {
				if( canBeModifiedOnFixedAA ) {
					pname += "/Cterm";
					massdelta += Constants.CTERM_FIX_MOD;
				}
				else continue;
			}

			double ac_delta = Double.parseDouble(elemPTM.getAttributeValue("massdiff"));
			if( ac_delta < Constants.minModifiedMass || ac_delta > Constants.maxModifiedMass ) continue;
			ac_delta -= massdelta;

			this.add(new PTM(
						id++,
						pname,
						"",
						ac_delta,
						0,
						residue,
						position
						)
					);
		}
		return 1;//constructPTMTable();
	}

	public int setFixedModificatinos(org.jdom.Element modifications, double[] fixedAA) throws Exception
	{
		int id=0;
		for(Object obj : modifications.getChildren("mod"))
		{
			org.jdom.Element elemPTM = (org.jdom.Element)obj;

			String residueStr = elemPTM.getAttributeValue("site");
			AminoAcid residue = null;

			if( residueStr == null ) continue;
			else if( residueStr.compareToIgnoreCase("N-term") == 0 ){
				Constants.NTERM_FIX_MOD += Double.parseDouble(elemPTM.getAttributeValue("massdiff"));
			}
			else if( residueStr.compareToIgnoreCase("C-term") == 0 ){
				Constants.CTERM_FIX_MOD += Double.parseDouble(elemPTM.getAttributeValue("massdiff"));
			}
			else residue = AminoAcid.getAminoAcid(residueStr.charAt(0));

			double diff = Double.parseDouble(elemPTM.getAttributeValue("massdiff"));

			PTMPosition position = null;
			String positionStr = elemPTM.getAttributeValue("position");
			if(positionStr.equalsIgnoreCase("ANYWHERE"))
				position = PTMPosition.ANYWHERE;
			else if(positionStr.equalsIgnoreCase("ANY_N_TERM"))
				position = PTMPosition.ANY_N_TERM;
			else if(positionStr.equalsIgnoreCase("ANY_C_TERM"))
				position = PTMPosition.ANY_C_TERM;
			else if(positionStr.equalsIgnoreCase("PROTEIN_N_TERM"))
				position = PTMPosition.PROTEIN_N_TERM;
			else if(positionStr.equalsIgnoreCase("PROTEIN_C_TERM"))
				position = PTMPosition.PROTEIN_C_TERM;
			else
				assert(false);

			this.add(new PTM(
						id++,
						elemPTM.getAttributeValue("name"),
						"",
						diff,
						diff,
						residue,
						position
						)
					);
			if( residue != null ){
				if( fixedAA[residueStr.charAt(0)-'A'] != 0 ) return 0;
				AminoAcid.modifiedAminoAcidMass( residueStr.charAt(0), diff );
				MSMass.modifiedAminoAcidMass( residueStr.charAt(0), diff );
				fixedAA[residueStr.charAt(0)-'A']= diff;
			}
		}
		return 1;
	}

	public double minimumModifiedMass( Sequence seq, int pos ){
		double min1=0, min2=0, extra = 0;
		for( AminoAcid aa : seq ){
			int aix= aa.getIndex();
			if( min1 > modRange[0][aix] ) {
				min2 = min1;
				min1 = modRange[0][aix];
			}
			else if( min2 > modRange[0][aix] ) min2 = modRange[0][aix];
		}

		if( pos == 1 ){
			int aix= seq.get(0).getIndex();
			if( min1 > ntermModRange[0][aix] ) {
				min2 = min1;
				min1 = ntermModRange[0][aix];
			}
			else if( min2 > ntermModRange[0][aix] ) min2 = ntermModRange[0][aix];
		}
		else if( pos == 2 ){
			int aix= seq.get(seq.size()-1).getIndex();
			if( min1 > ctermModRange[0][aix] ) {
				min2 = min1;
				min1 = ctermModRange[0][aix];
			}
			else if( min2 > ctermModRange[0][aix] ) min2 = ctermModRange[0][aix];
		}
		if( seq.size() == 1 ) return min1;
		return min1 + min2;
	}

	public double maximumModifiedMass(Sequence seq, int pos){
		double max1=0, max2=0;
		for( AminoAcid aa : seq ){
			int aix= aa.getIndex();
			if( max1 < modRange[1][aix] ) {
				max2 = max1;
				max1 = modRange[1][aix];
			}
			else if( max2 < modRange[1][aix] ) max2 = modRange[1][aix];
		}
		if( pos == 1 ){
			int aix= seq.get(0).getIndex();
			if( max1 < ntermModRange[1][aix] ) {
				max2 = max1;
				max1 = ntermModRange[1][aix];
			}
			else if( max2 < ntermModRange[1][aix] ) max2 = ntermModRange[1][aix];
		}
		else if( pos == 2 ){
			int aix= seq.get(seq.size()-1).getIndex();
			if( max1 < ctermModRange[1][aix] ) {
				max2 = max1;
				max1 = ctermModRange[1][aix];
			}
			else if( max2 < ctermModRange[1][aix] ) max2 = ctermModRange[1][aix];
		}
		if( seq.size() == 1 ) return max1;
		return max1 + max2;
	}

	public void setPTMDiagnosticIon() {
		for (PTM ptm : this) {

			//setting common modifications
			if( Math.abs(15.994915-ptm.getMassDifference()) < 0.01 ){
				if( ptm.getAbbAA() =='M' ) {
					ptm.setPenalty(0.5);
					ptm.setModCount(0);
					AminoAcid.canBeEasilyModified('M');
					ptm.setNeutralLoss( ptm.getMassDifference()+48.003371 );
				}
				if( ptm.getAbbAA() =='W' ) ptm.setPenalty(0.8);
			}

			if( ptm.getAbbAA() =='Q' && ptm.getPTMPosition() == PTMPosition.ANY_N_TERM && Constants.round(ptm.getMassDifference()+Constants.NTERM_FIX_MOD) == -17 ){
				ptm.setPenalty(0.5);
				ptm.setModCount(0);
			}

			//setting isobaric labeling
			if( Constants.reporterMassOfIsobaricTag != null && Math.abs(Constants.reporterMassOfIsobaricTag[0]-ptm.getMassDifference()) < 0.01 ){
				if( ptm.getAbbAA() == 'S' || ptm.getAbbAA() == 'T' ) {
					if( "itraq4plex".compareTo(Constants.isobaricTag) == 0 ) ptm.setDiagnosticIon(163.1199);
					ptm.setPenalty(0.5);
					ptm.setModCount(0);
					AminoAcid.canBeEasilyModified(ptm.getAbbAA());
				}
			}

			//setting specific modifications
			if( ptm.getName().startsWith("Acetyl") ) {
				if( ptm.getAbbAA() =='K' ) ptm.setDiagnosticIon(126.0913);
				if( ptm.getPTMPosition() == PTMPosition.PROTEIN_N_TERM ) ptm.setPenalty(0.5);

				if( "Acetyl".compareToIgnoreCase(Constants.enrichedModification) == 0 ) {
					ptm.setPenalty(0.5);
					if( ptm.getAbbAA() =='K' ) {
						ptm.setModCount(0);
						AminoAcid.canBeEasilyModified('K');
					}
					if( ptm.getPTMPosition() == PTMPosition.PROTEIN_N_TERM ) {
						ptm.setModCount(0);
					}
				}
			}

			if( "Phospho".compareTo(ptm.getName()) == 0 ){
				if( ptm.getAbbAA() =='S' || ptm.getAbbAA() =='T' ) {
					ptm.setNeutralLoss(ptm.getMassDifference()+Constants.H2O);
					if( "Phospho".compareToIgnoreCase(Constants.enrichedModification) == 0 ) {
						ptm.setPenalty(0.5);
						ptm.setModCount(0);
						AminoAcid.canBeEasilyModified(ptm.getAbbAA());
					}
				}
			}
		}
	}

}


