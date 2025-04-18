package modi;

import java.util.ArrayList;
import java.util.Stack;
import moda.ThreadPoolManager;

public class PTMRunBuilder {

    public static void runDFS(
            PTMDB db,
            Sequence seq,
            double massDiff,
            PTMPosition position,
            PTM[] occur,
            int[] numNextFixSite,
            int numMaxMods,
            ArrayList<PTMRun> result
    ) {
        class State {
            double mass;
            int pos, cnt, extra;
            PTM ptmAtPos;

            State(double mass, int pos, int cnt, int extra, PTM ptmAtPos) {
                this.mass = mass;
                this.pos = pos;
                this.cnt = cnt;
                this.extra = extra;
                this.ptmAtPos = ptmAtPos;
            }
        }

        Stack<State> stack = new Stack<>();
        stack.push(new State(0.0, 0, 0, 0, null));

        int slotIdx = ThreadPoolManager.getSlotIndex();

        while (!stack.isEmpty()) {
            State s = stack.pop();

            if (s.pos > seq.size()) continue;

            if (s.pos > 0 && s.pos - 1 < occur.length) {
                occur[s.pos - 1] = s.ptmAtPos;
            }

            if (s.pos == seq.size()) {
                double ierror = Math.abs(s.mass - massDiff);
                if (ierror <= Constants.gapTolerance[slotIdx] && Constants.isWithinAccuracy(ierror)) {
                    PTMRun run = new PTMRun();
                    for (int i = 0; i < seq.size(); i++) {
                        if (occur[i] != null) {
                            run.add(new PTMOccurrence(i, occur[i]));
                        }
                    }
                    if (run.size() > 0) {
                        run.setError(ierror);
                        result.add(run);
                    }
                }
                continue;
            }

            if (s.pos >= seq.size() || seq.get(s.pos) == null) {
                System.err.println("NULL AminoAcid at seq[" + s.pos + "], skipping...");
                continue;
            }

            stack.push(new State(s.mass, s.pos + 1, s.cnt, s.extra, null));

            if (s.cnt >= Constants.maxPTMPerGap) continue;
            if (s.extra == 1 && numMaxMods == 1 && (s.pos >= numNextFixSite.length || numNextFixSite[s.pos] == 0)) continue;

            int residueIndex = seq.get(s.pos).getIndex();
            if (residueIndex >= db.getPTMTable().length || db.getPTMTable()[residueIndex] == null) continue;

            for (int i = 1; i < PTMPosition.PTMPOSITION_COUNT.ordinal(); i++) {
                if (s.pos > 0 && s.pos < seq.size() - 1 && i != 1) continue;
                if ((i == 2 || i == 4) && s.pos != 0) continue;
                if ((i == 3 || i == 5) && s.pos != seq.size() - 1) continue;

                if ((i == 2) && position != PTMPosition.ANY_N_TERM && position != PTMPosition.PROTEIN_N_TERM) continue;
                if ((i == 3) && position != PTMPosition.ANY_C_TERM && position != PTMPosition.PROTEIN_C_TERM) continue;
                if ((i == 4) && position != PTMPosition.PROTEIN_N_TERM) continue;
                if ((i == 5) && position != PTMPosition.PROTEIN_C_TERM) continue;

                ArrayList<PTM> ptms = db.getPTMTable()[residueIndex][i];
                if (ptms == null) continue;

                for (PTM ptm : ptms) {
                    if (ptm == null) continue;
                    if (s.extra + ptm.getModCount() > numMaxMods) continue;

                    stack.push(new State(
                            s.mass + ptm.getMassDifference(),
                            s.pos + 1,
                            s.cnt + 1,
                            s.extra + ptm.getModCount(),
                            ptm
                    ));
                }
            }

            if (s.pos < occur.length) {
                occur[s.pos] = null;
            }
        }
    }
}
