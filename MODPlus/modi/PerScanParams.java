package modi;

/**
 * 스캔마다 달라지는 파라미터(불변) 묶음
 */
public final class PerScanParams {
	public final double precursorTolerance;   // e.g. 0.5
	public final double precursorAccuracy;    // e.g. 0.5
	public final double fragmentTolerance;    // e.g. 0.5
	public final double gapTolerance;         // e.g. 0.6
	public final double gapAccuracy;          // e.g. 1.6
	public final double nonModifiedDelta;     // de-novo 등에서 쓰던 massToleranceForDenovo
	public final int    maxNoOfC13;           // e.g. 0

	public PerScanParams(
		double precursorTolerance,
		double precursorAccuracy,
		double fragmentTolerance,
		double gapTolerance,
		double gapAccuracy,
		double nonModifiedDelta,
		int maxNoOfC13
	) {
		this.precursorTolerance = precursorTolerance;
		this.precursorAccuracy  = precursorAccuracy;
		this.fragmentTolerance  = fragmentTolerance;
		this.gapTolerance       = gapTolerance;
		this.gapAccuracy        = gapAccuracy;
		this.nonModifiedDelta   = nonModifiedDelta;
		this.maxNoOfC13         = maxNoOfC13;
	}
}
