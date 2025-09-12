package modi;

/**
 * ThreadLocal에 현재 스캔의 PerScanParams를 보관.
 * try-with-resources로 push/pop.
 */
public final class ScanContext implements AutoCloseable {
	private static final ThreadLocal<PerScanParams> TL = new ThreadLocal<>();

	private ScanContext() {}

	public static ScanContext push(PerScanParams p) {
		TL.set(p);
		return new ScanContext();
	}

	public static PerScanParams current() {
		PerScanParams p = TL.get();
		if (p == null) throw new IllegalStateException("ScanContext not bound for current thread");
		return p;
	}

	@Override
	public void close() {
		TL.remove();
	}
}
