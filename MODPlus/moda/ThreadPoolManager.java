package moda;

public final class ThreadPoolManager {
    public static final int numSlots = Runtime.getRuntime().availableProcessors();
    //public static final int numSlots = 1;

    private static final java.util.concurrent.atomic.AtomicInteger NEXT = new java.util.concurrent.atomic.AtomicInteger(0);
    private static final ThreadLocal<Integer> SLOT = new ThreadLocal<>();
    private static volatile boolean PARALLEL = false;

    private ThreadPoolManager() {}

    public static void enterParallelPhase() {
        PARALLEL = true;
        NEXT.set(0);
        SLOT.remove();
    }
    public static void exitParallelPhase() {
        PARALLEL = false;
        SLOT.remove();
    }
    public static void bindSlotForWorker() {
        final int assigned = NEXT.getAndIncrement();
        if (assigned >= numSlots)
            throw new IllegalStateException("More threads than slots: " + assigned + " >= " + numSlots);
        if (SLOT.get() != null)
            throw new IllegalStateException("Slot already bound on this thread");
        SLOT.set(assigned);
    }
    /** 병렬 구간 외에는 0을 반환(스펙 로딩 등 싱글 단계). 병렬 구간 중 미바인딩이면 예외 */
    public static int getSlotIndex() {
        Integer idx = SLOT.get();
        if (idx != null) return idx;
        if (!PARALLEL) return 0;
        throw new IllegalStateException("Slot not bound on this thread");
    }
}

