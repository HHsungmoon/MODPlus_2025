package modi;

public class ThreadPoolManager {
    public static final int MAX_THREAD = 2;

    // slot index를 각 스레드에 할당
    public static final ThreadLocal<Integer> SLOT_INDEX = new ThreadLocal<>();


    public static int getSlotIndex() {
        return SLOT_INDEX.get();
    }

}
