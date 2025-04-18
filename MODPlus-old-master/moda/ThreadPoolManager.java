package moda;


import java.util.concurrent.atomic.AtomicInteger;

public class ThreadPoolManager {

    // 실제 사용하는 스레드 슬롯 수
    public static final int numSlots = Runtime.getRuntime().availableProcessors();

    // 고유 슬롯 인덱스를 순차적으로 할당하기 위한 카운터
    private static final AtomicInteger counter = new AtomicInteger(0);

    // 각 스레드에 고유한 슬롯 인덱스를 부여하는 ThreadLocal
    private static final ThreadLocal<Integer> threadSlot = ThreadLocal.withInitial(() -> {
        int slot = counter.getAndIncrement() % numSlots;
        return slot;
    });

    // 현재 스레드에 할당된 슬롯 인덱스를 반환
    public static int getSlotIndex() {
        return threadSlot.get();
    }
}
