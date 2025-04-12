package moda;

import java.util.concurrent.atomic.AtomicInteger;

public class ThreadPoolManager {

    public static final int numSlots = Runtime.getRuntime().availableProcessors();
    // 전역 카운터 (모든 스레드에 대해 증가, 모듈러 연산을 통해 슬롯 인덱스 결정)
    private static final AtomicInteger counter = new AtomicInteger(0);

    private static final int MAX_SLOTS = 36;
    // ThreadLocal을 이용해 각 스레드마다 고유한 슬롯 인덱스를 할당
    private static final ThreadLocal<Integer> threadSlot = ThreadLocal.withInitial(() -> {
        // 카운터의 값을 MAX_SLOTS로 모듈러 연산
        return counter.getAndIncrement() % MAX_SLOTS;
    });

    /**
     * 현재 스레드에 할당된 슬롯 인덱스를 리턴합니다.
     */
    public static int getSlotIndex() {
        return threadSlot.get();
    }

}
