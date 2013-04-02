package uk.ac.ox.well.indiana.utils.performance;

public class PerformanceUtils {
    private PerformanceUtils() {}

    public static String getMemoryUsageStats() {
        int mb = 1024*1024;

        Runtime runtime = Runtime.getRuntime();

        long usedMemory = (runtime.totalMemory() - runtime.freeMemory()) / mb;
        long freeMemory = runtime.freeMemory() / mb;
        long totalMemory = runtime.totalMemory() / mb;
        long maxMemory = runtime.maxMemory() / mb;

        return String.format("Used: %dmb; Free: %dmb; Total: %dmb; Max: %dmb", usedMemory, freeMemory, totalMemory, maxMemory);
    }

    public static String getCompactMemoryUsageStats() {
        int mb = 1024*1024;

        Runtime runtime = Runtime.getRuntime();

        long usedMemory = (runtime.totalMemory() - runtime.freeMemory()) / mb;
        long freeMemory = runtime.freeMemory() / mb;
        long totalMemory = runtime.totalMemory() / mb;
        long maxMemory = runtime.maxMemory() / mb;

        return String.format("U:%dmb; F:%dmb; T:%dmb; M:%dmb", usedMemory, freeMemory, totalMemory, maxMemory);
    }
}
