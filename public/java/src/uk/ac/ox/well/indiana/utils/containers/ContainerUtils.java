package uk.ac.ox.well.indiana.utils.containers;

import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class ContainerUtils {
    private void ContainerUtils() {};

    public static void increment(Map<String, Integer> map, String key) {
        if (!map.containsKey(key)) {
            map.put(key, 1);
        } else {
            map.put(key, map.get(key) + 1);
        }
    }

    public static <K, V> void add(Map<K, Set<V>> map, K key, V value) {
        if (!map.containsKey(key)) {
            map.put(key, new HashSet<>());
        }
        map.get(key).add(value);
    }

    public static String mostCommonKey(Map<String, Integer> map) {
        String commonKey = map.keySet().iterator().next();
        int count = map.get(commonKey);

        for (String key : map.keySet()) {
            if (map.get(key) > count) {
                commonKey = key;
                count = map.get(key);
            }
        }

        return commonKey;
    }
}
