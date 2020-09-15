package ru.spb.arcadia.jnj.aam.utils;

import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;


public class ThreadSafeCache<K, V> implements Cache<K, V> {
    private final Map<K, V> map;

    private static final ThreadSafeCache SC = new ThreadSafeCache();

    public static ThreadSafeCache getInstance() {
        return SC;
    }

    private ThreadSafeCache() {
        map = new ConcurrentHashMap<>();
    }

    @Override
    public synchronized void put(K key, V value) {
        map.put(key, value);
    }

    @Override
    public synchronized V get(K key) {
        return map.get(key);
    }

    public synchronized boolean containsKey(K key) {
        return map.containsKey(key);
    }

    public synchronized void cleanup() {
        synchronized (map) {
            map.clear();
        }
        Thread.yield();
    }

    public synchronized int size() {
        return map.size();
    }

    public synchronized Set<K> keySet() {
        return map.keySet();
    }
}