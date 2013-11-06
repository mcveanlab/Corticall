package uk.ac.ox.well.indiana.utils.containers;

import java.util.*;

public class ApproximateHashSet<E> implements Set<E> {
    private Set<Integer> records = new HashSet<Integer>();

    @Override
    public int size() { return records.size(); }

    @Override
    public boolean add(E e) {
        if (!records.contains(e.hashCode())) {
            records.add(e.hashCode());

            return true;
        }

        return false;
    }

    @Override
    public boolean isEmpty() {
        return records.isEmpty();
    }

    @Override
    public boolean contains(Object o) {
        return records.contains(o.hashCode());
    }

    @Override
    public Iterator<E> iterator() {
        throw new UnsupportedOperationException("iterator() not implemented for ApproximateHashSet");
    }

    @Override
    public Object[] toArray() {
        throw new UnsupportedOperationException("toArray() not implemented for ApproximateHashSet");
    }

    @Override
    public <T> T[] toArray(T[] a) {
        throw new UnsupportedOperationException("toArray() not implemented for ApproximateHashSet");
    }

    @Override
    public boolean remove(Object o) {
        return records.remove(o.hashCode());
    }

    @Override
    public boolean containsAll(Collection<?> c) {
        throw new UnsupportedOperationException("containsAll() not implemented for ApproximateHashSet");
    }

    @Override
    public boolean addAll(Collection<? extends E> c) {
        for (E e : c) {
            add(e);
        }

        return true;
    }

    @Override
    public boolean retainAll(Collection<?> c) {
        throw new UnsupportedOperationException("retainAll() not implemented for ApproximateHashSet");
    }

    @Override
    public boolean removeAll(Collection<?> c) {
        for (Object o : c) {
            remove(o);
        }

        return true;
    }

    @Override
    public void clear() {
        records.clear();
    }
}
