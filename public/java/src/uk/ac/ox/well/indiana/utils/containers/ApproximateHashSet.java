package uk.ac.ox.well.indiana.utils.containers;

import java.util.Collection;
import java.util.Iterator;
import java.util.Set;

// Inspired by a blog post at http://ivory.idyll.org/blog/kmer-filtering.html
public class ApproximateHashSet<E> implements Set<E> {
    private char[] kmers;

    public ApproximateHashSet() {
        long freeMemory = Runtime.getRuntime().freeMemory();
        int tableLength = (int) (0.8f * (float) freeMemory);

        initialize(tableLength);
    }

    public ApproximateHashSet(float memoryFraction) {
        long freeMemory = Runtime.getRuntime().freeMemory();
        int tableLength = (int) (memoryFraction * (float) freeMemory);

        initialize(tableLength);
    }

    public ApproximateHashSet(int tableLength) {
        initialize(tableLength);
    }

    private void initialize(int tableLength) {
        kmers = new char[tableLength];
    }

    @Override
    public int size() {
        int size = 0;
        for (int i = 0; i < kmers.length; i++) {
            size += kmers[i];
        }

        return size;
    }

    @Override
    public boolean isEmpty() {
        throw new UnsupportedOperationException("isEmpty() not implemented for ApproximateHashSet");
    }

    @Override
    public boolean contains(Object o) {
        throw new UnsupportedOperationException("contains() not implemented for ApproximateHashSet");
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
    public boolean add(E e) {
        int index = e.hashCode() & (kmers.length - 1);

        kmers[index]++;

        return true;
    }

    @Override
    public boolean remove(Object o) {
        throw new UnsupportedOperationException("remove() not implemented for ApproximateHashSet");
    }

    @Override
    public boolean containsAll(Collection<?> c) {
        throw new UnsupportedOperationException("containsAll() not implemented for ApproximateHashSet");
    }

    @Override
    public boolean addAll(Collection<? extends E> c) {
        throw new UnsupportedOperationException("addAll() not implemented for ApproximateHashSet");
    }

    @Override
    public boolean retainAll(Collection<?> c) {
        throw new UnsupportedOperationException("retainAll() not implemented for ApproximateHashSet");
    }

    @Override
    public boolean removeAll(Collection<?> c) {
        throw new UnsupportedOperationException("removeAll() not implemented for ApproximateHashSet");
    }

    @Override
    public void clear() {
        throw new UnsupportedOperationException("clear() not implemented for ApproximateHashSet");
    }
}
