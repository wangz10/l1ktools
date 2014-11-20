package cmap.matrix;

import java.util.Collection;
import java.util.Comparator;

public class SlicedVector implements Vector {
    private int[] indices;

    private Vector vector;

    public SlicedVector(Vector vector, int[] indices) {
	this.vector = vector;
	if (indices == null) {
	    indices = new int[vector.size()];
	    for (int i = 0, length = indices.length; i < length; i++) {
		indices[i] = i;
	    }

	}
	this.indices = indices;
    }

    @Override
    public Class<?> getColumnClass() {
	return vector.getColumnClass();
    }

    @Override
    public Comparator<Object> getComparator() {
	return vector.getComparator();
    }

    public int[] getIndices() {
	return indices;
    }

    @Override
    public String getName() {
	return vector.getName();
    }

    @Override
    public Object getProperty(Object key) {
	return vector.getProperty(key);
    }

    @Override
    public Collection<Object> getPropertyKeys() {
	return vector.getPropertyKeys();
    }

    @Override
    public Object getValue(int index) {
	return vector.getValue(indices[index]);
    }

    public Vector getVector() {
	return vector;
    }

    @Override
    public void removeProperty(Object key) {
	vector.removeProperty(key);
    }

    @Override
    public void setColumnClass(Class<?> columnClass) {
	vector.setColumnClass(columnClass);
    }

    @Override
    public void setComparator(Comparator<Object> comparator) {
	vector.setComparator(comparator);
    }

    @Override
    public void setProperty(Object key, Object value) {
	vector.setProperty(key, value);
    }

    @Override
    public void setValue(int index, Object value) {
	vector.setValue(indices[index], value);
    }

    @Override
    public int size() {
	return indices.length;
    }

    @Override
    public String toString() {
	return vector.toString();
    }

}
