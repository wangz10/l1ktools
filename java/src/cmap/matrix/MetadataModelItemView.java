/**
 * Copyright 2009-2010 by The Broad Institute.
 * 
 
 */
package cmap.matrix;

public class MetadataModelItemView implements MetadataModel {
    private int[] indices;
    private MetadataModel metadata;

    public MetadataModelItemView(MetadataModel metaData, int[] indices) {
	this.metadata = metaData;
	this.indices = indices;
    }

    @Override
    public Vector add(String columnName, Class<?> columnClass) {
	Vector v = metadata.add(columnName, columnClass);
	return indices == null ? v : new SlicedVector(v, indices);
    }

    @Override
    public Vector get(int columnIndex) {
	return indices == null ? metadata.get(columnIndex) : new SlicedVector(
		metadata.get(columnIndex), indices);
    }

    @Override
    public Vector get(String name) {
	int index = metadata.getColumnIndex(name);
	if (index == -1) {
	    return null;
	}
	return indices == null ? metadata.get(index) : new SlicedVector(
		metadata.get(index), indices);
    }

    @Override
    public int getColumnIndex(String columnName) {
	return metadata.getColumnIndex(columnName);
    }

    @Override
    public String getColumnName(int columnIndex) {
	return metadata.getColumnName(columnIndex);
    }

    @Override
    public int getItemCount() {
	return indices != null ? indices.length : metadata.getItemCount();
    }

    @Override
    public int getMetadataCount() {
	return metadata.getMetadataCount();
    }

    @Override
    public Object getValue(int rowIndex, int columnIndex) {
	return metadata.getValue(getMappedRowIndex(rowIndex), columnIndex);
    }

    @Override
    public Object getValue(int rowIndex, String columnName) {
	return metadata.getValue(getMappedRowIndex(rowIndex), columnName);
    }

    @Override
    public void remove(int index) {
	metadata.remove(index);
    }

    @Override
    public void setValue(int rowIndex, int columnIndex, Object aValue) {
	metadata.setValue(getMappedRowIndex(rowIndex), columnIndex, aValue);
    }

    @Override
    public void setValue(int rowIndex, String columnName, Object value) {
	metadata.setValue(getMappedRowIndex(rowIndex), columnName, value);
    }

    private final int getMappedRowIndex(int row) {
	return indices != null ? indices[row] : row;
    }

}
