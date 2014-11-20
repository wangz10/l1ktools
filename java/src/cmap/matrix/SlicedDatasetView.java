/**
 * Copyright 2009-2010 by The Broad Institute.
 * 
 
 */
package cmap.matrix;

import java.util.Collection;

/**
 * Views a slice of a {@link Dataset}
 */
public class SlicedDatasetView implements Dataset {
    protected int[] columnIndices;

    protected Dataset dataset;

    protected int[] rowIndices;

    private MetadataModel columnMetadata;

    private MetadataModel rowMetadata;

    public SlicedDatasetView(Dataset dataset, int[] rowIndices, int[] columnIndices) {
	this.dataset = dataset;
	setIndices(rowIndices, columnIndices);
    }

    @Override
    public int getColumnCount() {
	return columnIndices == null ? dataset.getColumnCount() : columnIndices.length;
    }

    public int[] getColumnIndices() {
	return this.columnIndices;
    }

    @Override
    public MetadataModel getColumnMetadata() {
	return columnMetadata;
    }

    public Dataset getDataset() {
	return dataset;
    }

    @Override
    public String getName() {
	return dataset.getName();
    }

    @Override
    public Object getObjectValue(int row, int column, int seriesIndex) {
	return dataset.getObjectValue(convertRowIndexToModel(row), convertColumnIndexToModel(column), seriesIndex);
    }

    @Override
    public Object getProperty(Object key) {
	return dataset.getProperty(key);
    }

    @Override
    public Collection<Object> getPropertyKeys() {
	return dataset.getPropertyKeys();
    }

    @Override
    public int getRowCount() {
	return rowIndices == null ? dataset.getRowCount() : rowIndices.length;
    }

    public int[] getRowIndices() {
	return this.rowIndices;
    }

    @Override
    public MetadataModel getRowMetadata() {
	return rowMetadata;
    }

    @Override
    public Class<?> getSeriesClass(int index) {
	return dataset.getSeriesClass(index);
    }

    @Override
    public int getSeriesCount() {
	return dataset.getSeriesCount();
    }

    @Override
    public String getSeriesName(int index) {
	return dataset.getSeriesName(index);
    }

    @Override
    public float getValue(int rowIndex, int columnIndex) {
	return dataset.getValue(convertRowIndexToModel(rowIndex), convertColumnIndexToModel(columnIndex));
    }

    @Override
    public float getValue(int row, int column, int seriesIndex) {
	return dataset.getValue(convertRowIndexToModel(row), convertColumnIndexToModel(column), seriesIndex);
    }

    public void setColumnIndices(int[] columnIndices) {
	this.columnIndices = columnIndices;
	this.columnMetadata = columnIndices == null ? dataset.getColumnMetadata() : new MetadataModelItemView(
		dataset.getColumnMetadata(), columnIndices);
    }

    public void setIndices(int[] rowIndices, int[] columnIndices) {
	setRowIndices(rowIndices);
	setColumnIndices(columnIndices);
    }

    @Override
    public void setName(String name) {
	dataset.setName(name);
    }

    @Override
    public void setProperty(Object key, Object value) {
	dataset.setProperty(key, value);
    }

    public void setRowIndices(int[] rowIndices) {
	this.rowIndices = rowIndices;
	this.rowMetadata = rowIndices == null ? dataset.getRowMetadata() : new MetadataModelItemView(
		dataset.getRowMetadata(), rowIndices);
    }

    @Override
    public void setValue(int rowIndex, int columnIndex, float value) {
	dataset.setValue(convertRowIndexToModel(rowIndex), convertColumnIndexToModel(columnIndex), value);
    }

    private final int convertColumnIndexToModel(int column) {
	return columnIndices == null ? column : columnIndices[column];
    }

    private final int convertRowIndexToModel(int row) {
	return rowIndices == null ? row : rowIndices[row];
    }

}
