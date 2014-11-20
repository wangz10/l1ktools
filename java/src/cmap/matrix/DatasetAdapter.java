/**
 * Copyright 2009-2010 by The Broad Institute.
 * 
 
 */
package cmap.matrix;

import java.util.Collection;

public class DatasetAdapter implements Dataset {
    protected Dataset dataset;

    public DatasetAdapter(Dataset dataset) {
	this.dataset = dataset;
    }

    @Override
    public int getColumnCount() {
	return dataset.getColumnCount();
    }

    @Override
    public MetadataModel getColumnMetadata() {
	return dataset.getColumnMetadata();
    }

    public Dataset getDataset() {
	return dataset;
    }

    @Override
    public String getName() {
	return dataset.getName();
    }

    @Override
    public Object getObjectValue(int rowIndex, int columnIndex, int seriesIndex) {
	return dataset.getObjectValue(rowIndex, columnIndex, seriesIndex);
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
	return dataset.getRowCount();
    }

    @Override
    public MetadataModel getRowMetadata() {
	return dataset.getRowMetadata();
    }

    @Override
    public Class<?> getSeriesClass(int series) {
	return dataset.getSeriesClass(series);
    }

    @Override
    public int getSeriesCount() {
	return dataset.getSeriesCount();
    }

    @Override
    public String getSeriesName(int series) {
	return dataset.getSeriesName(series);
    }

    @Override
    public float getValue(int rowIndex, int columnIndex) {
	return dataset.getValue(rowIndex, columnIndex);
    }

    @Override
    public float getValue(int rowIndex, int columnIndex, int seriesIndex) {
	return dataset.getValue(rowIndex, columnIndex, seriesIndex);
    }

    @Override
    public void setName(String name) {
	dataset.setName(name);
    }

    @Override
    public void setProperty(Object key, Object value) {
	dataset.setProperty(key, value);
    }

    @Override
    public void setValue(int rowIndex, int columnINdex, float value) {
	dataset.setValue(rowIndex, columnINdex, value);
    }

}
