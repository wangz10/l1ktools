/**
 * Copyright 2009-2010 by The Broad Institute.
 * 
 
 */
package cmap.matrix;

import java.util.Collection;

/**
 * The interface for a dataset consisting of a two-dimensional matrix of float
 * values. A dataset may also optionally contain one or more series of
 * two-dimensional matrices. A dataset also has metadata associated with each
 * row and column.
 */
public interface Dataset {
    /**
     * Returns the number of columns in the dataset.
     * 
     * @return the number of columns
     */
    public int getColumnCount();

    /**
     * Gets the column metadata for this dataset.
     * 
     * @return The column metadata
     */
    public MetadataModel getColumnMetadata();

    /**
     * Returns the name of this dataset.
     * 
     * @return The dataset name
     */
    public String getName();

    /**
     * Returns the value at the given row and column for the given series.
     * Series can be used to store standard error of data points for example.
     * 
     * @param rowIndex
     *            the row index
     * @param columnIndex
     *            the column index
     * @param seriesIndex
     *            the series index
     * @return the value
     * @throws IndexOutOfBoundsException
     *             if <tt>rowIndex</tt> or <tt>columnIndex</tt> are out of
     *             range.
     */
    public Object getObjectValue(int rowIndex, int columnIndex, int seriesIndex);

    /**
     * Returns the value to which the specified key is mapped, or <tt>null</tt>
     * if this dataset does not contain the specified property.
     * 
     * @param key
     *            the key
     * @return the value the value for the specified key
     */
    public Object getProperty(Object key);

    /**
     * Returns a collection of all the property keys contained in this dataset.
     * 
     * @return the property keys
     */
    public Collection<Object> getPropertyKeys();

    /**
     * Returns the number of rows in the dataset.
     * 
     * @return the number of rows
     */
    public int getRowCount();

    /**
     * Gets the row metadata for this dataset.
     * 
     * @return the row metadata
     */
    public MetadataModel getRowMetadata();

    /**
     * Returns the most specific superclass for all the values in the series.
     * 
     * @param series
     *            the index of the series
     * @return the common ancestor class of the series values
     */
    public Class<?> getSeriesClass(int series);

    /**
     * Returns the number of matrix series, which is always at least 1 (the
     * dataset itself). Series can be used to store standard error of data
     * points for example.
     * 
     * @return the number of series
     */
    public int getSeriesCount();

    /**
     * Returns the name for the given series. Series can be used to store
     * standard error of data points for example.
     * 
     * @param series
     *            the series
     * @return the data series name
     * @throws IndexOutOfBoundsException
     *             if <tt>series</tt> is out of range
     *             <tt>(index &lt; 0 || index &gt;=
     *                                   getSeriesCount())</tt>
     */
    public String getSeriesName(int series);

    /**
     * Returns the value at the given row and column.
     * 
     * @param rowIndex
     *            the row index.
     * @param columnIndex
     *            the column index.
     * @return the value.
     * @throws IndexOutOfBoundsException
     *             if <tt>rowIndex</tt> or <tt>columnIndex</tt> are out of
     *             range.
     */
    public float getValue(int rowIndex, int columnIndex);

    /**
     * Returns the value at the given row and column for the given series as a
     * float. Series can be used to store standard error of data points for
     * example.
     * 
     * @param rowIndex
     *            the row index
     * @param columnIndex
     *            the column index
     * @param seriesIndex
     *            the series index
     * @return the value
     * @throws IndexOutOfBoundsException
     *             if <tt>rowIndex</tt> or <tt>columnIndex</tt> are out of range
     */
    public float getValue(int rowIndex, int columnIndex, int seriesIndex);

    /**
     * Sets the name of this dataset.
     * 
     * @throws UnsupportedOperationException
     *             if the <tt>setName</tt> operation is not supported by this
     *             Dataset
     * 
     * @param name
     *            The dataset name
     */
    public void setName(String name);

    /**
     * Sets the property for the specified key.
     * 
     * @param key
     *            the key
     * @param value
     *            the value
     */
    public void setProperty(Object key, Object value);

    /**
     * Sets the value at the given row and column indices.
     * 
     * @param rowIndex
     *            the row index
     * 
     * @param columnIndex
     *            the column index
     * @param value
     *            the value
     * @throws UnsupportedOperationException
     *             if the <tt>setValue</tt> operation is not supported by this
     *             Dataset
     */
    public void setValue(int rowIndex, int columnIndex, float value);

}
