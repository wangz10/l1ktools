/**
 * Copyright 2009-2010 by The Broad Institute.
 *
 */
package cmap.matrix;

/**
 * Stores annotations for the rows or columns of a dataset.
 */
public interface MetadataModel {

    /**
     * Appends the specified column name to the end of this meta data (optional
     * operation).
     * 
     * @param columnName
     *            The column name to be inserted into this meta data instance.
     * @param columnClass
     *            The class of the column to be inserted.
     * @return the added vector.
     * @throws UnsupportedOperationException
     *             if the <tt>addColumn</tt> method is not supported by this
     *             meta data instance.
     */
    public Vector add(String columnName, Class<?> columnClass);

    /**
     * Returns the vector at the specified metadata index.
     * 
     * @param metadataIndex
     *            the metadata index
     * @return the vector
     */
    public Vector get(int metadataIndex);

    /**
     * Returns the vector witht the specified name.
     * 
     * @param name
     *            the vector name
     * @return the vector
     */
    public Vector get(String name);

    /**
     * Returns the column index for the given column name or <tt>-1</tt> if the
     * name is not found.
     * 
     * @param columnName
     *            the column name
     * @return the column index
     */
    public int getColumnIndex(String columnName);

    /**
     * Returns the name of the column at <code>metadataIndex</code>. Note: this
     * name needs to be unique; two columns in a meta data instance can not have
     * the same name.
     * 
     * @param metadataIndex
     *            the index of the column
     * @return the name of the column
     */
    public String getColumnName(int metadataIndex);

    /**
     * Returns the number of items that a vector in this instance contains.
     * 
     * @return the item count
     */
    public int getItemCount();

    /**
     * Returns the number of columns in this meta data instance.
     * 
     * @return the number of columns.
     */
    public int getMetadataCount();

    /**
     * Returns the value for the element at <code>itemIndex</code> and
     * <code>metadataIndex</code>.
     * 
     * @param itemIndex
     *            the item whose value is to be queried
     * @param metadataIndex
     *            the column name whose value is to be queried
     * @return the value Object at the specified element.
     */
    public Object getValue(int itemIndex, int metadataIndex);

    /**
     * Returns the value for the element at <code>metadataIndex</code> and
     * <code>columnName</code>.
     * 
     * @param itemIndex
     *            the item whose value is to be queried
     * @param columnName
     *            the column name whose value is to be queried
     * @return the value Object at the specified element.
     */
    public Object getValue(int itemIndex, String columnName);

    /**
     * Removes the column at the specified position in this meta data instance
     * (optional operation). Shifts any subsequent columns to the left
     * (subtracts one from their indices).
     * 
     * @param metadataIndex
     *            the column index to remove.
     * @throws UnsupportedOperationException
     *             if the <tt>remove</tt> method is not supported by this list.
     * @throws IndexOutOfBoundsException
     *             if the index is out of range (metadataIndex &lt; 0 || index
     *             &gt;= getColumnCount()).
     */
    public void remove(int metadataIndex);

    /**
     * Sets the value in the cell at <code>columnName</code> and
     * <code>itemIndex</code> to <code>aValue</code>.
     * 
     * @param aValue
     *            the new value
     * @param itemIndex
     *            the item whose value is to be changed
     * @param metadataIndex
     *            the metadata index whose value is to be changed
     * @see #getValue
     */

    public void setValue(int itemIndex, int metadataIndex, Object aValue);

    /**
     * Sets the value in the cell at <code>columnName</code> and
     * <code>itemIndex</code> to <code>aValue</code>.
     * 
     * @param aValue
     *            the new value
     * @param itemIndex
     *            the item whose value is to be changed
     * @param columnName
     *            the column name whose value is to be changed
     * @see #getValue
     */

    public void setValue(int itemIndex, String columnName, Object aValue);
}
