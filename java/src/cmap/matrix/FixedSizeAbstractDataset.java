/**
 * Copyright 2009-2010 by The Broad Institute.
 * 
 
 */
package cmap.matrix;

/**
 * 
 * Fixed size abstract dataset.
 * 
 */
public abstract class FixedSizeAbstractDataset extends AbstractDataset {

    protected MetadataModel columnMetadata;

    protected int columns;

    protected MetadataModel rowMetadata;

    protected int rows;

    /**
     * Creates a new dataset
     * 
     * @param name
     *            the dataset name
     * @param rows
     *            the number of rows
     * @param columns
     *            the number of columns
     */
    protected FixedSizeAbstractDataset(String name, int rows, int columns) {
	this.name = name;
	this.rows = rows;
	this.columns = columns;
	rowMetadata = new DefaultMetadataModel(rows);
	columnMetadata = new DefaultMetadataModel(columns);
    }

    protected FixedSizeAbstractDataset(String name, MetadataModel rowMetadata, MetadataModel columnMetadata) {
	this.name = name;
	this.rows = rowMetadata.getItemCount();
	this.columns = columnMetadata.getItemCount();
	this.rowMetadata = rowMetadata;
	this.columnMetadata = columnMetadata;
    }

    /**
     * Inserts the specified series.
     * 
     * @param name
     *            series name to be inserted
     * @param values
     *            array of values
     * @param c
     *            class of values
     * @return the series index
     */
    public int addSeries(String name, Object[][] values, Class<?> c) {
	seriesNames.add(name);
	seriesClasses.add(c);
	matrices.add(new ObjectMatrix(values));
	return seriesNames.size();
    }

    @Override
    public int getColumnCount() {
	return columns;
    }

    @Override
    public MetadataModel getColumnMetadata() {
	return columnMetadata;
    }

    @Override
    public int getRowCount() {
	return rows;
    }

    @Override
    public MetadataModel getRowMetadata() {
	return rowMetadata;
    }

    public void setColumnMetadata(MetadataModel metadata) {
	if (metadata.getItemCount() != getColumnCount()) {
	    throw new IllegalArgumentException("Wrong number of items.");
	}
	this.columnMetadata = metadata;
    }

    public void setRowMetadata(MetadataModel metadata) {
	if (metadata.getItemCount() != getRowCount()) {
	    throw new IllegalArgumentException("Wrong number of items.");
	}
	this.rowMetadata = metadata;
    }

    protected Matrix createMatrix() {
	return new ObjectMatrix(getRowCount(), getColumnCount());
    }

    static class ObjectMatrix implements Matrix {
	private Object[][] matrix;

	public ObjectMatrix(Object[][] values) {
	    this.matrix = values;
	}

	private ObjectMatrix(int rows, int columns) {
	    matrix = new Object[rows][columns];
	}

	@Override
	public Object get(int row, int column) {
	    return matrix[row][column];
	}

	@Override
	public void set(int row, int column, Object value) {
	    matrix[row][column] = value;
	}
    }

}
