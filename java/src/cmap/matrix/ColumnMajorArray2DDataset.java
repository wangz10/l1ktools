package cmap.matrix;

public class ColumnMajorArray2DDataset extends FixedSizeAbstractDataset {
    private float[][] array;

    public ColumnMajorArray2DDataset(String name, int rows, int columns) {
	super(name, rows, columns);
	this.array = new float[columns][rows];
    }

    public ColumnMajorArray2DDataset(String name, int rows, int columns, float[][] array) {
	super(name, rows, columns);
	this.array = array;
    }

    public ColumnMajorArray2DDataset(String name, MetadataModel rowMetadata, MetadataModel columnMetadata,
	    float[][] array) {
	super(name, rowMetadata, columnMetadata);
	this.array = array;
    }

    @Override
    public float getValue(int rowIndex, int columnIndex) {
	return array[columnIndex][rowIndex];
    }

    @Override
    public void setValue(int rowIndex, int columnIndex, float value) {
	array[columnIndex][rowIndex] = value;
    }
}
