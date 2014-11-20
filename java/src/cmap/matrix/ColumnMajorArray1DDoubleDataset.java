package cmap.matrix;

public class ColumnMajorArray1DDoubleDataset extends FixedSizeAbstractDataset {
    private double[] array;

    public ColumnMajorArray1DDoubleDataset(String name, int rows, int columns,
	    double[] array) {
	super(name, rows, columns);
	this.array = array;
    }

    @Override
    public float getValue(int rowIndex, int columnIndex) {
	return (float) array[rowIndex + columnIndex * rows];
    }

    @Override
    public void setValue(int rowIndex, int columnIndex, float value) {
	array[rowIndex + columnIndex * rows] = value;
    }
}
