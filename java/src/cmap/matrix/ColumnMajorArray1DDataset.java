package cmap.matrix;

public class ColumnMajorArray1DDataset extends FixedSizeAbstractDataset {
    private float[] array;

    public ColumnMajorArray1DDataset(String name, int rows, int columns, float[] array) {
	super(name, rows, columns);
	this.array = array;
    }

    public float[] getArray() {
        return array;
    }

    @Override
    public float getValue(int rowIndex, int columnIndex) {
	try {
	    return array[rowIndex + columnIndex * rows];
	} catch (ArrayIndexOutOfBoundsException x) {
	    throw new RuntimeException("ArrayIndexOutOfBoundsException at " + rowIndex + ", " + columnIndex);
	}
    }

    @Override
    public void setValue(int rowIndex, int columnIndex, float value) {
	array[rowIndex + columnIndex * rows] = value;
    }
}
