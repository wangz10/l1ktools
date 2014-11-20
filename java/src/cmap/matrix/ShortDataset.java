package cmap.matrix;

public class ShortDataset extends FixedSizeAbstractDataset {
    private short[] array;

    public ShortDataset(String name, int rows, int columns, short[] array) {
	super(name, rows, columns);
	this.array = array;
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
	array[rowIndex + columnIndex * rows] = (short) value;
    }
}
