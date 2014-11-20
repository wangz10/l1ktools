package cmap.matrix;

public class BigColumnMajorArray1DDataset extends FixedSizeAbstractDataset {

    private float[][] arrays;
    private int ncolsPerChunk;

    public BigColumnMajorArray1DDataset(String name, int rows, int columns,
	    float[][] arrays, int ncolsPerChunk) {
	super(name, rows, columns);
	this.arrays = arrays;
	// totalColumnCount = new long[ncols.length];
	// long total = 0;
	// for (int i = 0, l = totalColumnCount.length; i < l; i++) {
	// total += ncols[i];
	// totalColumnCount[i] = total;
	// }
	this.ncolsPerChunk = ncolsPerChunk;
    }

    @Override
    public float getValue(int rowIndex, int columnIndex) {
	int arrayIndex = columnIndex / ncolsPerChunk;
	float[] array = arrays[arrayIndex];
	columnIndex = columnIndex % ncolsPerChunk;
	return array[rowIndex + columnIndex * rows];
    }

    @Override
    public void setValue(int rowIndex, int columnIndex, float value) {
	int arrayIndex = columnIndex / ncolsPerChunk;
	float[] array = arrays[arrayIndex];
	columnIndex = columnIndex % ncolsPerChunk;
	array[rowIndex + columnIndex * rows] = value;
    }
}
