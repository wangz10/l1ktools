package cmap.matrix;

/**
 * 
 * A tranposed view of a {@link Dataset}.
 *
 */

public class TransposedDatasetView extends DatasetAdapter {

    public TransposedDatasetView(Dataset dataset) {
	super(dataset);
    }

    @Override
    public int getColumnCount() {
	return dataset.getRowCount();
    }

    @Override
    public MetadataModel getColumnMetadata() {
	return dataset.getRowMetadata();
    }

    @Override
    public Object getObjectValue(int row, int column, int seriesIndex) {
	return dataset.getObjectValue(column, row, seriesIndex);
    }

    @Override
    public int getRowCount() {
	return dataset.getColumnCount();
    }

    @Override
    public MetadataModel getRowMetadata() {
	return dataset.getColumnMetadata();
    }

    @Override
    public float getValue(int row, int column) {
	return dataset.getValue(column, row);
    }

    @Override
    public float getValue(int row, int column, int seriesIndex) {
	return dataset.getValue(column, row, seriesIndex);
    }

    @Override
    public void setValue(int rowIndex, int columnIndex, float value) {
	dataset.setValue(columnIndex, rowIndex, value);
    }
}
