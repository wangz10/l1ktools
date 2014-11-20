package cmap.matrix;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public abstract class AbstractMetadataModel implements MetadataModel {
    private Map<String, Integer> columnNameToIndex = new HashMap<String, Integer>();

    private List<Vector> vectors = new ArrayList<Vector>();

    protected AbstractMetadataModel() {
    }

    @Override
    public Vector add(String name, Class<?> c) {
	if (name == null) {
	    throw new NullPointerException("Name is null");
	}
	if (c == null) {
	    throw new NullPointerException("Class is null");
	}
	Vector vector = createVector(name, c);
	add(vector);
	return vector;
    }

    @Override
    public Vector get(int columnIndex) {
	return vectors.get(columnIndex);
    }

    @Override
    public Vector get(String name) {
	int index = getColumnIndex(name);
	return index == -1 ? null : vectors.get(index);
    }

    public Class<?> getColumnClass(int columnIndex) {
	return vectors.get(columnIndex).getColumnClass();
    }

    @Override
    public int getColumnIndex(String name) {
	Integer index = columnNameToIndex.get(name);
	return index != null ? index : -1;
    }

    @Override
    public String getColumnName(int columnIndex) {
	return vectors.get(columnIndex).getName();
    }

    @Override
    public int getMetadataCount() {
	return vectors.size();
    }

    @Override
    public Object getValue(int rowIndex, int columnIndex) {
	Vector metadataColumn = vectors.get(columnIndex);
	return metadataColumn.getValue(rowIndex);
    }

    @Override
    public Object getValue(int rowIndex, String column) {
	Integer columnIndex = columnNameToIndex.get(column);
	if (columnIndex == null) {
	    throw new IllegalArgumentException("Column " + column + " not found.");
	}
	return getValue(rowIndex, columnIndex);

    }

    public Object getValueAt(int rowIndex, int columnIndex) {
	Vector metadataColumn = vectors.get(columnIndex);
	return metadataColumn.getValue(rowIndex);
    }

    @Override
    public void remove(int index) {
	vectors.remove(index);
	rebuildMap();
    }

    @Override
    public void setValue(int rowIndex, int columnIndex, Object value) {
	Vector metadataColumn = vectors.get(columnIndex);
	metadataColumn.setValue(rowIndex, value);
    }

    @Override
    public void setValue(int rowIndex, String name, Object value) {
	if (name == null) {
	    throw new NullPointerException("Name is null");
	}
	Integer columnIndex = columnNameToIndex.get(name);
	if (columnIndex == null) {
	    StringBuilder buf = new StringBuilder();
	    for (int j = 0, cols = getMetadataCount(); j < cols; j++) {
		if (j > 0) {
		    buf.append(", ");
		}
		buf.append(vectors.get(j).getName());
	    }
	    throw new NullPointerException("No column named '" + name + "' found. Column names are: " + buf.toString());
	}
	setValue(rowIndex, columnIndex, value);
    }

    @Override
    public String toString() {
	StringBuilder buf = new StringBuilder();
	for (int i = 0, nitems = getItemCount(); i < nitems; i++) {
	    for (int j = 0, mcount = getMetadataCount(); j < mcount; j++) {
		if (j > 0) {
		    buf.append("\t");
		}
		Object value = get(j).getValue(i);
		buf.append(value);
	    }
	    buf.append("\n");
	}
	return buf.toString();
    }

    protected void add(Vector vector) {
	if (columnNameToIndex.containsKey(vector.getName())) {
	    int existingIndex = columnNameToIndex.get(vector.getName());
	    vectors.remove(existingIndex);
	    columnNameToIndex.remove(vector.getName());
	    rebuildMap();
	}
	vectors.add(vector);
	columnNameToIndex.put(vector.getName(), vectors.size() - 1);
    }

    protected abstract Vector createVector(String name, Class<?> c);

    private void rebuildMap() {
	columnNameToIndex.clear();
	for (int i = 0, size = vectors.size(); i < size; i++) {
	    columnNameToIndex.put(vectors.get(i).getName(), i);
	}
    }

}
