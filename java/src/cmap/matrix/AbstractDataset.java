package cmap.matrix;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * 
 * Abstract dataset.
 * 
 */
public abstract class AbstractDataset implements Dataset {

    protected List<Matrix> matrices = new ArrayList<Matrix>();
    protected String name;
    protected Map<Object, Object> properties = new HashMap<Object, Object>();
    protected List<Class<?>> seriesClasses = new ArrayList<Class<?>>();
    protected List<String> seriesNames = new ArrayList<String>();

    public AbstractDataset() {
    }

    public int addSeries(String name, Class<?> c) {
	seriesNames.add(name);
	seriesClasses.add(c);
	matrices.add(createMatrix());
	return seriesNames.size();
    }

    @Override
    public final String getName() {
	return name;
    }

    @Override
    public Object getObjectValue(int rowIndex, int columnIndex, int seriesIndex) {
	if (seriesIndex == 0) {
	    return getValue(rowIndex, columnIndex);
	}
	Matrix matrix = matrices.get(seriesIndex - 1);
	return matrix.get(rowIndex, columnIndex);
    }

    @Override
    public Object getProperty(Object key) {
	return properties.get(key);
    }

    @Override
    public Collection<Object> getPropertyKeys() {
	return properties.keySet();
    }

    @Override
    public Class<?> getSeriesClass(int series) {
	if (series == 0) {
	    return Float.class;
	}
	return seriesClasses.get(series - 1);
    }

    @Override
    public int getSeriesCount() {
	return matrices.size() + 1;
    }

    @Override
    public String getSeriesName(int series) {
	if (series == 0) {
	    return getName();
	}
	return seriesNames.get(series - 1);
    }

    @Override
    public float getValue(int rowIndex, int columnIndex, int seriesIndex) {
	if (seriesIndex == 0) {
	    return getValue(rowIndex, columnIndex);
	}
	Matrix matrix = matrices.get(seriesIndex - 1);
	Number number = (Number) matrix.get(rowIndex, columnIndex);
	return number != null ? number.floatValue() : Float.NaN;
    }

    public void removeSeries(int index) {
	if (index > 0) {
	    seriesNames.remove(index - 1);
	    seriesClasses.remove(index - 1);
	    matrices.remove(index - 1);
	}
    }

    @Override
    public void setName(String name) {
	this.name = name;
    }

    /**
     * Sets the value at the given row and column.
     * 
     * @param rowIndex
     *            the row index.
     * @param columnIndex
     *            the column index.
     * @param seriesIndex
     *            the series index.
     * @throws IndexOutOfBoundsException
     *             if <tt>rowIndex</tt> or <tt>columnIndex</tt> are out of
     *             range.
     */
    public void setObjectValue(int rowIndex, int columnIndex, Object value, int seriesIndex) {
	if (seriesIndex == 0) {
	    setValue(rowIndex, columnIndex, ((Number) value).floatValue());
	} else {
	    Matrix m = matrices.get(seriesIndex - 1);
	    m.set(rowIndex, columnIndex, value);
	}
    }

    @Override
    public void setProperty(Object key, Object value) {
	if (value == null) {
	    properties.remove(key);
	} else {
	    properties.put(key, value);
	}
    }

    protected abstract Matrix createMatrix();

    protected static interface Matrix {

	public Object get(int row, int column);

	public void set(int row, int column, Object value);
    }

}
