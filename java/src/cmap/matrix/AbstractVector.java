package cmap.matrix;

import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;

public abstract class AbstractVector implements Vector {
    private Class<?> columnClass;
    private Comparator<Object> comparator;
    private String name;
    private Map<Object, Object> properties = new HashMap<Object, Object>();

    public AbstractVector(String columnName, Class<?> c) {
	this.name = columnName;
	this.columnClass = c;
	this.comparator = Comparators.getDefaultComparator(columnClass);
    }

    @Override
    public Class<?> getColumnClass() {
	return columnClass;
    }

    @Override
    public Comparator<Object> getComparator() {
	return comparator;
    }

    @Override
    public String getName() {
	return name;
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
    public void removeProperty(Object key) {
	properties.remove(key);
    }

    @Override
    public void setColumnClass(Class<?> c) {
	this.columnClass = c;
	this.comparator = null;
	this.comparator = Comparators.getDefaultComparator(columnClass);
    }

    @Override
    public void setComparator(Comparator<Object> comparator) {
	if (comparator == null) {
	    comparator = Comparators.getDefaultComparator(columnClass);
	}
	this.comparator = comparator;

    }

    @Override
    public void setProperty(Object key, Object value) {
	properties.put(key, value);

    }

    @Override
    public String toString() {
	return name;
    }
}
