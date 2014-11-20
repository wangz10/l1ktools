package cmap.matrix;

import java.util.Collection;
import java.util.Comparator;

/**
 * 
 * A collection of values of the same type.
 * 
 */
public interface Vector {
    /**
     * Returns the most specific superclass for all the values in this vector.
     * 
     * @return the common ancestor class of the object values in this instance
     */
    public Class<?> getColumnClass();

    /**
     * Returns the comparator used to compare values in this vector.
     * 
     * @return the comparator
     */
    public Comparator<Object> getComparator();

    /**
     * Returns the name of this vector.
     * 
     * @return the name
     */
    public String getName();

    /**
     * Gets the value for the specified key.
     * 
     * @param key
     *            the key
     * @return the value
     */
    public Object getProperty(Object key);

    /**
     * Returns a collection of all property keys in this vector.
     * 
     * @return the collection of keys
     */
    public Collection<Object> getPropertyKeys();

    /**
     * Returns the value at the specified index.
     * 
     * @param index
     *            the index
     * @return the value.
     */
    public Object getValue(int index);

    /**
     * Removes the specified property.
     * 
     * @param property
     *            The property key
     */
    public void removeProperty(Object key);

    /**
     * Sets the most specific superclass for all the values in this vector.
     * 
     * @param columnClass
     *            the common ancestor class of the object values in this
     *            instance
     */
    public void setColumnClass(Class<?> columnClass);

    /**
     * Sets the comparator used to compare values in this vector.
     * 
     * @param comparator
     *            the comparator
     */
    public void setComparator(Comparator<Object> comparator);

    /**
     * Sets the property for the specified key.
     * 
     * @param key
     *            the key
     * @param value
     *            the value
     */
    public void setProperty(Object key, Object value);

    /**
     * Sets the value at the specified index.
     * 
     * @param index
     *            the index
     * @param value
     *            the value
     */
    public void setValue(int index, Object value);

    /**
     * Returns the number of elements in this vector.
     * 
     * @return the size.
     */
    public int size();
}
