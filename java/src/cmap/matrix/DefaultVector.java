package cmap.matrix;

public class DefaultVector extends AbstractVector {
    private Object[] array;

    public DefaultVector(String columnName, Class<?> c, int size) {
	super(columnName, c);
	array = new Object[size];
    }

    public DefaultVector(String columnName, Class<?> c, Object[] array) {
	super(columnName, c);
	this.array = array;
    }

    /**
     * For use from R only
     * 
     * @param columnName
     * @param array
     */
    public DefaultVector(String columnName, Object obj) {
	super(columnName, String.class);
	array = createArray(obj);
	for (int i = 0, length = array.length; i < length; i++) {
	    if (array[i] != null) {
		setColumnClass(array[i].getClass());
		break;
	    }
	}

    }

    @Override
    public Object getValue(int index) {
	return array[index];
    }

    @Override
    public void setValue(int rowIndex, Object value) {
	array[rowIndex] = value;
    }

    @Override
    public int size() {
	return array.length;
    }

    private static Object[] createArray(Object obj) {
	Object[] array = null;
	if (obj instanceof double[]) {
	    double[] tmp = (double[]) obj;
	    array = new Object[tmp.length];
	    for (int i = 0, length = tmp.length; i < length; i++) {
		array[i] = tmp[i];
	    }
	} else if (obj instanceof String[]) {
	    String[] tmp = (String[]) obj;
	    array = new Object[tmp.length];
	    for (int i = 0, length = tmp.length; i < length; i++) {
		array[i] = tmp[i];
	    }
	} else if (obj instanceof int[]) {
	    int[] tmp = (int[]) obj;
	    array = new Object[tmp.length];
	    for (int i = 0, length = tmp.length; i < length; i++) {
		array[i] = tmp[i];
	    }
	} else if (obj instanceof boolean[]) {
	    boolean[] tmp = (boolean[]) obj;
	    array = new Object[tmp.length];
	    for (int i = 0, length = tmp.length; i < length; i++) {
		array[i] = tmp[i];
	    }
	} else {
	    System.err.println("Unknown type: " + obj.getClass());
	}
	return array;
    }

}