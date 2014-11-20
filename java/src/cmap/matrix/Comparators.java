package cmap.matrix;

import java.net.URL;
import java.util.Comparator;
import java.util.Date;
import java.util.LinkedHashMap;
import java.util.Map;

public class Comparators {
    private static final Map<Class<?>, Comparator> classToComparator;

    static {
	classToComparator = new LinkedHashMap<Class<?>, Comparator>();
	classToComparator.put(String.class, StringCaseInsensitiveComparator.INSTANCE);
	classToComparator.put(Number.class, ComparableComparator.INSTANCE);
	classToComparator.put(float.class, ComparableComparator.INSTANCE);
	classToComparator.put(double.class, ComparableComparator.INSTANCE);
	classToComparator.put(int.class, ComparableComparator.INSTANCE);
	classToComparator.put(Date.class, ComparableComparator.INSTANCE);
	classToComparator.put(URL.class, new ObjectComparator());

    }

    public static int compareToHandleNull(Object c1, Object c2) {
	if (c1 == null && c2 == null) {
	    return 0;
	} else if (c1 == null) {
	    return -1;
	} else if (c2 == null) {
	    return 1;
	}
	return 2;
    }

    public static Comparator getDefaultComparator(Class<?> c) {

	for (Class<?> test : classToComparator.keySet()) {
	    if (test.isAssignableFrom(c)) {
		return classToComparator.get(test);
	    }
	}
	return new ObjectComparator();
    }

    public static class ComparableComparator implements Comparator {
	public static final Comparator INSTANCE = new ComparableComparator();

	private ComparableComparator() {
	}

	public int compare(Object obj1, Object obj2) {
	    if ((obj1 == null) && (obj2 == null)) {
		return 0;
	    }
	    if (obj1 == null) {
		return -1;
	    }
	    if (obj2 == null) {
		return 1;
	    }
	    return ((Comparable) obj1).compareTo(obj2);
	}
    }

    public static class ObjectComparator implements Comparator<Object> {

	@Override
	public int compare(Object o1, Object o2) {
	    String s1 = String.valueOf(o1);
	    String s2 = String.valueOf(o2);
	    int value = compareToHandleNull(s1, s2);
	    return value == 2 ? String.CASE_INSENSITIVE_ORDER.compare(s1, s2) : value;
	}
    }

    public static class StringCaseInsensitiveComparator implements Comparator<String> {
	public static final Comparator<String> INSTANCE = new StringCaseInsensitiveComparator();

	@Override
	public int compare(String o1, String o2) {
	    int value = compareToHandleNull(o1, o2);
	    return value == 2 ? String.CASE_INSENSITIVE_ORDER.compare(o1, o2) : value;
	}
    }

}
