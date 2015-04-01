package cmap.io;

import static ncsa.hdf.hdf5lib.H5.H5Dclose;
import static ncsa.hdf.hdf5lib.H5.H5Dget_space;
import static ncsa.hdf.hdf5lib.H5.H5Dget_type;
import static ncsa.hdf.hdf5lib.H5.H5Dopen;
import static ncsa.hdf.hdf5lib.H5.H5Dread;
import static ncsa.hdf.hdf5lib.H5.H5Lexists;
import static ncsa.hdf.hdf5lib.H5.H5Sclose;
import static ncsa.hdf.hdf5lib.H5.H5Screate_simple;
import static ncsa.hdf.hdf5lib.H5.H5Sget_simple_extent_dims;
import static ncsa.hdf.hdf5lib.H5.H5Sget_simple_extent_ndims;
import static ncsa.hdf.hdf5lib.H5.H5Sselect_hyperslab;
import static ncsa.hdf.hdf5lib.H5.H5Tclose;
import static ncsa.hdf.hdf5lib.H5.H5Tcopy;
import static ncsa.hdf.hdf5lib.H5.H5Tequal;
import static ncsa.hdf.hdf5lib.H5.H5Tget_native_type;
import static ncsa.hdf.hdf5lib.H5.H5Tget_size;
import static ncsa.hdf.hdf5lib.H5.H5Tset_size;
import static ncsa.hdf.hdf5lib.HDF5Constants.H5P_DEFAULT;
import static ncsa.hdf.hdf5lib.HDF5Constants.H5S_SELECT_NOTB;
import static ncsa.hdf.hdf5lib.HDF5Constants.H5S_SELECT_OR;
import static ncsa.hdf.hdf5lib.HDF5Constants.H5S_SELECT_SET;
import static ncsa.hdf.hdf5lib.HDF5Constants.H5S_SELECT_XOR;
import static ncsa.hdf.hdf5lib.HDF5Constants.H5T_FORTRAN_S1;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import ncsa.hdf.hdf5lib.H5;
import ncsa.hdf.hdf5lib.HDF5Constants;
import ncsa.hdf.object.Dataset;
import ncsa.hdf.object.FileFormat;
import ncsa.hdf.object.Group;
import ncsa.hdf.object.h5.H5Datatype;
import cmap.matrix.BigColumnMajorArray1DDataset;
import cmap.matrix.ColumnMajorArray1DDataset;
import cmap.matrix.ColumnMajorArray1DDoubleDataset;
import cmap.matrix.ColumnMajorArray2DDataset;
import cmap.matrix.DatasetAdapter;
import cmap.matrix.DefaultMetadataModel;
import cmap.matrix.MetadataModel;
import cmap.matrix.ShortDataset;
import cmap.matrix.SlicedDatasetView;
import cmap.matrix.Vector;

/**
 * 
 * Read files in gctx format.
 *
 */
public class GctxReader {

    public static final String COL_METADATA_PATH = "/0/META/COL";
    public static final String DATA_MATRIX_PATH = "/0/DATA/0/matrix";
    public static final String ROW_METADATA_PATH = "/0/META/ROW";

    private FileFormat h5File;

    private String path;

    /**
     * Creates a new reader.
     * 
     * @param path
     *            The input file path.
     */
    public GctxReader(String path) {
	try {
	    if (!new File(path).exists()) {
		throw new IOException(path + " not found.");
	    }
	    this.path = path;
	    FileFormat fileFormat = FileFormat.getFileFormat(FileFormat.FILE_TYPE_HDF5);
	    if (fileFormat == null) {
		throw new RuntimeException("Cannot find HDF5 FileFormat.");
	    }
	    h5File = fileFormat.createInstance(path, FileFormat.READ);
	    h5File.open();
	} catch (Exception x) {
	    throw new RuntimeException(x);
	}

    }

    public void close() {
	try {
	    if (h5File != null) {
		h5File.close();
		h5File = null;
	    }
	    System.gc();
	} catch (Exception e) {
	    e.printStackTrace();
	}
    }

    public DatasetInfo getDatasetInfo() {
	String matrixPath = GctxReader.DATA_MATRIX_PATH;

	try {
	    Dataset tmpDataset = (Dataset) h5File.get(matrixPath);
	    if (tmpDataset == null) {
		throw new NullPointerException("No data found at " + matrixPath);
	    }
	    tmpDataset.init();

	    long[] dims = tmpDataset.getDims();
	    long[] chunk = tmpDataset.getChunkSize();
	    DatasetInfo info = new DatasetInfo();
	    info.nrows = (int) dims[1];
	    info.ncols = (int) dims[0];
	    info.rowChunkSize = chunk == null ? -1 : (int) chunk[1];
	    info.columnChunkSize = chunk == null ? -1 : (int) chunk[0];

	    return info;
	} catch (Exception x) {
	    throw new RuntimeException(x);

	}

    }

    public Set<String> getDatasets() {
	try {
	    Set<String> names = new HashSet<String>();
	    Group group = (Group) h5File.get("/DATA");
	    if (group != null) {
		for (Object obj : group.getMemberList()) {
		    if (obj instanceof Group) {
			Group datasetGroup = (Group) obj;
			names.add(datasetGroup.getName());
		    }
		}
	    }
	    return names;
	} catch (Exception x) {
	    throw new RuntimeException(x);

	}
    }

    public String getFile() {
	return path;
    }

    public List<String> getRowsOrColumnNames(boolean isRows) {
	String path = (!isRows ? GctxReader.COL_METADATA_PATH : GctxReader.ROW_METADATA_PATH);
	if (!pathExists(path)) {
	    path = isRows ? "/GIN" : "/SIN";
	}
	List<String> names = new ArrayList<String>();

	try {
	    int count = H5.H5Gn_members(h5File.getFID(), path);
	    if (count > 0) {
		String[] oname = new String[count];
		int[] otype = new int[count];
		int[] ltype = new int[count];
		long[] orefs = new long[count];
		H5.H5Gget_obj_info_all(h5File.getFID(), path, oname, otype, ltype, orefs, HDF5Constants.H5_INDEX_NAME);

		for (int indx = 0; indx < otype.length; indx++) {
		    String name = oname[indx];
		    names.add(name);
		}
	    }
	} catch (Exception e) {
	    e.printStackTrace();
	}

	return names;
    }

    public MetadataModel getRowsOrColumns(boolean isRows) {
	return getRowsOrColumns(isRows, null, null);
    }

    public MetadataModel getRowsOrColumns(boolean isRows, int[] indices, String[] fields) {
	String path = (!isRows ? GctxReader.COL_METADATA_PATH : GctxReader.ROW_METADATA_PATH);
	if (!pathExists(path)) {
	    path = isRows ? "/GIN" : "/SIN";
	}

	return getRowsOrColumns(path, indices, fields, isRows);
    }

    public boolean pathExists(String path) {
	try {
	    return h5File.get(path) != null;
	} catch (Exception x) {
	    throw new RuntimeException(x);
	}
    }

    public cmap.matrix.Dataset read() {
	return read(null, null, true, true, true, true, -1);
    }

    public cmap.matrix.Dataset read(int[] rowIndices, int[] columnIndices, boolean readRowMetadata,
	    boolean readColumnMetadata, boolean rowIndicesAreAscending, boolean columnIndicesAreAscending) {
	return read(rowIndices, columnIndices, readRowMetadata, readColumnMetadata, rowIndicesAreAscending,
		columnIndicesAreAscending, -1);
    }

    /**
     * Read a dataset.
     * 
     * @param rowIndices
     *            The row indices or <tt>null</tt>.
     * @param columnIndices
     *            The column indices or <tt>null</tt>
     * @param readRowMetadata
     *            Whether to read the row metadata.
     * @param readColumnMetadata
     *            Whether to read the column metadata.
     * @param indicesAreAscending
     *            Whether the row indices or column indices (if provided) are in
     *            ascending order. If this flag is <tt>false</tt>, then the
     *            indices in the returned dataset will be in the same order as
     *            the specified order, otherwise the indices are reshuffled so
     *            that they are in ascending order (the returned indices may not
     *            be in the same order).
     * @param memoryDataType
     *            the memory data type to read the data into or -1 to read it as
     *            the data type that is it stored on disk.
     * @return the dataset.
     */
    public cmap.matrix.Dataset read(int[] rowIndices, int[] columnIndices, boolean readRowMetadata,
	    boolean readColumnMetadata, boolean rowIndicesAreAscending, boolean columnIndicesAreAscending,
	    int memoryDataType) {
	try {
	    String matrixPath = GctxReader.DATA_MATRIX_PATH;
	    int[] rowIndexTableIndex = null;
	    if (!rowIndicesAreAscending && rowIndices != null) {
		rowIndexTableIndex = index(rowIndices);
	    }

	    int[] columnIndexTableIndex = null;
	    if (!columnIndicesAreAscending && columnIndices != null) {
		columnIndexTableIndex = index(columnIndices);
	    }
	    if (!pathExists(matrixPath)) {
		matrixPath = "/DATA/GCT/matrix";
	    }
	    if (!pathExists(matrixPath)) {
		matrixPath = "/0/DATA/matrix0";
	    }
	    if (!pathExists(matrixPath)) {
		throw new IOException("Path " + GctxReader.DATA_MATRIX_PATH + " not found.");
	    }
	    Dataset tmpDataset = (Dataset) h5File.get(matrixPath);
	    tmpDataset.init();

	    int rowCount = rowIndices != null ? rowIndices.length : (int) tmpDataset.getDims()[1];
	    int columnCount = columnIndices != null ? columnIndices.length : (int) tmpDataset.getDims()[0];
	    if (columnCount == 0) {
		throw new IllegalArgumentException("Zero columns specified.");
	    }

	    // if (H5Aexists_by_name(h5File.getFID(), matrixPath, "row_major",
	    // H5P_DEFAULT)) {
	    // int attribute_id = H5Aopen_by_name(h5File.getFID(), matrixPath,
	    // "row_major", H5P_DEFAULT, H5P_DEFAULT);
	    // int[] tmp = new int[1];
	    // H5Aread(attribute_id, HDF5Constants.H5T_NATIVE_INT, tmp);
	    // rowMajorOrder = tmp[0] == 1;
	    // H5Aclose(attribute_id);
	    // }
	    cmap.matrix.Dataset dataset = null;
	    Object matrix = null;
	    // matrix is stored in column major order
	    List<Hyperslab> slabs = new ArrayList<Hyperslab>();
	    if (columnIndices != null) {
		slabs.add(new Hyperslab(columnIndices, 1));
	    }
	    if (rowIndices != null) {
		slabs.add(new Hyperslab(rowIndices, 2));
	    }
	    try {
		matrix = read(h5File.getFID(), matrixPath, slabs.toArray(new Hyperslab[0]), memoryDataType);
	    } catch (ArrayTooBigException x) {
		if (columnIndices != null) {
		    throw new IllegalArgumentException("User chunking is only supported when columnIndices is null.");
		}
		// only support whole dataset for now
		int ncolsPerChunk = Integer.MAX_VALUE / rowCount;
		int nchunks = (int) Math.ceil(((float) columnCount) / ncolsPerChunk);
		System.out.println("Array is too big-chunking into " + nchunks + " arrays.");
		float[][] array2d = new float[nchunks][];

		Hyperslab[] _slabs = null;
		if (rowIndices != null) {
		    _slabs = new Hyperslab[2];
		    _slabs[1] = slabs.get(0);
		} else {
		    _slabs = new Hyperslab[1];
		}
		for (int j = 0, columnIndex = 0; j < nchunks; j++, columnIndex += ncolsPerChunk) {
		    int arraySize;
		    if (j == nchunks - 1) {
			arraySize = columnCount - columnIndex;
		    } else {
			arraySize = ncolsPerChunk;
		    }
		    int[] tmp = new int[arraySize];
		    for (int i = 0, k = columnIndex; i < arraySize; i++, k++) {
			tmp[i] = k;
		    }

		    _slabs[0] = new Hyperslab(tmp, 1);
		    array2d[j] = (float[]) read(h5File.getFID(), matrixPath, _slabs, -1);
		}
		dataset = new BigColumnMajorArray1DDataset(new File(path).getName(), rowCount, columnCount, array2d,
			ncolsPerChunk);

	    }

	    if (dataset == null) {
		if (matrix instanceof float[]) {
		    dataset = new ColumnMajorArray1DDataset(new File(path).getName(), rowCount, columnCount,
			    (float[]) matrix);
		} else if (matrix instanceof short[]) {
		    dataset = new ShortDataset(new File(path).getName(), rowCount, columnCount, (short[]) matrix);

		} else if (matrix instanceof float[][]) {
		    dataset = new ColumnMajorArray2DDataset(new File(path).getName(), rowCount, columnCount,
			    (float[][]) matrix);
		} else if (matrix instanceof double[]) {
		    dataset = new ColumnMajorArray1DDoubleDataset(new File(path).getName(), rowCount, columnCount,
			    (double[]) matrix);
		} else {
		    throw new IllegalArgumentException(matrix != null ? "Unknown matrix data type " + matrix.getClass()
			    + "." : "Matrix is null.");
		}
	    }

	    MetadataModel columnMetadata = !readColumnMetadata ? null : getRowsOrColumns(false, columnIndices, null);
	    if (columnMetadata == null) {
		columnMetadata = new DefaultMetadataModel(columnCount);
	    }
	    MetadataModel rowMetadata = !readRowMetadata ? null : getRowsOrColumns(true, rowIndices, null);
	    if (rowMetadata == null) {
		rowMetadata = new DefaultMetadataModel(rowCount);
	    }
	    final MetadataModel _columnMetadata = columnMetadata;
	    final MetadataModel _rowMetadata = rowMetadata;
	    dataset = new DatasetAdapter(dataset) {
		@Override
		public MetadataModel getColumnMetadata() {
		    return _columnMetadata;
		}

		@Override
		public MetadataModel getRowMetadata() {
		    return _rowMetadata;
		}
	    };
	    if (rowIndexTableIndex != null || columnIndexTableIndex != null) {
		dataset = new SlicedDatasetView(dataset, rowIndexTableIndex, columnIndexTableIndex);
	    }
	    return dataset;
	} catch (Exception x) {
	    throw new RuntimeException(x);
	}
    }

    public byte[] readBytes(String path) throws Exception {
	Dataset dataset = (Dataset) h5File.get(path);
	dataset.init();
	return dataset.readBytes();
    }

    private MetadataModel addToMetadata(MetadataModel metadata, String path, int[] indices, boolean rows)
	    throws Exception {
	Object obj = read(h5File.getFID(), path,
		indices == null ? null : new Hyperslab[] { new Hyperslab(indices, 1) }, -1);
	Vector vector = null;
	String name = path.substring(path.lastIndexOf('/') + 1);
	if (obj instanceof String[]) {
	    String[] array = (String[]) obj;
	    if (metadata == null) {
		metadata = new DefaultMetadataModel(array.length);
	    }
	    vector = metadata.add(name, String.class);
	    for (int i = 0, length = metadata.getItemCount(); i < length; i++) {
		String val = array[i];
		if (val != null && val.equals("")) {
		    val = null;
		}
		vector.setValue(i, val);
	    }
	} else if (obj instanceof float[]) {
	    float[] array = (float[]) obj;
	    if (metadata == null) {
		metadata = new DefaultMetadataModel(array.length);
	    }
	    vector = metadata.add(name, Float.class);
	    for (int i = 0, length = metadata.getItemCount(); i < length; i++) {

		vector.setValue(i, array[i]);
	    }
	} else if (obj instanceof double[]) {
	    double[] array = (double[]) obj;
	    if (metadata == null) {
		metadata = new DefaultMetadataModel(array.length);
	    }
	    vector = metadata.add(name, Double.class);
	    for (int i = 0, length = metadata.getItemCount(); i < length; i++) {
		vector.setValue(i, array[i]);
	    }
	} else if (obj instanceof int[]) {
	    int[] array = (int[]) obj;
	    if (metadata == null) {
		metadata = new DefaultMetadataModel(array.length);
	    }
	    vector = metadata.add(name, Integer.class);
	    for (int i = 0, length = metadata.getItemCount(); i < length; i++) {
		vector.setValue(i, array[i]);
	    }
	}
	return metadata;
    }

    private MetadataModel getRowsOrColumns(String path, int[] indices, String[] fields, boolean rows) {
	try {
	    if (H5Lexists(h5File.getFID(), path, H5P_DEFAULT)) {
		MetadataModel model = null;
		if (fields == null) {
		    int count = H5.H5Gn_members(h5File.getFID(), path);
		    if (count > 0) {
			String[] oname = new String[count];
			int[] otype = new int[count];
			int[] ltype = new int[count];
			long[] orefs = new long[count];
			H5.H5Gget_obj_info_all(h5File.getFID(), path, oname, otype, ltype, orefs,
				HDF5Constants.H5_INDEX_NAME);

			for (int indx = 0; indx < otype.length; indx++) {
			    String name = oname[indx];
			    model = addToMetadata(model, path + "/" + name, indices, rows);
			}
		    }
		} else {
		    for (String name : fields) {
			String metadataPath = path + "/" + name;
			if (H5Lexists(h5File.getFID(), metadataPath, H5P_DEFAULT)) {
			    model = addToMetadata(model, metadataPath, indices, rows);
			} else {
			    metadataPath = path + "/" + name.toLowerCase();
			    if (H5Lexists(h5File.getFID(), metadataPath, H5P_DEFAULT)) {
				model = addToMetadata(model, metadataPath, indices, rows);
			    }
			}

		    }
		}

		return model;
	    } else {
		System.out.println("empty " + (rows ? "row" : "column") + " metadata for " + this.path);
		return null;
	    }
	} catch (Exception x) {
	    x.printStackTrace();
	    throw new RuntimeException(x);

	}
    }

    public static Object read(int file_id, String path, Hyperslab[] slabs, int memTypeId) throws Exception {
	boolean isString = false;
	int dataset_id = -1;
	int dataspace_id = -1;
	int tid = -1;
	int mspace = -1;
	try {
	    dataset_id = H5Dopen(file_id, path, H5P_DEFAULT);
	    dataspace_id = H5Dget_space(dataset_id);
	    int rank = H5Sget_simple_extent_ndims(dataspace_id);
	    if (rank == 0) {
		// a scalar data point
		rank = 1;
	    }
	    long[] dims = new long[rank];
	    H5Sget_simple_extent_dims(dataspace_id, dims, null);

	    long[] memSpaceDimensions = dims.clone();
	    mspace = H5P_DEFAULT;
	    if (slabs != null && slabs.length > 0) {
		for (Hyperslab slab : slabs) {
		    memSpaceDimensions[slab.rank - 1] = slab.indices.length;
		}

		mspace = H5Screate_simple(rank, memSpaceDimensions, null);
		if (slabs.length == 1) {
		    Hyperslab slab = slabs[0];
		    long[] count = dims.clone();
		    long[] start = new long[rank];
		    count[slab.rank - 1] = 1;
		    int[] indices = slab.indices;
		    for (int j = 0, length = indices.length; j < length; j++) {
			start[slab.rank - 1] = slab.indices[j];
			H5Sselect_hyperslab(dataspace_id, j == 0 ? H5S_SELECT_SET : H5S_SELECT_OR, start, null, count,
				null);
		    }
		} else { // only handle slab length of 2

		    Hyperslab slab0 = slabs[0];
		    long[] count = dims.clone();
		    long[] start = new long[rank];
		    count[slab0.rank - 1] = 1;

		    for (int j = 0, length = slab0.indices.length; j < length; j++) {
			start[slab0.rank - 1] = slab0.indices[j];
			H5Sselect_hyperslab(dataspace_id, j == 0 ? H5S_SELECT_SET : H5S_SELECT_OR, start, null, count,
				null);
		    }

		    Hyperslab slab1 = slabs[1];
		    count = dims.clone();
		    start = new long[rank];
		    count[slab1.rank - 1] = 1;
		    for (int j = 0, length = slab1.indices.length; j < length; j++) {
			start[slab1.rank - 1] = slab1.indices[j];
			H5Sselect_hyperslab(dataspace_id, H5S_SELECT_NOTB, start, null, count, null);
		    }

		    count = dims.clone();
		    start = new long[rank];
		    count[slab0.rank - 1] = 1;
		    for (int j = 0, length = slab0.indices.length; j < length; j++) {
			start[slab0.rank - 1] = slab0.indices[j];
			H5Sselect_hyperslab(dataspace_id, H5S_SELECT_XOR, start, null, count, null);
		    }

		}

	    }
	    long lsize = 1;
	    for (int i = 0; i < rank; i++) {
		lsize *= memSpaceDimensions[i];
	    }
	    tid = H5Dget_type(dataset_id);
	    int sdim = H5Tget_size(tid);

	    isString = memTypeId == H5T_FORTRAN_S1;
	    if (isString) {
		memTypeId = H5Tcopy(H5T_FORTRAN_S1);
		H5Tset_size(memTypeId, sdim);
	    }

	    if (memTypeId == -1) {
		int nativeTypeId = H5Tget_native_type(tid);
		boolean isNativeDatatype = H5Tequal(tid, nativeTypeId);

		H5Tclose(nativeTypeId);

		if (!isNativeDatatype) {
		    int tmptid = tid;
		    try {
			tid = H5Tget_native_type(tmptid);
		    } finally {
			try {
			    H5Tclose(tmptid);
			} catch (Exception e) {
			}
		    }
		}
	    }
	    Object data = null;
	    if (isString) {
		byte[] bytes = new byte[(int) (lsize * sdim)];
		H5Dread(dataset_id, memTypeId != -1 ? memTypeId : tid, mspace, dataspace_id, H5P_DEFAULT, bytes);
		data = Dataset.byteToString(bytes, sdim);

	    } else {
		if (lsize > Integer.MAX_VALUE) {
		    throw new ArrayTooBigException();
		}
		if (data == null) {
		    data = H5Datatype.allocateArray(memTypeId != -1 ? memTypeId : tid, (int) lsize);
		}

		if (data == null) {
		    throw new IllegalArgumentException("Unable to allocate array.");
		}
		H5Dread(dataset_id, memTypeId != -1 ? memTypeId : tid, mspace, dataspace_id, H5P_DEFAULT, data);
		if (data instanceof byte[]) {
		    data = Dataset.byteToString((byte[]) data, sdim);
		}
	    }
	    return data;
	} finally {
	    if (dataset_id != -1) {
		H5Dclose(dataset_id);
	    }
	    if (dataspace_id != -1) {
		H5Sclose(dataspace_id);
	    }
	    if (tid != -1) {
		H5Tclose(tid);
	    }
	    if (slabs != null && mspace != -1 && mspace != H5P_DEFAULT) {
		H5Sclose(mspace);
	    }
	    if (isString) {
		H5Tclose(memTypeId);
	    }
	}

    }

    /**
     * Read a subset of rows and columns from a gctx file.
     * 
     * @param file
     *            The gctx file.
     * @param rowIds
     *            The row ids or <tt>null</tt>.
     * @param columnIds
     *            The column ids or <tt>null</tt>.
     * @return The dataset
     */
    public cmap.matrix.Dataset read(Collection<String> rowIds, Collection<String> columnIds) {
	int[] columnIndices = null;
	if (columnIds != null) {
	    MetadataModel columnModel = this.getRowsOrColumns(false);
	    if (columnModel == null) {
		throw new NullPointerException("Column metadata is null");
	    }
	    Vector columnVector = columnModel.get("id");
	    if (columnVector == null) {
		columnVector = columnModel.get(0);
	    }
	    columnIndices = getIndicesForValues(columnVector, columnIds);
	}
	int[] rowIndices = null;
	if (rowIds != null) {
	    MetadataModel rowModel = this.getRowsOrColumns(true);
	    if (rowModel == null) {
		throw new NullPointerException("Row metadata is null");
	    }
	    Vector rowVector = rowModel.get("pr_id");
	    if (rowVector == null) {
		rowVector = rowModel.get("id");
	    }
	    rowIndices = getIndicesForValues(rowVector, rowIds);
	}
	cmap.matrix.Dataset dataset = this.read(rowIndices, columnIndices, true, true, false, false);

	return dataset;
    }

    private static int[] index(int[] indices) {
	IntIntPair[] pairs = new IntIntPair[indices.length];
	for (int i = 0, length = indices.length; i < length; i++) {
	    pairs[i] = new IntIntPair(indices[i], i);
	}
	Arrays.sort(pairs);
	int[] indexTable = new int[pairs.length];

	for (int i = 0, length = pairs.length; i < length; i++) {
	    indexTable[i] = pairs[i].index;
	}
	return indexTable;
    }

    public static class DatasetInfo {
	public int nrows, ncols, rowChunkSize, columnChunkSize;

	@Override
	public String toString() {
	    return "[nrows=" + nrows + ", ncols=" + ncols + ", rowChunkSize=" + rowChunkSize + ", columnChunkSize="
		    + columnChunkSize + "]";
	}
    }

    private static class IntIntPair implements Comparable<IntIntPair> {
	private int index;

	private int value;

	public IntIntPair(int value, int index) {
	    this.value = value;
	    this.index = index;
	}

	@Override
	public int compareTo(IntIntPair o) {
	    int thisVal = this.value;
	    int anotherVal = o.value;
	    return (thisVal < anotherVal ? -1 : (thisVal == anotherVal ? 0 : 1));
	}

	public int getValue() {
	    return value;
	}

	@Override
	public String toString() {
	    return value + ", " + index;
	}
    }

    private static class ArrayTooBigException extends Exception {

    }

    private static final int[] getIndicesForValues(Vector vector, Collection<?> values) {
	if (vector == null) {
	    throw new NullPointerException("vector is null");
	}
	if (values == null) {
	    throw new NullPointerException("values is null");
	}
	if (values.size() == 0) {
	    return null;
	}
	Map<Object, Integer> valueToIndex = new HashMap<Object, Integer>();
	for (int i = 0, size = vector.size(); i < size; i++) {
	    Object val = vector.getValue(i);
	    if (valueToIndex.containsKey(val)) {
		throw new IllegalArgumentException("Duplicate value: " + val + " in " + vector.getName());
	    }
	    valueToIndex.put(val, i);
	}
	int[] indices = new int[values.size()];
	int counter = 0;
	for (Object value : values) {
	    Integer index = valueToIndex.get(value);
	    if (index != null) {
		indices[counter++] = index;
	    } else {
		throw new IllegalArgumentException(value + " not found.");
	    }
	}

	return indices;
    }

    private static class Hyperslab {
	int[] indices;
	int rank;

	public Hyperslab(int[] indices, int rank) {
	    this.indices = indices;
	    this.rank = rank;
	}
    }

}
