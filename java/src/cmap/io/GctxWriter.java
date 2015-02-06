package cmap.io;

import static ncsa.hdf.hdf5lib.H5.H5Aclose;
import static ncsa.hdf.hdf5lib.H5.H5Acreate;
import static ncsa.hdf.hdf5lib.H5.H5Aexists_by_name;
import static ncsa.hdf.hdf5lib.H5.H5Awrite;
import static ncsa.hdf.hdf5lib.H5.H5Dclose;
import static ncsa.hdf.hdf5lib.H5.H5Dcreate;
import static ncsa.hdf.hdf5lib.H5.H5Dget_space;
import static ncsa.hdf.hdf5lib.H5.H5Dwrite;
import static ncsa.hdf.hdf5lib.H5.H5Dwrite_float;
import static ncsa.hdf.hdf5lib.H5.H5Fclose;
import static ncsa.hdf.hdf5lib.H5.H5Fcreate;
import static ncsa.hdf.hdf5lib.H5.H5Fopen;
import static ncsa.hdf.hdf5lib.H5.H5Gclose;
import static ncsa.hdf.hdf5lib.H5.H5Gcreate;
import static ncsa.hdf.hdf5lib.H5.H5Gopen;
import static ncsa.hdf.hdf5lib.H5.H5Ldelete;
import static ncsa.hdf.hdf5lib.H5.H5Pclose;
import static ncsa.hdf.hdf5lib.H5.H5Pcreate;
import static ncsa.hdf.hdf5lib.H5.H5Pset_chunk;
import static ncsa.hdf.hdf5lib.H5.H5Pset_create_intermediate_group;
import static ncsa.hdf.hdf5lib.H5.H5Pset_deflate;
import static ncsa.hdf.hdf5lib.H5.H5Pset_libver_bounds;
import static ncsa.hdf.hdf5lib.H5.H5Sclose;
import static ncsa.hdf.hdf5lib.H5.H5Screate_simple;
import static ncsa.hdf.hdf5lib.H5.H5Sselect_hyperslab;
import static ncsa.hdf.hdf5lib.H5.H5Tclose;
import static ncsa.hdf.hdf5lib.H5.H5Tcopy;
import static ncsa.hdf.hdf5lib.H5.H5Tset_size;
import static ncsa.hdf.hdf5lib.HDF5Constants.H5F_ACC_EXCL;
import static ncsa.hdf.hdf5lib.HDF5Constants.H5F_ACC_RDWR;
import static ncsa.hdf.hdf5lib.HDF5Constants.H5F_LIBVER_LATEST;
import static ncsa.hdf.hdf5lib.HDF5Constants.H5P_DATASET_CREATE;
import static ncsa.hdf.hdf5lib.HDF5Constants.H5P_DEFAULT;
import static ncsa.hdf.hdf5lib.HDF5Constants.H5P_FILE_ACCESS;
import static ncsa.hdf.hdf5lib.HDF5Constants.H5P_LINK_CREATE;
import static ncsa.hdf.hdf5lib.HDF5Constants.H5S_ALL;
import static ncsa.hdf.hdf5lib.HDF5Constants.H5S_SELECT_SET;
import static ncsa.hdf.hdf5lib.HDF5Constants.H5S_UNLIMITED;
import static ncsa.hdf.hdf5lib.HDF5Constants.H5T_FORTRAN_S1;
import static ncsa.hdf.hdf5lib.HDF5Constants.H5T_IEEE_F32BE;
import static ncsa.hdf.hdf5lib.HDF5Constants.H5T_NATIVE_FLOAT;
import static ncsa.hdf.hdf5lib.HDF5Constants.H5T_NATIVE_INT;
import static ncsa.hdf.hdf5lib.HDF5Constants.H5T_STD_I32BE;

import java.io.File;
import java.net.URL;
import java.util.Date;

import ncsa.hdf.hdf5lib.HDF5Constants;
import ncsa.hdf.hdf5lib.exceptions.HDF5Exception;
import cmap.matrix.Dataset;
import cmap.matrix.MetadataModel;
import cmap.matrix.TransposedDatasetView;
import cmap.matrix.Vector;

/**
 * 
 * Write files in gctx format.
 *
 */
public class GctxWriter {

    private static final String COL_METADATA_PATH = "/0/META/COL";
    private static final String DATA_MATRIX_PATH = "/0/DATA/0/matrix";
    private static final String ROW_METADATA_PATH = "/0/META/ROW";
    private static final String VERSION = "GCTX1.0";
    private static final String VERSION_ATTR = "version";
    private long columnDataChunkSize = 200;
    private int file_id = -1;
    private boolean gzip = false;
    private int matrixDataType = H5T_IEEE_F32BE;
    private long rowDataChunkSize = 100;
    private boolean writeInSlabs = true;

    public GctxWriter(String file, Dataset dataset) throws Exception {
	try {
	    int fapl = H5Pcreate(H5P_FILE_ACCESS);
	    H5Pset_libver_bounds(fapl, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST);

	    if (!new File(file).exists()) {
		file_id = H5Fcreate(file, H5F_ACC_EXCL, H5P_DEFAULT, fapl);
	    } else {
		file_id = H5Fopen(file, H5F_ACC_RDWR, fapl);
	    }

	    if (!H5Aexists_by_name(file_id, "/", VERSION_ATTR, H5P_DEFAULT)) {
		int SDIM = VERSION.length();
		int filetype = H5Tcopy(H5T_FORTRAN_S1);
		H5Tset_size(filetype, SDIM);
		int memtype = H5Tcopy(H5T_FORTRAN_S1);
		H5Tset_size(memtype, SDIM);

		long[] dims = new long[] { 1 };
		int space = H5Screate_simple(1, dims, null);

		int dset = H5Acreate(file_id, VERSION_ATTR, filetype, space, H5P_DEFAULT, H5P_DEFAULT);
		byte[] byteArray = VERSION.getBytes();
		H5Awrite(dset, memtype, byteArray);
		H5Aclose(dset);
		H5Sclose(space);
		H5Tclose(filetype);
		H5Tclose(memtype);
	    }
	    writeDataset(DATA_MATRIX_PATH, dataset);

	} finally {
	    close();
	}
    }

    private void close() {
	try {
	    if (file_id != -1) {
		H5Fclose(file_id);
		// H5close();
		file_id = -1;
	    }
	} catch (Exception e) {
	}
    }

    private boolean deleteIfExists(String name) {
	try {
	    H5Ldelete(file_id, name, HDF5Constants.H5P_DEFAULT);
	    return true;
	} catch (Exception e) {
	    return false;
	}
    }

    private void write1DArray(String path, Object array, int length, int fileType, int memType, long chunkSize)
	    throws HDF5Exception {
	int mid = H5Screate_simple(1, new long[] { length }, null);
	int dcpl = H5Pcreate(H5P_DATASET_CREATE);
	int dset = H5Dcreate(file_id, path, fileType, mid, H5P_DEFAULT, dcpl, H5P_DEFAULT);
	if (chunkSize > 0) {
	    H5Pset_chunk(dcpl, 1, new long[] { chunkSize });
	    H5Pset_deflate(dcpl, 6);
	}

	int fid = H5Dget_space(dset);

	H5Dwrite(dset, memType, mid, fid, H5P_DEFAULT, array);

	H5Dclose(dset);
	H5Sclose(mid);
	H5Sclose(fid);
	H5Pclose(dcpl);
    }

    private void writeDataMatrix(String path, final cmap.matrix.Dataset dataset, boolean rowMajorOrder,
	    long rowDataChunkSize, long columnDataChunkSize) throws Exception {

	deleteIfExists(path);
	int gcpl = H5Pcreate(H5P_LINK_CREATE);
	H5Pset_create_intermediate_group(gcpl, true);
	int slash = path.lastIndexOf('/');
	String groupPath = slash == -1 ? path : path.substring(0, slash);
	int group = H5Gcreate(file_id, groupPath, gcpl, H5P_DEFAULT, H5P_DEFAULT);

	H5Pclose(gcpl);
	H5Gclose(group);

	cmap.matrix.Dataset datasetForDataVariable = rowMajorOrder ? dataset : new TransposedDatasetView(dataset);
	long[] dims = new long[] { datasetForDataVariable.getRowCount(), datasetForDataVariable.getColumnCount() };
	long[] maxDims = new long[] { H5S_UNLIMITED, datasetForDataVariable.getColumnCount() }; // columns
												// can
												// grow

	long[] chunks = rowMajorOrder ? new long[] { rowDataChunkSize, columnDataChunkSize } : new long[] {
		columnDataChunkSize, rowDataChunkSize };

	int mid = H5Screate_simple(2, dims, maxDims);
	int dcpl = H5Pcreate(H5P_DATASET_CREATE);
	H5Pset_chunk(dcpl, 2, chunks);
	if (gzip) {
	    H5Pset_deflate(dcpl, 6);
	}

	int dset = H5Dcreate(file_id, path, matrixDataType, mid, H5P_DEFAULT, dcpl, H5P_DEFAULT);

	int fid = H5Dget_space(dset);

	if (writeInSlabs) {
	    long[] stride = new long[] { 1, 1 };
	    long[] count = new long[] { 1, 1 };
	    long[] block = new long[] { 1, 1 };
	    long[] start = new long[] { 0, 0 };
	    float[] array = new float[1];

	    H5Sselect_hyperslab(mid, H5S_SELECT_SET, start, stride, count, block);

	    for (int i = 0, nrows = datasetForDataVariable.getRowCount(); i < nrows; i++) {
		start[0] = i;

		for (int j = 0, ncols = datasetForDataVariable.getColumnCount(); j < ncols; j++) {
		    start[1] = j;
		    array[0] = datasetForDataVariable.getValue(i, j);
		    H5Sselect_hyperslab(fid, H5S_SELECT_SET, start, stride, count, block);
		    H5Dwrite_float(dset, H5T_NATIVE_FLOAT, mid, fid, H5P_DEFAULT, array);
		    // H5Dwrite_float(dset, H5T_NATIVE_FLOAT, mid, fid,
		    // H5P_DEFAULT, array);

		}
	    }
	} else {
	    float[] array = new float[datasetForDataVariable.getRowCount() * datasetForDataVariable.getColumnCount()];
	    for (int i = 0, nrows = datasetForDataVariable.getRowCount(); i < nrows; i++) {
		for (int j = 0, ncols = datasetForDataVariable.getColumnCount(); j < ncols; j++) {
		    array[i * ncols + j] = datasetForDataVariable.getValue(i, j);
		}
	    }

	    H5Dwrite_float(dset, H5T_NATIVE_FLOAT, mid, fid, H5P_DEFAULT, array);
	}

	H5Dclose(dset);
	H5Sclose(mid);
	H5Sclose(fid);
	H5Pclose(dcpl);

    }

    private void writeDataset(String path, final cmap.matrix.Dataset dataset) throws Exception {

	long columnChunk = columnDataChunkSize;
	if (columnDataChunkSize <= 0 || columnDataChunkSize > dataset.getColumnCount()) {
	    columnChunk = dataset.getColumnCount();
	}
	long rowChunk = rowDataChunkSize;
	if (rowDataChunkSize <= 0 || rowDataChunkSize > dataset.getRowCount()) {
	    rowChunk = dataset.getRowCount();
	}
	writeDataMatrix(path, dataset, false, rowChunk, columnChunk);

	MetadataModel columnMetadata = dataset.getColumnMetadata();
	if (columnMetadata.getMetadataCount() > 0) {
	    for (int i = 0, size = columnMetadata.getMetadataCount(); i < size; i++) {
		Vector v = columnMetadata.get(i);
		writeMetadata(true, columnMetadata.getColumnName(i), v, v.getColumnClass());
	    }
	}

	MetadataModel rowMetadata = dataset.getRowMetadata();
	if (rowMetadata.getMetadataCount() > 0) {
	    for (int i = 0, size = rowMetadata.getMetadataCount(); i < size; i++) {
		Vector v = rowMetadata.get(i);
		writeMetadata(false, rowMetadata.getColumnName(i), v, v.getColumnClass());
	    }
	}
    }

    private void writeMetadata(boolean isColumnMetadata, String annotationName, Vector values, Class<?> cl)
	    throws Exception {
	String path = (isColumnMetadata ? COL_METADATA_PATH : ROW_METADATA_PATH) + "/";
	try {
	    int gid = H5Gopen(file_id, path, HDF5Constants.H5P_DEFAULT);
	    H5Gclose(gid);
	} catch (Exception x) {
	    int gcpl = H5Pcreate(H5P_LINK_CREATE);
	    H5Pset_create_intermediate_group(gcpl, true);
	    int group = H5Gcreate(file_id, path, gcpl, H5P_DEFAULT, H5P_DEFAULT);
	    H5Pclose(gcpl);
	    H5Gclose(group);
	}

	if ("".equals(annotationName)) {
	    annotationName = "annot";
	}
	annotationName = annotationName.replace('/', '_');

	if (String.class.isAssignableFrom(cl)) {
	    writeString(path + annotationName, values, values.size());
	} else if (Float.class.isAssignableFrom(cl)) {
	    float[] tmp = new float[values.size()];
	    for (int i = 0, nrows = tmp.length; i < nrows; i++) {
		Number val = (Number) values.getValue(i);
		tmp[i] = val == null ? Float.NaN : val.floatValue();
	    }

	    write1DArray(path + annotationName, tmp, tmp.length, H5T_IEEE_F32BE, H5T_NATIVE_FLOAT, tmp.length);
	} else if (Integer.class.isAssignableFrom(cl)) {
	    int[] tmp = new int[values.size()];
	    for (int i = 0, nrows = tmp.length; i < nrows; i++) {
		Integer val = (Integer) values.getValue(i);
		tmp[i] = val == null ? -Integer.MAX_VALUE : val;
	    }
	    write1DArray(path + annotationName, tmp, tmp.length, H5T_STD_I32BE, H5T_NATIVE_INT, tmp.length);

	} else {
	    writeString(path + annotationName, values, values.size());
	}

    }

    private void writeString(String path, Vector vector, long chunkSize) throws NullPointerException, HDF5Exception {
	int SDIM = maxStringLength(vector);
	if (SDIM == 0) {
	    System.out.println("All values are empty for " + vector.getName() + " vector");
	    return;
	}
	int filetype = H5Tcopy(H5T_FORTRAN_S1);
	H5Tset_size(filetype, SDIM);
	int memtype = H5Tcopy(H5T_FORTRAN_S1);
	H5Tset_size(memtype, SDIM);
	// H5Tset_strpad(memtype, H5T_STR_NULLPAD);
	int dcpl = H5Pcreate(H5P_DATASET_CREATE);

	long[] dims = new long[] { vector.size() };
	int space = H5Screate_simple(1, dims, null);
	if (chunkSize > 0) {
	    if (chunkSize > vector.size()) {
		chunkSize = vector.size();
	    }
	    H5Pset_chunk(dcpl, 1, new long[] { chunkSize });
	    H5Pset_deflate(dcpl, 6);
	}

	int dset = H5Dcreate(file_id, path, filetype, space, H5P_DEFAULT, dcpl, H5P_DEFAULT);

	byte[][] byteArray = new byte[vector.size()][];
	for (int i = 0, nrows = vector.size(); i < nrows; i++) {
	    // start[0] = i * SDIM;
	    // H5Sselect_hyperslab(space, H5S_SELECT_SET, start, stride, count,
	    // block);
	    byteArray[i] = getBytes((String) vector.getValue(i), new byte[SDIM]);

	}
	H5Dwrite(dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, byteArray);
	H5Dclose(dset);
	H5Sclose(space);
	H5Tclose(filetype);
	H5Tclose(memtype);
	H5Pclose(dcpl);

    }

    private static byte[] getBytes(String s, byte[] array) {
	if (s == null) {
	    s = "";
	}
	int slength = s.length();
	for (int i = 0; i < slength; i++) {
	    array[i] = (byte) s.charAt(i);
	}
	for (int i = slength, length = array.length; i < length; i++) {
	    array[i] = 0;
	}
	return array;
    }

    private static ValueEncoder getEncoder(Class<?> c) {
	if (Date.class.isAssignableFrom(c)) {
	    return new DateEncoder();
	} else if (URL.class.isAssignableFrom(c)) {
	    return new URLEncoder();
	}
	return new NullEncoder();
    }

    private static int maxStringLength(Vector vector) {
	int maxLength = 0;
	ValueEncoder encoder = getEncoder(vector.getColumnClass());
	for (int i = 0, nrows = vector.size(); i < nrows; i++) {
	    Object val = vector.getValue(i);
	    if (val != null) {
		val = encoder.encode(val);
		int length = val.toString().length();
		if (length > maxLength) {
		    maxLength = length;
		}
	    }
	}
	return maxLength;
    }

    private static class DateEncoder implements ValueEncoder {

	@Override
	public Object encode(Object value) {
	    return value != null ? ((Date) value).getTime() : null;
	}

    }

    private static class NullEncoder implements ValueEncoder {

	@Override
	public Object encode(Object value) {
	    return value;
	}

    }

    private static class URLEncoder implements ValueEncoder {

	@Override
	public Object encode(Object value) {
	    return value != null ? ((URL) value).toString() : null;
	}

    }

    private static interface ValueEncoder {
	public Object encode(Object value);
    }

}
