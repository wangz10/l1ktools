import java.util.ArrayList;
import java.util.Collection;

import cmap.io.GctxReader;
import cmap.matrix.Dataset;
import cmap.matrix.Vector;

public class ReadGctxExample {

    public static void main(String[] args) throws Exception {
	GctxReader reader = new GctxReader("../data/modzs_n272x978.gctx");
	try {
	    // read the full dataset
	    Dataset dataset = reader.read();
	    System.out.println(dataset.getRowCount() + " rows, " + dataset.getColumnCount() + " columns");
	    Vector columnIds = dataset.getColumnMetadata().get("id");

	    // read the 1st 10 column ids
	    Collection<String> columnIdsToExtract = new ArrayList<String>();

	    for (int j = 0; j < 10; j++) {
		columnIdsToExtract.add((String) columnIds.getValue(j));
	    }
	    dataset = reader.read(null, columnIdsToExtract);
	    System.out.println(dataset.getRowCount() + " rows, " + dataset.getColumnCount() + " columns");

	} finally {
	    reader.close();
	}

    }
}
