import l1ktools.cmap.GctxDataset;
import l1ktools.cmap.GctxReader;

public class ReadGctxExample {

    public static void main(String[] args) {
	GctxReader reader = new GctxReader("../data/modzs_n272x978.gctx");
	try {
	    GctxDataset dataset = reader.read();
	    System.out.println(dataset.getRowCount() + " rows, " + dataset.getColumnCount() + " columns");
	    for (String key : dataset.getColumnMetadata().keySet()) {
		System.out.println(key + "=" + dataset.getColumnMetadata().get(key));
	    }
	} catch (Exception e) {
	    e.printStackTrace();
	} finally {
	    reader.close();
	}

    }
}
