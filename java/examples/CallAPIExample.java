import l1ktools.cmap.GctxDataset;
import l1ktools.cmap.GctxReader;
import java.io.*;
import java.net.*;

public class CallAPIExample {
    /**
     *  return a list of cell lines with a TP53 mutation:
     */
    final static String CELL_LINES_WITH_TP53_MUTATION = "http://api.lincscloud.org/a2/cellinfo?q={%22mutations%22:%22TP53%22}&d=cell_id&user_key=lincsdemo";

    public static void main(String[] args) {
       BufferedReader reader = null;
       try {
            final URL url = new URL(CELL_LINES_WITH_TP53_MUTATION);

            final URLConnection urlConnection = url.openConnection();

            urlConnection.setAllowUserInteraction(false);
            reader = new BufferedReader(new InputStreamReader(urlConnection
                    .getInputStream()));
            String currentLine = null;
            while ((currentLine=reader.readLine())!=null) {
                System.out.println(currentLine);
            }
          } catch (Exception e) {
            e.printStackTrace();
        } finally {
           if(reader != null){
               try{
               reader.close();
               }catch (Exception ee) {

               }
           }
        }

    }
}
