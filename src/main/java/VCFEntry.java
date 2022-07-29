import java.util.zip.GZIPInputStream;
import java.io.*;
import java.util.*;

public class VCFEntry {
    String chrom;
    int pos;
    String id;
    String ref;
    String alt;
    String qual;
    String filter;
    LinkedHashMap<String, String> info = new LinkedHashMap<String, String>();
    String format;
    String[] genotypes;

    public VCFEntry(){
    }

    public VCFEntry(String[] src){
        chrom = src[0];
        pos = Integer.parseInt(src[1]);
        id = src[2];
        ref = src[3];
        alt = src[4];
        qual = src[5];
        filter = src[6];
        info = new LinkedHashMap<String, String>();
        for(String inf : src[7].split(";")){
            String[] x = inf.split("=");
            if(x.length==2){
                info.put(x[0], x[1]);
            }else {
                info.put(x[0], null);
            }
        }
        format = src[8];
        genotypes = new String[src.length-9];
        for(int i = 0; i<genotypes.length; i++){
            genotypes[i] = src[9+i];
        }
    }
}
