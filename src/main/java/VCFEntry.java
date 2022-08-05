import java.util.LinkedHashMap;

public class VCFEntry {
    String chrom;
    int pos;
    String id;
    String ref;
    String alt;
    String qual;
    String filter;
    LinkedHashMap<String, String> info = new LinkedHashMap<String, String>();
    String rawInfo;
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
        rawInfo = src[7];
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
    public void appendInfo(String app){
        rawInfo = rawInfo + ";" + app;
    }
    public String toString(){
        StringBuilder buf = new StringBuilder();
        buf.append(chrom);
        buf.append("\t");
        buf.append(String.valueOf(pos));
        buf.append("\t");
        buf.append(id);
        buf.append("\t");
        buf.append(ref);
        buf.append("\t");
        buf.append(alt);
        buf.append("\t");
        buf.append(qual);
        buf.append("\t");
        buf.append(filter);
        buf.append("\t");
        buf.append(rawInfo);
        buf.append("\t");
        buf.append(format);
        for(int i = 0; i<genotypes.length; i++){
            buf.append("\t");
            buf.append(String.valueOf(genotypes[i]));
        }
        return buf.toString();
    }
}
