import java.util.*;
import java.io.*;
import java.util.zip.GZIPInputStream;

import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.assertEquals;

public class VCFTest {
    @Test
    public void testVCFParser(){
        try {
            String vcfPath = "../test_data/ML-BA001.Maryland_DES_wave1.hg38.BQSR.merge.markdup.tagfix.mutect.GENCODEv32_RefSeq.slop100.PASS.allannotation.vcf.gz";
            BufferedReader br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(vcfPath))));
            String raw = null;
            while(null != (raw = br.readLine())){
                if(raw.startsWith("#")){
                    continue;
                }
                VCFEntry ent = new VCFEntry(raw.split("\t"));
            }
            br.close();
        }catch(Exception e){
            e.printStackTrace();
        }
    }
}
