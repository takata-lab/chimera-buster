import java.util.*;
import java.io.*;
import java.util.zip.GZIPInputStream;

import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.assertEquals;

public class JointTest {
    @Test
    public void testMain(){
        try {
            String vcfPath = "../test_data/ML-BA001.Maryland_DES_wave1.hg38.BQSR.merge.markdup.tagfix.mutect.GENCODEv32_RefSeq.slop100.PASS.allannotation.vcf.gz";
            String bamPath = "../test_data/ML-BA001.Maryland_DES_wave1.hg38.BQSR.merge.markdup.tagfix.mutect.GENCODEv32_RefSeq.slop100.PASS.slop300.bamout.bam";
            String[] args = new String[4];
            args[0] = "--vcf";
            args[1] = vcfPath;
            args[2] = "--bam";
            args[3] = bamPath;
            // run main()
            Main.main(args);
        }catch(Exception e){
            e.printStackTrace();
        }
    }
}
