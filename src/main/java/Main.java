import java.io.File;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.zip.GZIPInputStream;
import java.util.ArrayList;

import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.SamLocusIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;


public class Main {
    public static void analyseVCF(ArrayList<VCFEntry> vcf, String bamPath) throws IOException{
        SamLocusIterator locus = null;
        final SamReader reader = SamReaderFactory.makeDefault().open(new File(bamPath));

        for(VCFEntry snp: vcf){
            // TODO: validate snp
            System.err.println(snp.chrom + ":" + snp.pos + " " + snp.ref + " -> " + snp.alt);
            IntervalList region = createTargetRegionInterval(reader, snp);
            locus = new SamLocusIterator(reader, region);

        }
    }
    private static IntervalList createTargetRegionInterval(SamReader reader, VCFEntry snp){
        IntervalList list = new IntervalList(reader.getFileHeader());
        int start = snp.pos - 10;
        int end = snp.pos + 10;
        
        Interval interval = new Interval(snp.chrom, start, end);
        list.add(interval);

        return list;
    }
    public static void main(String[] argv){
        String vcfPath = null;
        String bamPath = null;
        for(int i = 0; i<argv.length-1; i++){
            if(argv[i].equals("--vcf")){
                vcfPath = argv[i+1];
            }else if(argv[i].equals("--bam")){
                bamPath = argv[i+1];
            }
        }
        if(vcfPath == null){
            System.err.println("No VCF file is given");
            printUsage();
            System.exit(1);
        }
        if(bamPath == null){
            System.err.println("No BAM file is given");
            printUsage();
            System.exit(1);
        }
        try {
            ArrayList<VCFEntry> vcf = loadVCF(vcfPath);
            analyseVCF(vcf, bamPath);
        }catch(Exception e){
            e.printStackTrace();
        }
    }
    public static ArrayList<VCFEntry> loadVCF(String vcfPath) throws IOException{
        ArrayList<VCFEntry> entries = new ArrayList<>();
        BufferedReader br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(vcfPath))));
        String raw = null;
        while(null != (raw = br.readLine())){
            if(raw.startsWith("#")){
                continue;
            }
            VCFEntry ent = new VCFEntry(raw.split("\t"));
            entries.add(ent);
        }
        br.close();
        return entries;
    }
    public static final void printUsage(){
        System.err.println("java -jar ChimeraBuster.jar --vcf [vcf file] --bam [bam file]");
    }
}
