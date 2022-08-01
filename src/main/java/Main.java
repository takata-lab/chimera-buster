import java.io.File;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.FileNotFoundException;
import java.util.zip.GZIPInputStream;
import java.util.ArrayList;
import java.util.List;

import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.SamLocusIterator;
import htsjdk.samtools.util.AbstractLocusInfo;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReaderFactory;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;


public class Main {
    static IndexedFastaSequenceFile refGenome;
    private static char getRefBase(String chrom, int pos){
        ReferenceSequence seq = refGenome.getSubsequenceAt(chrom, pos, pos);
        return seq.getBaseString().charAt(0);
    }
    private static List<SamLocusIterator.RecordAndOffset> collectAltAlignment(List<SamLocusIterator.RecordAndOffset> src, String ref, String alt){
        ArrayList<SamLocusIterator.RecordAndOffset> result = new ArrayList<>();
        for(SamLocusIterator.RecordAndOffset alignment: src){
            SAMRecord rec = alignment.getRecord();
            if(alignment.getReadBase() == alt.charAt(0)){
                result.add(alignment);
            }
        }
        return result;
    }
    private static List<SamLocusIterator.RecordAndOffset> collectNonRefAlignment(List<SamLocusIterator.RecordAndOffset> src, String ref, String alt){
        ArrayList<SamLocusIterator.RecordAndOffset> result = new ArrayList<>();
        for(SamLocusIterator.RecordAndOffset alignment: src){
            SAMRecord rec = alignment.getRecord();
            if(alignment.getReadBase() != ref.charAt(0)){
                result.add(alignment);
            }
        }
        return result;
    }
    // filer out ArtificialHaplotypeRG ReadGroup created by GATK
    private static List<SamLocusIterator.RecordAndOffset> filterAlignment(List<SamLocusIterator.RecordAndOffset> src){
        ArrayList<SamLocusIterator.RecordAndOffset> result = new ArrayList<>();
        for(SamLocusIterator.RecordAndOffset alignment: src){
            SAMRecord rec = alignment.getRecord();
            // System.err.println("RGID:" + rec.getReadGroup().getId());
            if(rec.getReadGroup().getId().equals("ArtificialHaplotypeRG")){
                continue;
            }else {
                result.add(alignment);
            }
        }
        return result;
    }
    public static void analyseVCF(ArrayList<VCFEntry> vcf, String bamPath) throws IOException{
        SamLocusIterator locus = null;
        final SamReader reader = SamReaderFactory.makeDefault().open(new File(bamPath));

        for(VCFEntry snp: vcf){
            // TODO: validate snp
            System.err.println(snp.chrom + ":" + snp.pos + " " + snp.ref + " -> " + snp.alt);
            IntervalList region = createTargetRegionInterval(reader, snp);
            locus = new SamLocusIterator(reader, region);
            System.err.println("checking "+ snp.chrom + ":" + snp.pos);
            while(locus.hasNext()){
                AbstractLocusInfo cur = locus.next();
                List<SamLocusIterator.RecordAndOffset> alignments = filterAlignment(cur.getRecordAndOffsets());
                // System.err.println("pos: " + cur.getSequenceName() + ":" + cur.getPosition() + " " + (snp.pos == cur.getPosition()));
                if(cur.getPosition() == snp.pos){
                    List<SamLocusIterator.RecordAndOffset> alts = collectAltAlignment(alignments, snp.ref, snp.alt);
                }else {
                    byte refBase = (byte)getRefBase(snp.chrom, snp.pos);
                    List<SamLocusIterator.RecordAndOffset> alts = collectAltAlignment(alignments, snp.ref, snp.alt);
                }
            }
            locus.close();
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
            }else if(argv[i].equals("--ref")){
                try {
                    refGenome = new IndexedFastaSequenceFile(new File(argv[i+1]));
                }catch(FileNotFoundException e){
                    System.err.println("Reference genome file: " + argv[i+1] + " was not found");
                    System.exit(1);
                }
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
        if(refGenome == null){
            System.err.println("No Reference genome file is given");
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
        System.err.println("java -jar ChimeraBuster.jar --vcf [vcf file] --bam [bam file] --ref [ref genome fasta file]");
    }
}
