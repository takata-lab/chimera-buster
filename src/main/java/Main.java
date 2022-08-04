import java.io.File;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.HashSet;
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
    private static void classifyAlignment(List<SamLocusIterator.RecordAndOffset> src, String refBase, String altBase, HashMap<String, SamLocusIterator.RecordAndOffset> refs, HashMap<String, SamLocusIterator.RecordAndOffset> alts){
        for(SamLocusIterator.RecordAndOffset alignment: src){
            SAMRecord rec = alignment.getRecord();
            if(alignment.getReadBase() == altBase.charAt(0)){
                alts.put(rec.getReadName(), alignment);
            }else if(alignment.getReadBase() == refBase.charAt(0)){
                refs.put(rec.getReadName(), alignment);
            }
        }
    }
    /*
    private static List<SamLocusIterator.RecordAndOffset> collectNonRefAlignment(List<SamLocusIterator.RecordAndOffset> src, String ref, String alt){
        ArrayList<SamLocusIterator.RecordAndOffset> result = new ArrayList<>();
        for(SamLocusIterator.RecordAndOffset alignment: src){
            SAMRecord rec = alignment.getRecord();
            if(alignment.getReadBase() != ref.charAt(0)){
                result.add(alignment);
            }
        }
        return result;
    }*/
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
    private static String calcScore(Object[] ref, Object[] alt){
        HashMap<String, SamLocusIterator.RecordAndOffset> ref1 = (HashMap<String, SamLocusIterator.RecordAndOffset>)ref[0];
        HashMap<String, SamLocusIterator.RecordAndOffset> ref2 = (HashMap<String, SamLocusIterator.RecordAndOffset>)ref[1];
        HashMap<String, SamLocusIterator.RecordAndOffset> alt1 = (HashMap<String, SamLocusIterator.RecordAndOffset>)alt[0];
        HashMap<String, SamLocusIterator.RecordAndOffset> alt2 = (HashMap<String, SamLocusIterator.RecordAndOffset>)alt[1];
        // filter uncommon reads: remove reads not covering both position
        HashSet<String> exclude = new HashSet<>();
        for(String name: alt1.keySet()){
            if(!ref2.containsKey(name) && !alt2.containsKey(name)){
                exclude.add(name);
            }
        }
        for(String name: alt2.keySet()){
            if(!ref1.containsKey(name) && !alt1.containsKey(name)){
                exclude.add(name);
            }
        }
        for(String name: exclude){
            alt1.remove(name);
            alt2.remove(name);
        }
        // forward
        int hit1 = 0;
        int hit2 = 0;
        for(String name: alt1.keySet()){
            if(alt2.containsKey(name)){
                hit1++;
            }
        }
        for(String name: alt2.keySet()){
            if(alt1.containsKey(name)){
                hit2++;
            }
        }
        double rate1 = ((double)hit1)/alt1.size();
        double rate2 = ((double)hit2)/alt2.size();

        String result = join("Forward=", String.valueOf(hit1), "/", String.valueOf(alt1.size()), ":", String.valueOf(rate1));
        result = join(result, ";Reverse=", String.valueOf(hit2), "/", String.valueOf(alt2.size()), ":", String.valueOf(rate2));
        return result;
    }
    public static String analyzeSNPPair(VCFEntry[] snps, SamReader reader){
        // classify ref reads and alt reads
        Object[] ref = new Object[2];
        Object[] alt = new Object[2];
        for(int i = 0; i<2; i++){
            ref[i] = new HashMap<String, SamLocusIterator.RecordAndOffset>();
            alt[i] = new HashMap<String, SamLocusIterator.RecordAndOffset>();
            VCFEntry snp = snps[i];
            IntervalList region = createTargetRegionInterval(reader, snp);
            SamLocusIterator locus = new SamLocusIterator(reader, region);
            while(locus.hasNext()){
                AbstractLocusInfo cur = locus.next();
                List<SamLocusIterator.RecordAndOffset> alignments = filterAlignment(cur.getRecordAndOffsets());
                // System.err.println("pos: " + cur.getSequenceName() + ":" + cur.getPosition() + " " + (snp.pos == cur.getPosition()));
                // System.err.print("   " + cur.getSequenceName() + " " + cur.getPosition());
                if(cur.getPosition() == snp.pos){
                    classifyAlignment(alignments, snp.ref, snp.alt, (HashMap<String, SamLocusIterator.RecordAndOffset>) ref[i], (HashMap<String, SamLocusIterator.RecordAndOffset>) alt[i]);
                }
            }
            locus.close();
        }
        return calcScore(ref, alt);
    }
    public static void analyseVCF(ArrayList<VCFEntry> vcf, String bamPath) throws IOException{
        final SamReader reader = SamReaderFactory.makeDefault().open(new File(bamPath));

        VCFEntry prev = null;
        for(VCFEntry snp: vcf){
            System.err.println(snp.chrom + ":" + snp.pos + " " + snp.ref + " -> " + snp.alt);
            if(prev != null && (snp.pos - prev.pos) < 10){
                System.err.println("checking "+ snp.chrom + ":" + snp.pos);
                VCFEntry[] pair = {prev, snp};
                String scores = analyzeSNPPair(pair, reader);
                snp.appendInfo(scores);
            }
            prev = snp;
            System.out.println(snp.toString());
        }
    }
    private static IntervalList createTargetRegionInterval(SamReader reader, VCFEntry snp){
        IntervalList list = new IntervalList(reader.getFileHeader());
        int start = snp.pos;
        int end = snp.pos;
        
        Interval interval = new Interval(snp.chrom, start, end);
        list.add(interval);

        return list;
    }
    public static String join(String... list){
        StringBuilder buf = new StringBuilder(list[0]);
        for(int i = 1; i<list.length; i++){
            buf.append(list[i]);
        }
        return buf.toString();
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
