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
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.filter.SecondaryAlignmentFilter;
import htsjdk.samtools.filter.DuplicateReadFilter;
import htsjdk.samtools.CigarOperator;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;


public class Main2 {
    static IndexedFastaSequenceFile refGenome;
    static boolean debug = false; // 83520057
    private static char getRefBase(String chrom, int pos){
        ReferenceSequence seq = refGenome.getSubsequenceAt(chrom, pos, pos);
        return seq.getBaseString().charAt(0);
    }
    private static void classifyAlignment(List<SamLocusIterator.RecordAndOffset> src, String refBase, String altBase, HashMap<String, SamLocusIterator.RecordAndOffset> refs, HashMap<String, SamLocusIterator.RecordAndOffset> alts){
        int total = 0;
        int alt = 0;
        int ref = 0;
        int suppl = 0;
        for(SamLocusIterator.RecordAndOffset alignment: src){
            total++;
            if(alignment.getRecord().isSecondaryOrSupplementary()){
            // if(alignment.getRecord().getSupplementaryAlignmentFlag()){
                suppl++;
                System.err.println("suppl:" + alignment.getRecord().getReadName() + " " + alignment.getReadBase());
            }
            SAMRecord rec = alignment.getRecord();
            if(alignment.getReadBase() == altBase.charAt(0)){
                alts.put(rec.getReadName(), alignment);
                alt++;
            }else if(alignment.getReadBase() == refBase.charAt(0)){
                refs.put(rec.getReadName(), alignment);
                ref++;
            }
        }
        if(debug){
            System.err.println("refbase: " + refBase);
            System.err.println("altbase: " + altBase);
            System.err.println("total: " + total);
            System.err.println("ref: " + ref);
            System.err.println("alt: " + alt);
            System.err.println("suppl: " + suppl);
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
    /*
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
    }*/
    private static String calcScore(HashMap<String, SamLocusIterator.RecordAndOffset> ref, HashMap<String, SamLocusIterator.RecordAndOffset> alt){
        // filter uncommon reads: remove reads not covering both position
        int ref_count = ref.size();
        int alt_count = alt.size();

        int alt_softclip_count = 0;
        int alt_hardclip_count = 0;
        for(SamLocusIterator.RecordAndOffset cur: alt.values()){
            if( cur.getRecord().getCigar().getFirstCigarElement().getOperator().equals(CigarOperator.S)){
                alt_softclip_count++;
            }else if( cur.getRecord().getCigar().getFirstCigarElement().getOperator().equals(CigarOperator.H)){
                alt_hardclip_count++;
            }
        }

        double rate_soft = ((double)alt_softclip_count)/alt_count;
        double rate_hard = ((double)alt_hardclip_count)/alt_count;

        String result = join("SoftClips=", String.valueOf(alt_softclip_count), "/", String.valueOf(alt_count), ":", String.valueOf(rate_soft));
        result = join(result, ";HardClips=", String.valueOf(alt_hardclip_count), "/", String.valueOf(alt_count), ":", String.valueOf(rate_hard));
        return result;
    }
    public static String analyzeSNP(VCFEntry snp, SamReader reader){
        // classify ref reads and alt reads
        HashMap<String, SamLocusIterator.RecordAndOffset> ref = new HashMap<String, SamLocusIterator.RecordAndOffset>();
        HashMap<String, SamLocusIterator.RecordAndOffset> alt = new HashMap<String, SamLocusIterator.RecordAndOffset>();

        IntervalList region = createTargetRegionInterval(reader, snp);
        SamLocusIterator locus = new SamLocusIterator(reader, region);
        final List<SamRecordFilter> filters = java.util.Arrays.asList(new SecondaryAlignmentFilter());
        // final List<SamRecordFilter> filters = java.util.Arrays.asList(new SecondaryAlignmentFilter(), new DuplicateReadFilter());
        // locus.setSamFilters(filters);
        locus.setIncludeNonPfReads(true);
        locus.setEmitUncoveredLoci(true);
        locus.setIncludeIndels(true);
        locus.setSamFilters(null);
        while(locus.hasNext()){
            AbstractLocusInfo cur = locus.next();
            List<SamLocusIterator.RecordAndOffset> alignments = cur.getRecordAndOffsets();
            // System.err.println("pos: " + cur.getSequenceName() + ":" + cur.getPosition() + " " + (snp.pos == cur.getPosition()));
            // System.err.print("   " + cur.getSequenceName() + " " + cur.getPosition());
            if(cur.getPosition() == snp.pos){
                if(snp.chrom.equals("chr5") && snp.pos > 83520050){
                    debug = true;
                }
                classifyAlignment(alignments, snp.ref, snp.alt, (HashMap<String, SamLocusIterator.RecordAndOffset>) ref, (HashMap<String, SamLocusIterator.RecordAndOffset>) alt);
                if(snp.chrom.equals("chr5") && debug && snp.pos > 83520062){
                    System.exit(0);
                }
                debug = false;
            }
        }
        locus.close();
        return calcScore(ref, alt);
    }
    public static void analyseVCF(ArrayList<VCFEntry> vcf, String bamPath) throws IOException{
        final SamReader reader = SamReaderFactory.makeDefault().open(new File(bamPath));

        VCFEntry prev = null;
        for(VCFEntry snp: vcf){
            System.err.println("checking "+ snp.chrom + ":" + snp.pos);
            String score = analyzeSNP(snp, reader);
            snp.appendInfo(score);
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
