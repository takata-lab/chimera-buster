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
import htsjdk.samtools.Cigar;
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
    public static class ChimeraInfo {
        int alts;
        int refs;
        int chimeras;
        int alt_and_chimera;
        int unknown;
        public ChimeraInfo(int r, int a, int c, int ac, int u){
            alts = a;
            refs = r;
            chimeras = c;
            alt_and_chimera = ac;
            unknown = u;
        }
    }
    private static boolean isSecondaryOrSupplementary(SAMRecord r){
        if(r.isSecondaryOrSupplementary()){
            return true;
        }
        Cigar cigar = r.getCigar();
        if(cigar.getFirstCigarElement().getOperator().isClipping() || cigar.getLastCigarElement().getOperator().isClipping()){
            return true;
        }
        return false;
    }
    private static String getReadName(SAMRecord r){
      try {
        boolean isFirst = r.getFirstOfPairFlag();
        if(isFirst){
            return r.getReadName() + "_1";
        }
      }catch(java.lang.IllegalStateException e){
          System.err.println("uncertain read: " + r.getReadName());
          return r.getReadName() + "_0";
      }
        return r.getReadName() + "_2";
    }
    private static ChimeraInfo classifyAlignment(List<SamLocusIterator.RecordAndOffset> src, String refBase, String altBase, HashMap<String, SamLocusIterator.RecordAndOffset> refs, HashMap<String, SamLocusIterator.RecordAndOffset> alts, int snpPos, int pos, HashSet<String> chimeras){
        int total = 0;
        int unknown = 0;
        HashSet<String> alt_chimera = new HashSet<>();
        for(SamLocusIterator.RecordAndOffset alignment: src){
            total++;
            SAMRecord rec = alignment.getRecord();
            if(isSecondaryOrSupplementary(rec)){
            // if(alignment.getRecord().getSupplementaryAlignmentFlag()){
                chimeras.add(rec.getReadName());
        //        System.err.println("chimera:" + alignment.getRecord().getReadName() + " [" + ((char)alignment.getReadBase()) + "] " + rec.getCigarString() + " " + rec.getFlags());
            }
        //  if(rec.getReadName().equals("A00120:278:HHFTVDSX3:2:2304:11939:16689")){
        //      System.err.println("###A00120:278:HHFTVDSX3:2:2304:11939:16689 [" + ((char)alignment.getReadBase()) + "] sec_supp:" + rec.isSecondaryOrSupplementary() + " read2: " + rec.getSecondOfPairFlag()  + " strand:" + (!rec.getReadNegativeStrandFlag()) + " cigar:" + rec.getCigarString() 
        //          + " offset:" + alignment.getOffset());
        //      System.err.println("    getReferencePositionAtReadPosition: " + rec.getReferencePositionAtReadPosition(22));
        //      System.err.println("    getReferencePositionAtReadPosition: " + rec.getReferencePositionAtReadPosition(21));
        //      System.err.println("    getReadString(): " + rec.getReadString().charAt(alignment.getOffset()-1));
        //      System.err.println("    getReadString(): " + rec.getReadString().charAt(alignment.getOffset()));
        //  }
            String bases = rec.getReadString();
            int beginIndex = alignment.getOffset() + (snpPos - pos);
            int endIndex = beginIndex + altBase.length();
            if(endIndex > bases.length()){
                endIndex = bases.length();
            }
        //  System.err.println("altBase: " + altBase);
        //  System.err.println("refBase: " + altBase);
        //  if(alignment.getOffset() + (snpPos - pos) >= 0){
        //      System.err.println("read: " + bases.substring(beginIndex, endIndex));
        //  }
            if(alignment.getOffset() + (snpPos - pos) >= 0 && bases.substring(beginIndex, endIndex).equals(altBase)){
                refs.remove(getReadName(rec));
                System.err.println("alt");
                alts.put(getReadName(rec), alignment);
            }else if(alignment.getOffset() + (snpPos - pos) >= 0 && bases.substring(beginIndex, endIndex).equals(refBase)){
                alts.remove(getReadName(rec));
                refs.put(getReadName(rec), alignment);
            }else {
                unknown++;
            }
        }
        for(SamLocusIterator.RecordAndOffset alignment: src){
            SAMRecord rec = alignment.getRecord();
            String bases = rec.getReadString();
            int beginIndex = alignment.getOffset() + (snpPos - pos);
            int endIndex = beginIndex + altBase.length();
            if(endIndex > bases.length()){
                endIndex = bases.length();
            }
            if(isSecondaryOrSupplementary(rec) && chimeras.contains(rec.getReadName()) && alts.containsKey(getReadName(rec))){
                System.err.println("#chimera");
                alt_chimera.add(getReadName(rec));
        //  }else if(!isSecondaryOrSupplementary(rec) && alignment.getOffset() + (snpPos - pos) >= 0 && bases.substring(beginIndex, endIndex).equals(altBase)){
        //      if(pos > 52897350 && pos < 52897364){
        //          System.err.println("Non chimeria alt read: " + getReadName(rec));
        //      }
            }else {
                System.err.println("---");
                System.err.println("secondary: "+ isSecondaryOrSupplementary(rec));
                System.err.println("chimera: "+ chimeras.contains(rec.getReadName()));
                System.err.println("alts.containsKey: "+ alts.containsKey(getReadName(rec)));
            }
        }
        if(debug){
            System.err.println("refbase: " + refBase);
            System.err.println("altbase: " + altBase);
            System.err.println("total: " + total);
            System.err.println("ref: " + refs.size());
            System.err.println("alt: " + alts.size());
            System.err.println("chimera: " + chimeras.size());
            System.err.println("alt_and_chimera: " + alt_chimera.size());
        }
        return new ChimeraInfo(refs.size(), alts.size(), chimeras.size(), alt_chimera.size(), unknown);
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
    private static String calcScore(HashMap<String, SamLocusIterator.RecordAndOffset> ref, HashMap<String, SamLocusIterator.RecordAndOffset> alt, ChimeraInfo info){
        int total = info.refs + info.alts + info.unknown;
        float alt_rate = ((float)info.alts)/total;
        String result = join("Alt=", String.valueOf(info.alts), "/", String.valueOf(total), ":", String.valueOf(alt_rate));
        result = join(result, ";ChimeraAlt=", String.valueOf(info.alt_and_chimera), "/", String.valueOf(info.alts), ";ChimeraAlt_ratio=", String.valueOf(((float)info.alt_and_chimera)/(info.alts)));
        result = join(result, ";StrictChimeraAlt=", String.valueOf(info.alt_and_chimera), "/", String.valueOf(info.alts+info.unknown), ";StrictChimeraAlt_ratio=", String.valueOf(((float)info.alt_and_chimera)/(info.alts+info.unknown)));
        return result;
    }
    /*
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
    }*/
    public static String analyzeSNP(VCFEntry snp, SamReader reader){
        // classify ref reads and alt reads
        HashMap<String, SamLocusIterator.RecordAndOffset> ref = new HashMap<String, SamLocusIterator.RecordAndOffset>();
        HashMap<String, SamLocusIterator.RecordAndOffset> alt = new HashMap<String, SamLocusIterator.RecordAndOffset>();
        HashSet<String> chimeras = new HashSet<String>();

        IntervalList region = createTargetRegionInterval(reader, snp);
        SamLocusIterator locus = new SamLocusIterator(reader, region);
        final List<SamRecordFilter> filters = java.util.Arrays.asList(new SecondaryAlignmentFilter());
        // final List<SamRecordFilter> filters = java.util.Arrays.asList(new SecondaryAlignmentFilter(), new DuplicateReadFilter());
        // locus.setSamFilters(filters);
        locus.setIncludeNonPfReads(true);
        locus.setEmitUncoveredLoci(true);
        locus.setIncludeIndels(true);
        locus.setSamFilters(null);
        ChimeraInfo info = null;
        while(locus.hasNext()){
            SamLocusIterator.LocusInfo cur = (SamLocusIterator.LocusInfo)locus.next();
            // AbstractLocusInfo cur = locus.next();
            List<SamLocusIterator.RecordAndOffset> alignments = cur.getRecordAndOffsets();
            if(cur.getPosition() > snp.pos){
                if(snp.chrom.equals("chr5") && snp.pos > 83520050){
                    debug = true;
                }
                info = classifyAlignment(alignments, snp.ref, snp.alt, (HashMap<String, SamLocusIterator.RecordAndOffset>) ref, (HashMap<String, SamLocusIterator.RecordAndOffset>) alt, snp.pos, cur.getPosition(), chimeras);
                /*
                if(snp.chrom.equals("chr5") && debug && snp.pos > 83520062){
                    System.exit(0);
                }*/
                debug = false;
            }
        }
        locus.close();
        return calcScore(ref, alt, info);
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
        int start = snp.pos+1;
        int end = snp.pos+1;
        
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
