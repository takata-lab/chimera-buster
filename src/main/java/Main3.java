import java.io.File;
import java.io.BufferedReader;
import java.io.FileReader;
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
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.filter.SecondaryAlignmentFilter;
import htsjdk.samtools.filter.DuplicateReadFilter;
import htsjdk.samtools.CigarOperator;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;


public class Main3 {
    static IndexedFastaSequenceFile refGenome;
    static boolean debug = false; // 83520057
    static ArrayList<String> header = new ArrayList<>();
    static int unknown = 0;
    static int READ_LENGTH = 151;
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
    private static boolean hasClipping(SAMRecord r){
        Cigar cigar = r.getCigar();
        if(cigar.getFirstCigarElement().getOperator().isClipping() || cigar.getLastCigarElement().getOperator().isClipping()){
            // System.err.println("READNAME: " + r.getReadName() + " true");
            return true;
        }
        return false;
    }
    /*
    private static boolean isSecondaryOrSupplementary(SAMRecord r){
        System.err.println("isSecondaryOrSupplementary()");
        if(r.isSecondaryOrSupplementary()){
            return true;
        }
        return hasClipping(r);
    }*/
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
   public static void analyseVCF(ArrayList<VCFEntry> vcf, String bamPath) throws IOException{
        VCFEntry prev = null;
        for(String line: header){
            if(line.startsWith("#CHROM")){
                printCustomHeader();
            }
            System.out.println(line);
        }
        for(VCFEntry snp: vcf){
            System.err.println("checking "+ snp.chrom + ":" + snp.pos);
            String score = analyzeSNP(snp, bamPath);
            snp.appendInfo(score);
            System.out.println(snp.toString());
        }
    }
    public static String analyzeSNP(VCFEntry snp, String bamPath) throws IOException{
        final SamReader reader1 = SamReaderFactory.makeDefault().open(new File(bamPath));
        final SamReader reader2 = SamReaderFactory.makeDefault().open(new File(bamPath));
        // prepare ovarlapping reads 
        int snp_len = (snp.ref.length() > snp.alt.length())? snp.ref.length(): snp.alt.length();
        // Scan wide and collect reads containing snp position
        HashMap<String, ExtendedSAMRecord> reads = new HashMap<>();
        try (SAMRecordIterator it = reader1.queryOverlapping(snp.chrom, snp.pos-READ_LENGTH, snp.pos +snp_len + READ_LENGTH)){
            while(it.hasNext()){
                SAMRecord read  = it.next();
                if(read.getReadUnmappedFlag()){
                    continue;
                }
                ExtendedSAMRecord eread = null;
                try {
                    eread = new ExtendedSAMRecord(read);           
                    String suffix = read.getFirstOfPairFlag()? "_1": "_2";
                    if(eread.contains(snp.pos, snp_len)){ // check snp position including clipping 
                        reads.put(read.getReadName() + suffix, eread);
                    }
                }catch(RuntimeException hard_clip_ex){
                }
            }
        }
        // scan alignment
        IntervalList region = createWideRegionInterval(reader1, snp);
        try (SamLocusIterator locus = new SamLocusIterator(reader2, region)){
            locus.setIncludeNonPfReads(true); //  no-filter
            locus.setEmitUncoveredLoci(true); //
            locus.setIncludeIndels(true);     // 
            locus.setSamFilters(null);        //
            // setup multiple alignment
            while(locus.hasNext()) {
                SamLocusIterator.LocusInfo cur = (SamLocusIterator.LocusInfo)locus.next();
                // AbstractLocusInfo cur = locus.next();
                List<SamLocusIterator.RecordAndOffset> alignments = cur.getRecordAndOffsets();
                for(SamLocusIterator.RecordAndOffset align: alignments){
                    SAMRecord rec = align.getRecord();
                    String name = rec.getReadName() + (rec.getFirstOfPairFlag()? "_1": "_2");
                    byte[] base = new byte[1];
                    base[0] = align.getReadBase();
                    if(!reads.containsKey(name)){ 
                        try {
                            ExtendedSAMRecord eread = null;
                            eread = new ExtendedSAMRecord(rec);
                            reads.put(name, eread);
                        }catch(RuntimeException e){
                            continue;
                        }
                    }
                    if(!reads.get(name).hasClip()){ // clipped-reads have more priority
                        reads.get(name).init(cur.getPosition(), new String(base));
                    }
                }
                List<SamLocusIterator.RecordAndOffset> deleted = cur.getDeletedInRecord();
                for(SamLocusIterator.RecordAndOffset align: deleted){
                    SAMRecord rec = align.getRecord();
                    String name = rec.getReadName() + (rec.getFirstOfPairFlag()? "_1": "_2");
                    String bases = get_bases(rec, align.getOffset(), align.getLength());
                    if(!reads.containsKey(name)){
                        ExtendedSAMRecord eread = null;
                        try {
                            eread = new ExtendedSAMRecord(rec);
                            reads.put(name, eread);
                        }catch(Exception e){
                            continue;
                        }
                    }
                    if(!reads.get(name).hasClip()){
                        reads.get(name).init(cur.getPosition(), bases);
                    }
                }
                List<SamLocusIterator.RecordAndOffset> inserted = cur.getInsertedInRecord();
                for(SamLocusIterator.RecordAndOffset align: inserted){
                    SAMRecord rec = align.getRecord();
                    String name = rec.getReadName() + (rec.getFirstOfPairFlag()? "_1": "_2");
                    String bases = get_bases(rec, align.getOffset(), align.getLength());
                    if(!reads.containsKey(name)){
                        ExtendedSAMRecord eread = null;
                        try {
                            eread = new ExtendedSAMRecord(rec);
                            reads.put(name, eread);
                        }catch(Exception e){
                            continue;
                        }
                    }
                    if(!reads.get(name).hasClip()){
                        reads.get(name).init(cur.getPosition(), bases);
                    }
                }
            }
        }
        reader2.close();
        reader1.close();
        
        int alt_count = 0;
        int alt_and_chimera = 0;
        int ref_count = 0;
        int ref_and_chimea = 0;
        int alt_allele_length = snp.getAltAllele().length();
        int ref_allele_length = snp.getRefAllele().length();
        for(ExtendedSAMRecord rec: reads.values()){
            int pos = snp.getPosition();
            // char refAllele = getRefBase(snp.getChromosome(), pos);
            if(rec.hasBaseAt(pos + alt_allele_length)){
                String assumedRefAllele = makeAllele(rec, pos, ref_allele_length);
                String assumedAltAllele = makeAllele(rec, pos, alt_allele_length);
                // System.out.println(refAllele + "\t" + allele + "\t" + rec.getReadName() + "\t" + rec.hasClip());
                if(assumedAltAllele != null && snp.getAltAllele().equals(assumedAltAllele)){
                    alt_count++;
                    if(rec.hasClip()){
                        alt_and_chimera++;
                    }
                }else if(assumedRefAllele != null && snp.getRefAllele().equals(assumedRefAllele)){
                    ref_count++;
                    if(rec.hasClip()){
                        ref_and_chimea++;
                    }
                }
            }
        }

        return "Alt=" + alt_count + "/" + (alt_count+ref_count) + ";"
            + "Alt_ratio=" + divide(alt_count, alt_count+ref_count) + ";"
            + "ChimeraAlt=" + alt_and_chimera + "/" + alt_count + ";"
            + "ChimeraAlt_ratio=" + divide(alt_and_chimera, alt_count); 
    }
    private static String makeAllele(ExtendedSAMRecord rec, int pos, int len){
        try {
            StringBuilder allele_buf = new StringBuilder();
            for(int i = pos; i<pos + len; i++){
                allele_buf.append(rec.getBaseAt(i));
            }
            return allele_buf.toString();
        }catch(RuntimeException e){
            return null;
        }
    }
    private static float divide(int a, int b){
        return (float)a/(float)b;
    }
    private static String get_bases(SAMRecord rec, int offset, int len){
        byte[] bases = rec.getReadBases();
        return new String(bases, offset, len);
    }
    private static IntervalList createJustRegionInterval(SamReader reader, VCFEntry snp){
        IntervalList list = new IntervalList(reader.getFileHeader());
        int start = snp.pos;
        int end = snp.pos+snp.alt.length();
        
        Interval interval = new Interval(snp.chrom, start, end);
        list.add(interval);

        return list;
    }
    private static IntervalList createWideRegionInterval(SamReader reader, VCFEntry snp){
        IntervalList list = new IntervalList(reader.getFileHeader());
        int start = snp.pos-10;
        int end = snp.pos+10;
        
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
        System.err.print("testing");
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
        BufferedReader br = null;
        if(vcfPath.endsWith("gz")){
            br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(vcfPath))));
        }else {
            br = new BufferedReader(new FileReader(vcfPath));
        }
        String raw = null;
        while(null != (raw = br.readLine())){
            if(raw.startsWith("#")){
                header.add(raw);
                continue;
            }
            VCFEntry ent = new VCFEntry(raw.split("\t"));
            entries.add(ent);
        }
        br.close();
        return entries;
    }
    public static final void printCustomHeader(){
        System.out.println("##INFO=<ID=Alt,Number=2,Type=Integer,Description=\"Alt carring read count in format 'Alt=A/B' where A is alt read count and B is total read count\">");
        System.out.println("##INFO=<ID=Alt_ratio,Number=1,Type=Float,Description=\"Alt read ratio in format 'Alt_ratio=X' where X is the ratio of alt carring read count to total read count at the point(i.e. depth)\">");
        System.out.println("##INFO=<ID=ChimeraAlt,Number=2,Type=Integer,Description=\"ChimeraAlt read count in format 'ChimeraAlt=A/B' where A is chmeric read count and B is total alt carring read count\">");
        System.out.println("##INFO=<ID=ChimeraAlt_ratio,Number=1,Type=Float,Description=\"ChimeraAlt read ratio in format 'ChimeraAlt_ratio=X' where X is the ratio of chimeric read count to total alt carring read count\">");

    }
    public static final void printUsage(){
        System.err.print("Usage:\n\t\t");
        System.err.println("java -jar ChimeraBuster.jar --vcf [vcf file] --bam [bam file] --ref [ref genome fasta file]");
    }
}
/*
    private static void classifyAlignment(List<SamLocusIterator.RecordAndOffset> src, String refBase, String altBase, HashMap<String, SamLocusIterator.RecordAndOffset> refs, HashMap<String, SamLocusIterator.RecordAndOffset> alts, int snpPos, int pos, HashSet<String> chimeras){
        int total = 0;
        unknown = 0;
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
            // System.err.println("trying classification: " + bases.substring(beginIndex, endIndex));
            try {
            //  if(alignment.getOffset() + (snpPos - pos) >= 0) {
                System.err.println("CHECK: " + snpPos + " " + pos + " " +  bases.substring(beginIndex, endIndex));
            }catch(Exception e){
            }
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
        return;
    }
    public static ChimeraInfo buildInfo(List<SamLocusIterator.RecordAndOffset> src, String refBase, String altBase, HashMap<String, SamLocusIterator.RecordAndOffset> refs, HashMap<String, SamLocusIterator.RecordAndOffset> alts, int snpPos, int pos, HashSet<String> chimeras){
        HashSet<String> alt_chimera = new HashSet<>();
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
         //     System.err.println("---");
         //     System.err.println("secondary: "+ isSecondaryOrSupplementary(rec));
         //     System.err.println("chimera: "+ chimeras.contains(rec.getReadName()));
         //     System.err.println("alts.containsKey: "+ alts.size());
         //     System.err.println("alts.containsKey: "+ alts.containsKey(getReadName(rec)));
            }
        }
        if(true){
            System.err.println("refbase: " + refBase);
            System.err.println("altbase: " + altBase);
            // System.err.println("total: " + total);
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
    }
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
    }
    private static String calcScore(HashMap<String, SamLocusIterator.RecordAndOffset> ref, HashMap<String, SamLocusIterator.RecordAndOffset> alt, ChimeraInfo info){
        int total = info.refs + info.alts + info.unknown;
        float alt_rate = ((float)info.alts)/total;
        String result = join("Alt=", String.valueOf(info.alts), "/", String.valueOf(total), ":", String.valueOf(alt_rate));
        result = join(result, ";ChimeraAlt=", String.valueOf(info.alt_and_chimera), "/", String.valueOf(info.alts), ";ChimeraAlt_ratio=", String.valueOf(((float)info.alt_and_chimera)/(info.alts)));
        result = join(result, ";StrictChimeraAlt=", String.valueOf(info.alt_and_chimera), "/", String.valueOf(info.alts+info.unknown), ";StrictChimeraAlt_ratio=", String.valueOf(((float)info.alt_and_chimera)/(info.alts+info.unknown)));
        return result;
    }

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
        HashSet<String> chimeras = new HashSet<String>();

        IntervalList region = createWideRegionInterval(reader, snp);
        SamLocusIterator locus = new SamLocusIterator(reader, region);
        locus.setIncludeNonPfReads(true);
        locus.setEmitUncoveredLoci(true);
        locus.setIncludeIndels(true);
        locus.setSamFilters(null);
        while(locus.hasNext()){
            SamLocusIterator.LocusInfo cur = (SamLocusIterator.LocusInfo)locus.next();
            // AbstractLocusInfo cur = locus.next();
            List<SamLocusIterator.RecordAndOffset> alignments = cur.getRecordAndOffsets();
            if(cur.getPosition() > snp.pos){
                if(snp.chrom.equals("chr5") && snp.pos > 83520050){
                    debug = true;
                }
                classifyAlignment(alignments, snp.ref, snp.alt, (HashMap<String, SamLocusIterator.RecordAndOffset>) ref, (HashMap<String, SamLocusIterator.RecordAndOffset>) alt, snp.pos, cur.getPosition(), chimeras);
                debug = false;
            }
        }
        locus.close();
        region = createTargetRegionInterval(reader2, snp);
        locus = new SamLocusIterator(reader2, region);
        locus.setIncludeNonPfReads(true);
        locus.setEmitUncoveredLoci(true);
        locus.setIncludeIndels(true);
        locus.setSamFilters(null);
        ChimeraInfo info = null;
        if(locus.hasNext()){
            SamLocusIterator.LocusInfo cur = (SamLocusIterator.LocusInfo)locus.next();
            List<SamLocusIterator.RecordAndOffset> alignments = cur.getRecordAndOffsets();
            info = buildInfo(alignments, snp.ref, snp.alt, (HashMap<String, SamLocusIterator.RecordAndOffset>) ref, (HashMap<String, SamLocusIterator.RecordAndOffset>) alt, snp.pos, cur.getPosition(), chimeras);
        }
        locus.close();
        return calcScore(ref, alt, info);
    }
    private static IntervalList createTargetRegionInterval(SamReader reader, VCFEntry snp){
        IntervalList list = new IntervalList(reader.getFileHeader());
        int start = snp.pos+1;
        int end = snp.pos+1;
        
        Interval interval = new Interval(snp.chrom, start, end);
        list.add(interval);

        return list;
    }
 
 */
