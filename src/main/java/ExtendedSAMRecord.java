import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarOperator;
import java.util.HashMap;

public class ExtendedSAMRecord {
    SAMRecord rec = null;
    int start = 0;
    int end = 0;
    byte[] leadingSoftClipBases = null;
    byte[] trailingSoftClipBases = null;
    HashMap<Integer, Byte> alignment = new HashMap<>();

    public ExtendedSAMRecord(SAMRecord r){
        this.rec = r;
        this.start = r.getUnclippedStart();
        this.end = r.getUnclippedEnd();

        Cigar cigar = r.getCigar();
        // leading clip
        if(cigar.getFirstCigarElement().getOperator().isClipping()){
            byte[] bases = this.rec.getReadBases();
            int clip_len = cigar.getFirstCigarElement().getLength();
            leadingSoftClipBases = java.util.Arrays.copyOfRange(bases, 0, clip_len);
        }
        // trailing clip
        if(cigar.getLastCigarElement().getOperator().isClipping()){
            byte[] bases = this.rec.getReadBases();
            int clip_len = cigar.getFirstCigarElement().getLength();
            trailingSoftClipBases = java.util.Arrays.copyOfRange(bases, bases.length - clip_len, bases.length);
        }
    }
    public void init(int pos, byte base){
        this.alignment.put(pos, base);
    }
    public void init(){
        this.alignment = alignment;
        if(leadingSoftClipBases != null){
            for(int i = 0; i<leadingSoftClipBases.length; i++){
                this.alignment.put(this.start + i, leadingSoftClipBases[i]);
            }
        }
        if(trailingSoftClipBases != null){
            for(int i = 0; i<trailingSoftClipBases.length; i++){
                this.alignment.put(this.rec.getEnd() + i, trailingSoftClipBases[i]);
            }
        }
    }
    public SAMRecord get(){
        return rec;
    }
    public int getStart(){
        return start;
    }
    public int getEnd(){
        return end;
    }
    public boolean contains(int pos, int len){
        // len is enough small than read(151bp)
        if(
           (pos >= start && pos <= end) || (pos+len >= start && pos+len<=end)){
            return true;
        }
        return false;
    }
}

