HTSJDK=/data10/system/src/IGV_Linux_2.12.3/lib/htsjdk-2.24.1.jar
JAVA_HOME=/usr/java/jdk-17.0.2


package:
	PATH=$(JAVA_HOME)/bin:$(PATH); mvn compile
	PATH=$(JAVA_HOME)/bin:$(PATH); mvn package

build:
	PATH=$(JAVA_HOME)/bin:$(PATH); mvn compile

fulltest:
	mvn test

install:
	-mkdir /usr/local/lib/chimera-buster
	cp target/chimera_buster-1.0-SNAPSHOT-jar-with-dependencies.jar /usr/local/lib/chimera-buster
	echo 'java -jar /usr/local/lib/chimera-buster/chimera_buster-1.0-SNAPSHOT-jar-with-dependencies.jar $$@' > /usr/local/bin/chimera-buster2
	chmod 755 /usr/local/bin/chimera-buster

test:
	chimera-buster --vcf ../test_data/mini.vcf.gz --bam ../test_data/ML-BA001.Maryland_DES_wave1.hg38.BQSR.merge.markdup.tagfix.mutect.GENCODEv32_RefSeq.slop100.PASS.slop300.bamout.bam --ref ../test_data/Homo_sapiens.GRCh38.dna_sm.primary_assembly_fix.fa
	# echo "check ML-BA001_HMCN1_chr1_185984194.png"
	# @echo "---"
	# java -cp target/chimera_buster-1.0-SNAPSHOT-jar-with-dependencies.jar Main --vcf ../test_data/mini.vcf.gz --bam ../test_data/ML-BA001.Maryland_DES_wave1.hg38.BQSR.merge.markdup.tagfix.mutect.GENCODEv32_RefSeq.slop100.PASS.slop300.bamout.bam \
#		--ref /disk1/2020_digenome_seek_july/Homo_sapiens.GRCh38.dna_sm.primary_assembly_fix.fa

