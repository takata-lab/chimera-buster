# chimera-buster
chimeric read buster

We aim to eliminate false positive single nucleotide polymorphisms (SNPs) from our analysis of somatic mutations.

Unfortunately, enzymatic fragmentation can generate artificial short reads during its shearing process.
This program can detect these artifacts and compute corresponding statistics for each SNP.

* Installation (on CentOS or RHEL):

Install Java 17

```
bash install-java.sh
```
- Set paths
```
export HTSJDK=/path/to/htsjdk-2.24.1.jar
export JAVA_HOME=/path/to/java17
export PATH=/usr/local/bin:${JAVA_HOME}/bin:$PATH
```
-Build and install
```
make
sudo make install
```
