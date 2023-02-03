# chimera-buster
chimeric read buster

We would like to eliminate false positive SNPs from our 
somatic mutations analysis.

Some of enzymatic fragmentation can bears some fraction of 
artificial short read during its process.
This program can detect such artifact and caluculate some stats
for each SNPs.

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
