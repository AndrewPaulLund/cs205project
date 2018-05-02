# Genomic Sequencing Analysis Parallelization

Team Members: Andrew Lund, Divyam Misra, Nripsuta Saxena, Kar-Tong Tan

## Project Website
**This README is our project website and serves as the final report for
the work we did throughout the semester.**

- Samtools source code is found [here](https://github.com/samtools)
- Evaluation data, modified source code, batch scripts, and visualization
notebooks are found in this repository

## Add makefiles for gprof and OpenMP
## Add OpenMP trial files
## Add R script for plots
## Add MPI work
## Add load balancing work

---

### Description of problem and the need for HPC and/or Big Data

[Genomic sequencing](https://en.wikipedia.org/wiki/DNA_sequencing) has many uses
in today's world. Some are as follows:
- Quicker diagnosis of mysterious diseases
- Finding patients with the same disease
- Very important for extremely rare disorders
- Testing for hereditary disorders
- In utero and carrier testing
- Predictive (presymptomatic) testing
- Faster pharmacogenetic testing
- Testing how someone will respond to a certain medication
- Used for certain kinds of cancers

The cost of genomic sequencing as dramatically decreased in the last decade as
evidenced in the following plot from the [National Human Genome Research Institute](https://www.genome.gov/).

<img src="report_images/cost.png" width=500>

As seen above, around 2008 the cost of sequencing has significantly decreased.
Sequencing the first genome is thought to have cost more than $2.7 billion and
took almost 15 years. Today, it can be as low at $1,000 and take as little as
one week.

The primary overhead for genomic sequencing is now computation. The
algorithms used are not easily parallelized, and do not scale linearly.

It is important for sequencing results to be returned in a timely manner
to clinicians and researchers for clinical applications. Greater than one
week runtimes in some applications is too slow for timely diagnostics.

One of the principal application of genomic analysis is identifying single
nucleotide polymorphisms (SNP).

[SNPs](https://en.wikipedia.org/wiki/Single-nucleotide_polymorphism) are
differences in single nucleotides that occur at a specific position in
the genome. The below image illustrates SNPs between three individuals.

![SNP](report_images/snp.png)

We describe our project as both HTC (high-throughput computing) and Big Data.
- HTC - in that there are almost a billion of high-frequency sequence reads for a
person's genome with today's techniques.
- Big Data - in that the genomic alignment files can exceed 200GB.

To that end, the principle goal of our project is:

**To speed up the identification of single nucleotide polymorphisms (SNP) in DNA
and RNA through big data and high throughput computing parallelization
 techniques.**

---

### Description of solution and comparison with existing work

Our project evaluated four primary speedup techniques for speeding up SNP
analysis:
1. Binning - distributing DNA & RNA reads (sequence
strings) into “bins” amongst cores
2. OpenMP - shared-memory technique to reduce execution time
3. MPI - distributed-memory technique to reduce execution timely
4. Load balancing - reduce heterogeneity by sorting genome alignment index files

Each of these techniques and their results are described in detail in the
following sections.

A source we used for our initial profiling analysis is Weeks and Luecke's 2017
paper, ["Performance Analysis and Optimization of SAMtools Sorting"](papers/samtoolsPaper.pdf). They used OpenMP to try and optimize
alignment file sorting.

---

### Model and Data

#### Model:

The SNP analysis software we use throughout the project is
 [SAMtools](http://www.htslib.org/). It is an open-source suite of programs
 for interacting with high-throughput sequencing data, and can be downloaded at
 the link above.

 There are three separate reposititories that make up SAMtools and are required
 for SNP analysis:
 - [Samtools](https://github.com/samtools/samtools) - used for reading, writing, editing, indexing, and viewing SAM/BAM/CRAM alignment files
 - [BCFtools](https://github.com/samtools/bcftools) - used for reading, writing, BCF2/VCF/gVCF files and calling, filtering, summarising SNP and short indel
 (insertion plus deletion) sequence variants
 - [HTSlib](https://github.com/samtools/htslib) - a C library used for reading and writing high-throughput sequencing data


#### Data:
- Two individuals' DNA and RNA alignment and index files: Sample 1 and Sample 2
- Files are public genomes from
the [1000 Genomes Project](http://www.internationalgenome.org/).
- Each alignment file (.bam) is about 10GB, and each index file (.bai) is about 5MB.

**Downloads:**

DNA1:

ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/exome_alignment/HG00096.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/exome_alignment/HG00096.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam.bai

DNA2:

ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00117/exome_alignment/HG00117.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00117/exome_alignment/HG00117.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam.bai

RNA1:


RNA2:


A key attribute of the data, and all genomic alignment files is that it is
extremely heterogeneous. Each chunk of data can have orders of magnitude different
number of reads (genome sequence strings). This makes sequentially processing
the alignment files very uneven. This unpredictable sizing is illustrated well for
both DNA and RNA of Sample 1 below:

|  Index File Distribution  | Heterogeneity |
|:---:|:---:|
|![DNA1](report_images/DNA_Distribution.png)  |  ![DNA2](report_images/DNA_Heterogeneity.png)|
|![RNA1](report_images/RNA_Distribution.png)  |  ![RNA2](report_images/RNA_Heterogeneity.png)|

---

### Infrastructure
<img src="report_images/hms.png" width=100>

Our team used the [Harvard Medical School Research Computing](https://rc.hms.harvard.edu/) (HMSRC) cluster for all our testing and analysis.

HMSRC description:

# ADD MORE HERE

- 8,000 cores with several PB network storage
- Nodes support up to 32 cores, but are capped at 30
- This is a known problem in parallelization of related algorithms

---

### Installing, Running, & Profiling SAMtools

(Technical description of the software design, code baseline, dependencies, how to use the code, and system and environment needed to reproduce your tests)

#### Installing the SAMtools suite
SAMtools is already available on the HMSRC cluster. In order to download,
install and run  and  local version on the login node, we performed the following steps:

```Bash
# clone repositories, and install programs
# HTSlib
$ git clone https://github.com/samtools/htslib.git
$ cd htslib
$ autoheader
$ autoconf
$ ./configure
$ make
$ make install
# SAMtools
$ git clone https://github.com/samtools/samtools.git
$ cd samtools
$ autoheader
$ autoconf
$ ./configure
$ make
$ make install
# BCFtools
$ git clone https://github.com/samtools/bcftools.git
$ cd bcftools
$ autoheader
$ autoconf
$ ./configure
$ make
$ make install
```

#### Running SAMtools on a sample alignment file
A sample alignment and index file for 10 million reads is found in this
repository's ```data/sample_data``` directory.

[mpileup](http://www.htslib.org/doc/samtools.html) is the function we identified
as the primary overhead in the profiling section. From the documentation, mpileup "Generate[s] VCF, BCF or pileup for one or multiple BAM files. Alignment records are grouped by sample (SM) identifiers in @RG header lines. If sample identifiers are absent, each input file is regarded as one sample."

In order to run ```mpileup``` with associated output timing follow these steps:

```Bash
# navigate to the local samtools installation
$ cd .../samtools
$ ./samtools mpileup
$ time ./samtools mpileup .../HG00096.mapped.ILLUMINA.bwa.GBR.exome.20120522.10mil.bam > /dev/null
```

Running the above command without ```time``` or ```> /dev/null``` will output
the read sequences directly to the terminal window.

#### Running ```mpileup``` batch jobs
In order to automate speedup analysis of multiple samples with different
parameters, we used perl and batch scripts on the HMSRC cluster.
Sample files are found in the ```data/batch_scripts```

Here is a sample from one of the
perl files used for binning:

```bash
#!/usr/bin/perl

use strict;
use warnings;

my @cores = (1,2,4,8,16,20);
my @binSize = (10000,100000,1000000,10000000,100000000,1000000000);
my $bamFile = "/n/scratch2/kt184/data/DNA/HG00117.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam";
my $refGenome = "/n/scratch2/kt184/genome/dna/hs37d5.fa";
my $chrom = "1";
my $outputJobScriptPrefix = "/home/kt184/scratch/results/testMultiCore/generateScripts/analysisRuntime/sampleDNA2";

for my $core(@cores){
for my $bin(@binSize){

my $jobScriptFile = $outputJobScriptPrefix . ".cores" . $core . ".binSize" . $bin . ".sh";
my $outputDirectory = $outputJobScriptPrefix . ".cores" . $core . ".binSize" . $bin . "/";
my $outputPrefix = $outputDirectory . "runDNA2." . ".cores" . $core . ".binSize" . $bin;
system("mkdir $outputDirectory");

open(my $OUTPUT, ">", $jobScriptFile) || die $!;

# Lazy to tab in
print $OUTPUT "#!/bin/bash\n";
print $OUTPUT "#SBATCH -p short #partition\n";
print $OUTPUT "#SBATCH -t 0-12:00 #time days-hr:min\n";
print $OUTPUT "#SBATCH -c $core #number of cores\n";
print $OUTPUT "#SBATCH -N 1 #confine cores to 1 node, default\n";
print $OUTPUT "#SBATCH --mem=8G #memory per job (all cores), GB\n";
print $OUTPUT "#SBATCH -o %j.out #out file\n";
print $OUTPUT "#SBATCH -e %j.err #error file\n";
print $OUTPUT "\n";
print $OUTPUT "\n";
print $OUTPUT "module load gcc/6.2.0\n";
print $OUTPUT "module load samtools/1.3.1\n";
print $OUTPUT "module load bcftools\n";
print $OUTPUT "module load perl/5.24.0\n";
print $OUTPUT "eval \$(perl -I\$HOME/perl5/lib/perl5 -Mlocal::lib)\n";


print $OUTPUT "bin_size=$bin\n";
print $OUTPUT "bamFile=\"$bamFile\"\n";
print $OUTPUT "refGenome=\"$refGenome\"\n";
print $OUTPUT "numProcessors=$core\n";
print $OUTPUT "outputPrefix=$outputPrefix\n";
print $OUTPUT "chrom=\"$chrom\"\n";


print $OUTPUT "/n/scratch2/kt184/results/testMultiCore/runParallelJobScript.pl \$bin_size \$bamFile \$refGenome \$numProcessors \$outputPrefix \$chrom\n";

close($OUTPUT);

print "sbatch " .  $jobScriptFile . "\n";
}
}
```

To Profile....


To run our OpenMP jobs...

To run our MPI jobs...

To run load balancing jobs...

---

### Speedup Techniques

(Technical description of the parallel application and programming models used)

#### Explain each technique in detail here

1. Binning
2. OpenMP
3. MPI
4. Load Balancing

---

### Results

(Performance evaluation (speed-up, throughput, weak and strong scaling) and discussion about overheads and optimizations done)

**1. Binning**

|  DNA  | RNA |
|:---:|:---:|
|![bin1](report_images/binning1.png)  |  ![bin3](report_images/binning2.png)|
|![bin2](report_images/binning3.png)  |  ![bin4](report_images/binning4.png)|


---

### Description of advanced features like models/platforms not explained in class, advanced functions of modules, techniques to mitigate overheads, challenging parallelization or implementation aspects...

We initially had a very hard timing profiling the SAMtools library. We spent more
than two weeks trying to get the gprof profiler to run, and finally had a
breakthrough when we realized we needed to include

Use of load balancing is an advanced feature in my opinion.

---

### Discussion about goals achieved, improvements suggested, lessons learnt, future work, interesting insights…

We are very happy with the speedup we achieved through binning.

---

### References

---
