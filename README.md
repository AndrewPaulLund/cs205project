### Harvard University - CS205 Computing Foundations for Computational Science - Spring 2018 - Final Project

# Genomic Sequencing Analysis Parallelization

Team Members: Andrew Lund, Divyam Misra, Nripsuta Saxena, Kar-Tong Tan

---
### Introduction

**What is genomic sequencing?**

[Genomic sequencing](https://en.wikipedia.org/wiki/DNA_sequencing) is the process
of determining the order of nucleotides in an individual. A human genome has
3 billion nucleotides in the form of four letters (A, C, G, and T).

One of the principal application of genomic sequencing analysis is identifying single
nucleotide polymorphisms (SNP).

[SNPs](https://en.wikipedia.org/wiki/Single-nucleotide_polymorphism) are
differences in single nucleotides that occur at a specific position in
the genome. The below image illustrates SNPs between three individuals.

![SNP](report_images/snp.png)

**What are the uses of genomic sequencing?**

Genomic sequencing has many uses in today's world. Some are as follows:
- Quicker diagnosis of mysterious diseases
- Finding patients with the same disease
- Very important for extremely rare disorders
- Testing for hereditary disorders
- In utero and carrier testing
- Predictive (presymptomatic) testing
- Faster pharmacogenetic testing
- Testing how someone will respond to a certain medication
- Used for certain kinds of cancers

**Explosion of cheap data outweighs compute infrastructure**

Today, genomic data is being generated at a much faster and efficient rate than the
compute infrastructure that supports its analysis. As evidenced in the following
plot from the [National Human Genome Research Institute](https://www.genome.gov/), the cost of genomic sequencing as dramatically decreased in the last decade.

As seen above, around 2008 the cost of sequencing has significantly decreased.
Sequencing the first genome is thought to have cost more than $2.7 billion and
took almost 15 years. Today, it can be as low at $1,000 and take as little as
one week.

The primary overhead for genomic sequencing is now computation. The
algorithms used are not easily parallelized, and do not scale linearly.

<img src="report_images/cost.png" width=500>

---

### Need for parallelization in computational analysis of genomic data

Given the huge amount of genomic data being produced today and the huge number
of individuals being sequenced. For example, both [GenomeAsia100K](http://www.genomeasia100k.com/) and the UK's [National Health
Service](http://www.sciencemag.org/news/2012/12/uk-unveils-plan-sequence-whole-genomes-100000-patients) are trying to sequence 100,000 individual genomes for
medical and population studies. For a single individual the genomic data is
approximately 100-200 GB. The total size of both these project is in the range
of 1-2 PB and is too large to practically manage on a single machine.

Analysis of a single individual's genome can take up to [~10 days](https://www.intel.com/content/dam/www/public/us/en/documents/white-papers/deploying-gatk-best-practices-paper.pdf) on a single threaded process. To analyze 100,000
such samples, it would take 24,000,000 core hours of computing power, or 2700
years on a single core. This is obviously impractical to run on a single core,
hence the need for parallelization.

In clinical applications, this long analysis is too slow. For example, if we could achieve a 100x speedup from parallelization, we reduce the analysis time from 10
days to 1-2 hours, which in a clinical setting could result in a life saved versus lost.

Given the sheer size of the data and computing power required for analysis, we
describe our project as both "Big Data" and "Big Compute."

---

### Need for "efficient" parallelization

Parallelization can result in speedup of analysis, but is often in a non-linear
manner. For example, we may be able to use 10 cores to achieve a 5x speedup. That is why we need to focus on "efficient" parallelization in order to try and achieve a linear speedup and reduce both time and cost. To that end:

1. If the analyses could be completed in-house or on an institutional high-performance computing infrastructure, where the number of nodes is limited and current analysis can take 2-3 months, a more linear speedup of ~20% could result in a time savings of ~2 weeks.

2. If you run analysis on the cloud, these analyses can be quite costly due to
cloud compute pricing. If we used Google Compute cloud infrastructure with a
single core [m1-standard-1](https://cloud.google.com/compute/pricing), running
the 100,000 samples mentioned above, it would cost about $1.14 million. Increasing the
efficiency of the analysis by 10-20% could allow a cost-savings around $100-200K.

The argument for efficient parallelization makes both timely and monetary sense.
To that end, the principle goal of our project is:

**To efficiently speed up SNP analysis through big data and big compute parallelization techniques.**

---

### Model and Data

#### SNP Analysis Model:

The SNP analysis software we used throughout the project is
 [SAMtools](http://www.htslib.org/). It is an open-source suite of programs
 for interacting with high-throughput sequencing data, and can be downloaded at
 the link above. It is a single-threaded program that does not natively support
 parallelization and consists of ~93,000 lines of mostly C code.

 There are three separate repositories that make up SAMtools and are required
 for SNP analysis:
 - [Samtools](https://github.com/samtools/samtools) - used for reading, writing, editing, indexing, and viewing SAM/BAM/CRAM alignment files

 - [BCFtools](https://github.com/samtools/bcftools) - used for reading, writing, BCF2/VCF/gVCF files and calling, filtering, summarizing SNP and short indel
 (insertion plus deletion) sequence variants
 - [HTSlib](https://github.com/samtools/htslib) - a C library used for reading and writing high-throughput sequencing data

Software usage is described below.

#### This Project's Genome Data:
- Two individuals' DNA and RNA alignment and index files - referred to as
"Sample 1", "Sample 2", "DNA1", "DNA2", "RNA1", or "RNA2" throughout this
project.
- Files are public genomes from the
[1000 Genomes Project](http://www.internationalgenome.org/).
- Each alignment file (.bam) is ~10GB, and each index file (.bai) is ~5MB.

Index files are much smaller than alignment files and act like a table of contents to
the alignment file, allowing programs like SAMtools to jump to specific sections
in the alignment file without having to read sequentially through it.

**Alignment & Index Data Downloads:**

DNA1: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/exome_alignment/HG00096.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam

ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/exome_alignment/HG00096.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam.bai

DNA2: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00117/exome_alignment/HG00117.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam

ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00117/exome_alignment/HG00117.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam.bai

RNA1: https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/HG00096.1.M_111124_6.bam

RNA2: https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/HG00117.1.M_111124_2.bam

Index files can be generated using ```$ samtools index sample.bam```

Details about the infrastructure we used are in the "Infrastructure" section
toward the bottom of this report.

---

### Description of solution and comparison with existing work

Our project evaluated four primary parallelization techniques for speeding up SNP
analysis:
1. Binning - distributing DNA & RNA reads (sequence
strings) into “bins” amongst cores
2. OpenMP - shared-memory technique to reduce execution time. This parallelization
technique has potential to speed up analysis time on single nodes with multiple
cores.
3. MPI - distributed-memory technique to reduce execution time across multiple
nodes. This process would spread analysis across independent nodes, and combine
their results into a overall SNP file.
4. Load balancing - reduce analysis heterogeneity by sorting genome alignment
index files. We developed a load balancing simulator described in detail below.

Each of these techniques and their results are described in detail in the
following sections.

A literature source we used for our initial profiling analysis is Weeks and Luecke's 2017
paper, ["Performance Analysis and Optimization of SAMtools Sorting"](papers/samtoolsPaper.pdf).
They used OpenMP to try and optimize alignment file sorting.

---

### Installing, Running, & Profiling SAMtools

(Technical description of the software design, code baseline, dependencies, how to use the code, and system and environment needed to reproduce your tests)

#### Installing the SAMtools analysis suite
SAMtools is already available on the HMSRC cluster. In order to download,
install and run  and  local version for us to modify on the HMSRC login node, we
performed the following steps:

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

All program suite makefiles use the ```-O2``` optimization flag by default.
We did not update this flag throughout the project.

#### SAMtools & BCFtools Dataflow

The flow of data from alignment and index files is illustrated below.

<img src="report_images/dataflow.png" width=500>

SAMtools reads the index and alignment files, pipes standard output to BCFtools,
which identifies SNPs in a data pipeline. The area of parallelization our
project focused on is the ```mpileup``` function within the SAMtools program
suite.

#### Running SAMtools on a sample alignment file
Sample 1 DNA alignment and index files for 10 million reads is found in this
repository's ```data/sample_data``` directory. We used this file extensively
for testing due to it's reasonable serial ```mpileup``` runtimes (~30 seconds).
All our analyses can be easily extended to larger alignment and index files.

The [mpileup](http://www.htslib.org/doc/samtools.html) function was identified
as the primary overhead in the profiling section below. From source code documentation, mpileup "Generate[s] VCF, BCF or pileup for one or multiple BAM files. Alignment records are grouped by sample (SM) identifiers in @RG header lines. If sample identifiers are absent, each input file is regarded as one sample."

In order to run ```mpileup``` with associated output timing follow these steps:

```Bash
# navigate to the local SAMtools installation
$ cd .../samtools
$ ./samtools mpileup
$ time ./samtools mpileup .../HG00096.mapped.ILLUMINA.bwa.GBR.exome.20120522.10mil.bam > /dev/null
```

Running the above command without ```time``` or ```> /dev/null``` will output
the read sequences directly to the terminal window and give you an idea of
what ```mpileup```'s standard output looks like.

#### Running ```mpileup``` batch jobs
In order to automate parallelization speedup analysis of multiple samples with different
parameters, we used [Perl](https://www.perl.org/) and batch scripts on the HMSRC cluster.
Sample batch files are found in the ```binning_tests/batch_scripts``` directory.

Here is a sample from one of the
perl files used for the binning technique described below:

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
#### Profiling
In order to evaluate specific functions within the SAMtools and BCFtools libraries
to focus our parallelization speedup techniques on, we completed profiling.

As seen below, we profiled both SAMtools and BCFtools. Both results are shown,
but we focused more on SAMtools than BCFtools, since SAMtools typically runs more
than 10-times longer than BCFtools for a given sample. Sample timings for
SAMtools and BCFtools on the same sample are found in the
```profiling/runTimes.txt``` file.

#### Profiling SAMtools

Profling was completed using the [gprof](https://sourceware.org/binutils/docs/gprof/)
profiler. In order to get the profiler to run, we had to include both the
```-pg``` flag in both the ```CFLAGS``` and ```LDFLAGS``` line items in the
SAMtools ```Makefile```.

Below is a sample of our SAMtools profiler results. The
full profiler text can be found in the ```profiling``` directory.

```bash
Flat profile:

Each sample counts as 0.01 seconds.
  %   cumulative   self              self     total           
 time   seconds   seconds    calls   s/call   s/call  name    
 51.10      1.16     1.16        1     1.16     2.27  mpileup
 18.06      1.57     0.41  4519010     0.00     0.00  bam_plp_next
 12.11      1.85     0.28 26195998     0.00     0.00  pileup_seq
  9.69      2.07     0.22 26896170     0.00     0.00  resolve_cigar2
  5.73      2.20     0.13  4213250     0.00     0.00  bam_mplp_auto
  0.44      2.21     0.01   353719     0.00     0.00  kh_get_olap_hash

index % time    self  children    called     name
                0.00    2.27       1/1           main [3]
[1]    100.0    0.00    2.27       1         bam_mpileup [1]
                1.16    1.11       1/1           mpileup [2]
                0.00    0.00       1/1           sam_global_args_init [114]
-----------------------------------------------
                1.16    1.11       1/1           bam_mpileup [1]
[2]    100.0    1.16    1.11       1         mpileup [2]
                0.13    0.69 4213250/4213250     bam_mplp_auto [4]
                0.28    0.01 26195998/26195998     pileup_seq [7]
                0.01    0.00 4213249/4213249     mplp_get_ref [20]
                0.00    0.00       1/1           bam_mplp_destroy [21]
                0.00    0.00       1/1           bam_smpl_init [76]
                0.00    0.00       1/1           hts_open_format [106]
                0.00    0.00       1/1           hts_set_opt [107]
                0.00    0.00       1/1           sam_hdr_read [115]
                0.00    0.00       1/1           bam_smpl_add [74]
                0.00    0.00       1/1           bcf_call_add_rg [77]
                0.00    0.00       1/1           bam_mplp_init [69]
                0.00    0.00       1/1           bam_mplp_init_overlaps [70]
                0.00    0.00       1/1           bcf_init [83]
                0.00    0.00       1/1           bam_mplp_set_maxcnt [71]
                0.00    0.00       1/1           bcf_destroy [80]
                0.00    0.00       1/1           bam_smpl_destroy [75]
                0.00    0.00       1/1           bcf_call_del_rghash [78]
                0.00    0.00       1/1           bam_hdr_destroy [65]
                0.00    0.00       1/1           hts_close [103]
-----------------------------------------------
                                                 <spontaneous>
[3]    100.0    0.00    2.27                 main [3]
                0.00    2.27       1/1           bam_mpileup [1]
-----------------------------------------------
                0.13    0.69 4213250/4213250     mpileup [2]
[4]     36.1    0.13    0.69 4213250         bam_mplp_auto [4]
                0.00    0.69 4213250/4213250     bam_plp_auto [5]
-----------------------------------------------
                0.00    0.69 4213250/4213250     bam_mplp_auto [4]
[5]     30.4    0.00    0.69 4213250         bam_plp_auto [5]
                0.41    0.24 4519010/4519010     bam_plp_next [6]
                0.00    0.02  305760/305760      bam_plp_push [9]
                0.01    0.01  305760/305760      mplp_func [10]
-----------------------------------------------
                0.41    0.24 4519010/4519010     bam_plp_auto [5]
[6]     28.6    0.41    0.24 4519010         bam_plp_next [6]
                0.22    0.00 26896170/26896170     resolve_cigar2 [8]
                0.01    0.00  305759/305760      mp_free [15]
                0.00    0.01  305759/305759      overlap_remove [19]
-----------------------------------------------
                0.28    0.01 26195998/26195998     mpileup [2]
[7]     12.6    0.28    0.01 26195998         pileup_seq [7]
                0.01    0.00    2670/2670        printw [18]
-----------------------------------------------
                0.22    0.00 26896170/26896170     bam_plp_next [6]
[8]      9.7    0.22    0.00 26896170         resolve_cigar2 [8]
-----------------------------------------------
```
As seen in the above SAMtools profiling statistics, the primary time overheads
occur in:
1. mpileup
2. bam_mpileup
3. pileup_seq

These three specific functions are the focus of our OpenMP parallelization attempt
outlined below.

#### Profiling BCFtools
To run ```gprof``` on BCFtools, we used the following commands inside the
local BCFtools directory:

```bash
autoheader
autoconf -Wno-syntax
LDFLAGS='-lprofiler -pg' CFLAGS='-lprofiler -pg' ./configure
CPUPROFILE=./prof.out cat ./DNA.10mil_for_bcftools.file | ./bcftools call -O b -c -v > test.bcf
gprof ./bcftools ./gmon.out
```
Below is a sample of our BCFtools profiling. The full output file can be
viewed in the ```profiling``` directory.

```bash
Flat profile:

Each sample counts as 0.01 seconds.
  %   cumulative   self              self     total           
 time   seconds   seconds    calls  ms/call  ms/call  name    
 23.08      0.21     0.21                             bcf_dec_size_safe
 15.38      0.35     0.14                             bcf_fmt_sized_array
 12.09      0.46     0.11                             bcf_fmt_array
  8.79      0.54     0.08        1    80.00    90.00  main_vcfcall
  7.69      0.61     0.07                             bgzf_read
  6.59      0.67     0.06                             _reader_next_line
  5.49      0.72     0.05                             bcf_dec_typed_int1_safe
  5.49      0.77     0.05                             bcf_sr_sort_next
  5.49      0.82     0.05                             kputc
  4.40      0.86     0.04                             bcf_clear
  2.20      0.88     0.02                             bcf_sr_next_line
  1.10      0.89     0.01    20067     0.00     0.00  mc_cal_afs
  1.10      0.90     0.01                             bcf_read
  1.10      0.91     0.01                             bcf_unpack

granularity: each sample hit covers 2 byte(s) for 1.10% of 0.91 seconds

index % time    self  children    called     name
                                                 <spontaneous>
[1]     23.1    0.21    0.00                 bcf_dec_size_safe [1]
-----------------------------------------------
                                                 <spontaneous>
[2]     15.4    0.14    0.00                 bcf_fmt_sized_array [2]
-----------------------------------------------
                                                 <spontaneous>
[3]     12.1    0.11    0.00                 bcf_fmt_array [3]
-----------------------------------------------
                0.08    0.01       1/1           main [5]
[4]      9.9    0.08    0.01       1         main_vcfcall [4]
                0.00    0.01   20067/20067       ccall [14]
                0.00    0.00   20067/20067       set_ploidy [32]
                0.00    0.00       1/1           ploidy_init_string [78]
                0.00    0.00       1/1           init_data [65]
                0.00    0.00       1/1           destroy_data [63]
-----------------------------------------------
                                                 <spontaneous>
[5]      9.9    0.00    0.09                 main [5]
                0.08    0.01       1/1           main_vcfcall [4]
-----------------------------------------------
                                                 <spontaneous>
[6]      7.7    0.07    0.00                 bgzf_read [6]
-----------------------------------------------

```

---

### Speedup Techniques

(Technical description of the parallel application and programming models used)

#### Data Heterogeneity
A key attribute of the data, and all genomic alignment files is that it is
extremely heterogeneous. Each chunk of data can have orders of magnitude different
number of reads (genome sequence strings). This makes sequentially processing
the alignment files very uneven. This unpredictable sizing is illustrated well for
both DNA and RNA of Sample 1 below:

|  Index File Distribution  | Heterogeneity |
|:---:|:---:|
|![DNA1](report_images/DNA_Distribution.png)  |  ![DNA2](report_images/DNA_Heterogeneity.png)|
|![RNA1](report_images/RNA_Distribution.png)  |  ![RNA2](report_images/RNA_Heterogeneity.png)|


#### Parallelization Techniques

**1. Binning** - Our first method of speeding up the SNP analysis involved what
we termed "binning" or distributing the DNA and RNA reads (sequence strings)
into bins amongst cores of a CPU. We distributed them into sequential chunks
from $10^4$ to $10^9$. This was done initially with one chromosome on one core
to determine the ideal bin size. After determining the ideal bin size, we tested
that bin size across a number of cores from 2 to 20. A sample batch script for
binning is outlined in the "Running ```mpileup``` batch jobs" section above, and
Results are discussed in the "Results" section below.

**2. OpenMP** - Based on the profiling outlined above, we focused our OpenMP
parallelization on the mpileup, bam_mpileup, pileup_seq functions. All three of
these functions are found in the SAMtools library, specifically in the
```bam_plcmd.c``` file.

It is important to note that when recompiling the SAMtools library, one must
include the ```-fopenmp``` flag on both the ```CFLAGS``` and ```LDFLAGS```
lines within the ```Makefile```. It is also important to add
```#include omp.h``` in the module header.

We moved sequentially through the functions, adding
```#pragma``` lines on for loops, recompiling and making using
```$ make clean``` and ```$ make all``` commands, and running SAMtools on the
HG00096.mapped.ILLUMINA.bwa.GBR.exome.20120522.10mil.bam sample found in the
```data/sample_data``` directory.

To run a timed OpenMP test on our small chromosome sample we used the following
generic command:

```Bash
$ time ./samtools mpileup sample.bam > dev/null
```
The number of OpenMP threads was adjusted using the ```$ export OMP_NUM_THREADS="n"```
command.

Results for our OpenMP tests are found in the
"Results" section below.


**3. MPI** - FILL IN HERE!

We performed the MPI runs on the Harvard Medical School cluster. As the cluster is
a shared compute cluster used by thousand of users, we do not have root access to
install the packages into the default system file paths for MPI. As such, we built
a custom environment based on Anaconda.

The anaconda package that was pre-installed on the cluster was loaded using the following
command:

```Bash
module load conda2/4.2.13
```

We then installed pip, a custom python package manager as follows:

```Bash
conda install pip
```

A custom python environment and the mpi4py library was then installed using
the following command:

```
pip install virtualenv
virtualenv cs205
source cs205 activate
pip install mpi4py
```

As we had compatability issues using the openmpi package provided
on the cluster with the mpi4py package that we had installed. We
manually installed openmpi-3.0.1 into a custom local path via the
following command aftere downloading the package into a local directory.

```Bash
./configure --prefix=$PWD/cs205_mpi/
make
make install
```

The version of mpi in the current conda environment was then utilized
for running mpi rather than the default openmpi installed on the
cluster
```Bash
/home/kt184/.conda/envs/cs205/bin/mpirun
```

Having set up openmpi with mpi4py in our environment, we then wrote a
job script and then submitted it onto the cluster for the run. The
MPI job was sent to an 'MPI' queue that was set up on the cluster
and for which we had to request special permission to use.

A sample jobscript for running the mpi python script we had written
for the project is as follows:

```Bash
#!/bin/bash
#SBATCH -p mpi #partition
#SBATCH -t 0-24:01 #time days-hr:min
#SBATCH -n 3 #number of tasks
#SBATCH --mem=8G #memory per job (all cores), GB
#SBATCH -o %j.out #out file
#SBATCH -e %j.err #error file

source activate cs205

extension=/home/kt184/cs205/mpi_examples/mpi4py-examples/results/DNA1.core2
/home/kt184/.conda/envs/cs205/bin/mpirun -n 3 python ./runMPIpileup.py 1000000 /n/scratch2/kt184/data/DNA/HG00096.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam /home/kt184/scratch/data/genome/KT_package/icgc_genome_broadvariant/Homo_sapiens_assembly19.fasta 3 $extension 1 249250621
```

In the above jobscript, we requested 3 tasks to represent a single master node and 2 worker nodes. The master node would decide upon the job needed to be performed by the slave node, collect back the results of analysis by the slave nodes and compile the results.

The sample python script we used for MPI4py can be found in the link ```mpi_tests/code```


The `runMPIpileup.py` python code starts a master node and allow us to vary the number of worker nodes available by MPI. Results are compiled
by the master node after all the workers have finished their jobs. We then tested how the speedup for the DNA and RNA files varies as the number of cores (processes) available for MPI varies.



**4. Load Balancing** - We often see a non-linear speed-up due to the data's
heterogeneity as outlined above.

For every alignment (.bam) file, an index (.bai) file can be generated. This index file does not contain any sequence data; it essentially acts like a table of contents for the alignment file. It is typically used to directly jump to specific parts of the alignment file without needing to read the entire fie sequentially, which can be incredibly helpful since an alignment file is typically quite big (for example, our alignment files are ~10GB). 

Our plan was to read the index file by converting it into text format, so that we could analyze it to try and determine the optimal way to divide the data during the analysis computation. 

An alignment (.bam) file can be easily converted from binary to human-readable text format (.sam, which stands for Sequence Alignment Map), as viewing the sequence in a human-readble form has many uses. However, since index files are mostly used only to be able to jump to the relavant sections of the alignment file, converting an index file into a human-readable text format is not a straightforward task. We realized just how difficult a task this was once we started working on this project. We consistently kept running into roadblocks while trying to read the index file. Multiple approaches were implemented; unfortunately none of them worked out the way we had hoped. 

Our first approach was to understand where the samtools code was reading the index file. We tried to append the code so that it would read the index file and try to convert it to a text format. That did not work. We also tried to use an open-source library for samtools written in Java (htsjdk) so that we could read the metadata from the index files. That did not work either. Lack of domain knowledge about genomic sequening data hindered us from writing our own code to convert the index file.

After we realized that we really needed to move on and focus on other aspects of the project, we decided to try different techniques to balance the load and examine which technique did better.


In order to process the heterogeneous data we
developed a load balancing simulator. The simulator
(```simulateLoadBalance.py```), batch script, sample input
and output timing files are found in the ```load_balance_simulator``` directory. It
simulates four different sorting
techniques to parallelize the data across a range of cores.

1. Ascending data size processing
2. Original data order processing
3. Random order processing
4. Descending data size processing

Results for these four sorting techniques are discussed in the "Results" section
below.

**Running Load Balancing Simulations**

To run the simulator, use a command similar to the following:
```Bash
$ pypy simulateLoadBalance.py input_sample.txt > output.simulate.txt
```


---

### Results

(Performance evaluation (speed-up, throughput, weak and strong scaling) and
discussion about overheads and optimizations done)

**1. Binning**

There are two principal steps to our binning parallelization technique.
1. Determine the minimum bin size for sequence reads.
2. Spread bins across cores from 2-20

As seen below, we determined that at least 1,000,000 reads per bin is the optimal
number. For both DNA and RNA, the one chromosome analysis runtime decreases
precipitously from 10,000 to 1,000,000 then levels off. We also decided to use
1,000,000 as the optimal size as increasing the bin size further will reduce the
total number of bins for a sample which reduces the benefit of parallelization.

As seen in the time and speedup plots, we see excellent speedup for our binning
technique and are pleased with the results.

|  DNA  | RNA |
|:---:|:---:|
|![bin1](report_images/binning_dna_bins.png)  |  ![bin3](report_images/binning_rna_bins.png)|
|![bin2](report_images/binning_dna_time.png)  |  ![bin4](report_images/binning_rna_time.png)|
|![bin5](report_images/binning_dna_speedup.png)  |  ![bin6](report_images/binning_rna_speedup.png)|
|![bin7](report_images/dna_binning_table.png) | ![bin8](report_images/rna_binning_table.png) |

**2. OpenMP**
As previously noted, we focused our OpenMP parallelization on three functions
within the ```bam_plcmd.c``` module of SAMtools. There were no for loops that
were parallizable in the bam_mpileup function. Module files for each OpenMP
attempt are in the ```open_mp_tests``` directory. The results for our
10 million read sample follow:

| Function | 1 thread time (seconds) | 8 threads time (seconds)
|---|---|---|
| no pragmas | 22.3 | --- |  
| mpileup (8 pragmas) | 37.5 | 37.4 |
| pileup_seq (2 pragmas) | 103.7 | 50.5 |

As seen above, all attempts to speed up the ```mpileup``` module using OpenMP
resulted in slower execution times. We suspect this is due to complex
interdependencies between SAMtools and HTSlib that we did not have time to
investigate within the scope of this project.

We were also ultimately unable to properly compile BCFtools with the ```-fopenmp```
flag on the HMSRC cluster, and due to
its relatively slow runtime when compared to SAMtools we did not further attempt
to parallelize BCFtools with OpenMP.

**3. MPI**

MPI results follow:

|Execution Time|Speedup|
|:--:|:--:|
|![](report_images/mpi_time.png)|![](report_images/mpi_speedup.png)
||![](report_images/mpi_table.png)

**4. Load Balancing**

|  DNA  | RNA |
|:---:|:---:|
|![sim1](report_images/dna_sim_idle.png)  |  ![](report_images/rna_sim_idle.png)|
|![sim2](report_images/dna_sim_speedup.png)  |  ![](report_images/rna_sim_speedup.png)|
|![](report_images/dna_sim_table.png)|![](report_images/rna_sim_table.png)

---

### Description of advanced features like models/platforms not explained in class, advanced functions of modules, techniques to mitigate overheads, challenging parallelization or implementation aspects...
- Our primary advanced features we developed is the load balancing simulation.
This module can be used to simulate balancing blah blah blah.
- We initially had a very hard timing profiling the SAMtools library. We spent more
than two weeks attempting to run the gprof profiler, and finally had a
breakthrough when we learned we needed to include the ```-pg``` flag in both the
CGLAFS and LDFLAGS sections of the associated ```Makefile```. This same technique
was used to compile our OpenMP parallelization attempts.
- OpenMP was not trivial, and unsuccessful in speeding up execution time

Use of load balancing is an advanced feature in my opinion.

---

### Discussion about goals achieved, improvements suggested, lessons learnt, future work, interesting insights…

We are very happy with the speedup we achieved through binning and simulated load balancing. After determining
the optimal bin size of 1 million
- Discuss trying to read the index file with JAVA, etc.
- We would have liked to include Spark. ADD COMMENTARY HERE

---

### Infrastructure
<img src="report_images/hms.png" width=100>

Our team used the [Harvard Medical School Research Computing](https://rc.hms.harvard.edu/) (HMSRC) cluster for all testing and analysis.

**HMSRC description:**
- 8,000 cores
  - 32 or 28 cores per node
    - Nodes support up to 32 cores, but are capped at 20
    - A known problem in parallelization of related genomic algorithms
  - 256 GB RAM per node (8-9 GB/core)
- 24 GPUs (8 M40 / 16 K80) - Not utilized for this project
- Login/load balancer 5 VM (8 cores/16 GB memory)
- InfiniBandconnectivity between nodes
- CentOS 7
- SLURM batch system
- Supports OpenMP, MPI
- Currently does not support Spark (we initially intended on using Spark for
  load balancing implementation)

**HMSRC Cluster Architecture:**

<img src="report_images/cluster.png" width=500>

---

### References

Throughout this report sources are cited through inline links.
They appear in the following order:
- ADD SOURCES HERE

---

**This README is our project website and serves as its final report.**

- Samtools source code is found [here](https://github.com/samtools).
- Evaluation data, modified source code, simulation and batch scripts, and visualization
notebooks are all found in this repository

# TODO:
- Add makefiles for gprof and OpenMP
- Add OpenMP trial files
- Add R script for plots
- Add MPI work
- Add load balancing work
- mpileup vs. bcftools runtimes (a lot more time for mpileup)
