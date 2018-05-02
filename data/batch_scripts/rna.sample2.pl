#!/usr/bin/perl

use strict;
use warnings;

my @cores = (1,2,4,8,16,20);
my @binSize = (10000,100000,1000000,10000000,100000000,1000000000);
my $bamFile = "/n/scratch2/kt184/data/RNA/HG00117.1.M_111124_2.bam";
my $refGenome = "/n/scratch2/kt184/data/genome/KT_package/hg19_kt/hg19.fa";
my $chrom = "chr1";
my $outputJobScriptPrefix = "/home/kt184/scratch/results/testMultiCore/generateScripts/analysisRuntime/RNA/sampleRNA2";

for my $core(@cores){
for my $bin(@binSize){

my $jobScriptFile = $outputJobScriptPrefix . ".cores" . $core . ".binSize" . $bin . ".sh";
my $outputDirectory = $outputJobScriptPrefix . ".cores" . $core . ".binSize" . $bin . "/";
my $outputPrefix = $outputDirectory . "runRNA2" . ".cores" . $core . ".binSize" . $bin;
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
