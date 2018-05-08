#!/usr/bin/perl

use strict;
use warnings;
use local::lib;
use Parallel::ForkManager;


my $chr1_size 	= 249250621;
#my $chr1_size  = 200000;
my $bin_size 	= $ARGV[0];
my $bamFile 	= $ARGV[1];
my $refGenome 	= $ARGV[2];
my $numProcessors = $ARGV[3];
my $outputPrefix = $ARGV[4];
my $chrom	= $ARGV[5];
my $outputRuntime = $outputPrefix . ".runtime.txt";
my @outputFiles;
my @resultAOA;

my $numberBins = int($chr1_size / $bin_size) + 1;


my $pm = Parallel::ForkManager->new($numProcessors, "/home/kt184/scratch/results/testMultiCore/tempDirParallel/");


$pm->run_on_finish(sub{
    my ($pid,$exit_code,$ident,$exit_signal,$core_dump,$data)=@_;
    my $label = $data->[0];
    my $start = $data->[1];
    my $duration = $data->[2];
    push(@resultAOA, [$label, $start, $duration]);
});


my $progStartTime = `date -u +%s.%N`;
for(my $i=1; $i<=$numberBins; $i++){
	my $start	= 1 + ($i - 1) * $bin_size;
	my $end		= $i * $bin_size;
	my $region 	= $chrom . ":" . $start . "-" . $end;
	my $outputFile 	= $outputPrefix . "_" . $region . ".bcf";
	push(@outputFiles, $outputFile);
	
	$pm->start and next;
	my $startTime = `date -u +%s.%N`;
	system("samtools mpileup -d 100000000000 -Buf $refGenome -r $region $bamFile | bcftools call -O b -c -v > $outputFile");
	my $endTime = `date -u +%s.%N`;
	my $timeElapse = $endTime - $startTime;

	$pm->finish(0, [$region, $start, $timeElapse]);
}

$pm->wait_all_children;
my $parallelEndTime = `date -u +%s.%N`;

my @outputFilesFil;
foreach my $file(@outputFiles){
	# Check if the file is empty. Add to arr if not.
	if(-s $file){
		push(@outputFilesFil, $file);
	}
}

my $fileStr = join(" ", @outputFiles);
my $filesStrFil = join(" ", @outputFilesFil);
my $outputBcf = $outputPrefix . ".combined.bcf";
system("bcftools concat $filesStrFil > $outputBcf");
system("rm $fileStr");

my $progEndTime = `date -u +%s.%N`;
my $parallelDuration = $parallelEndTime - $progStartTime;
my $progDuration = $progEndTime - $progStartTime;
push(@resultAOA, ["parallelTime", 0, $parallelDuration]);
push(@resultAOA, ["progTime", 0, $progDuration]);


open(my $OUTPUT, ">", $outputRuntime) || die $!;
for my $arrRef (@resultAOA){
	my @result = @{ $arrRef};
	print $OUTPUT join("\t", @result) . "\n";
}
close($OUTPUT);
