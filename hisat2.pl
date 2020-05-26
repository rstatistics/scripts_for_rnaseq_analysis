#!/usr/bin/env perl
use warnings;
use strict;
use 5.010;
# use Data::Printer;
use File::Basename;
use IPC::Cmd qw[can_run run];
use File::Copy;
use Getopt::Long;

# Author: Wei Zhang
# Email: admin@ncrna.net
# Date: 2020-05-27

my $usage = <<USAGE;

SYSNOPSIS

This program is used to batch run hisat2.pl

----------------------------------------------------------
----------------------------------------------------------

hisat2.pl [options] hisat2-index file|glob

eg. hisat2.pl /db/Human/hisat2-index/Homo_sapiens.GRCh38.p10 *.fastq.gz

 Options:
   -r|--reads         'pe' or 'se' end reads, default 'pe'
   -s|--strandness    library strandness, default '--fr'
   -a|--args          argments used for hisat2
   -o|--output        output folder for removed reads 
   -p|--progress      progress, default 16

NOTE: the version of samtools must be 1.9 or latter.

USAGE
my $reads      = 'pe';
my $out_folder = dirname './';
my $strandness = "--fr";
my $progress   = 16;
my $args_file  = '';
die $usage
  unless GetOptions(
    "r|reads:s"     => \$reads,
    "a|args:s"      => \$args_file,
    "s|strandness"  => \$strandness,
    "o|output:s"    => \$out_folder,
    "p|progress=i"  => \$progress,
  );
##########################################################################################
#Checking the parameter infomation
##########################################################################################  
#check the genome database
my $hisat2_index = shift;
die "\nERROR: Does not specifed hisat2 index database!\n$usage" unless $hisat2_index;
$hisat2_index =~ s/\/$//;
-s "$hisat2_index.1.ht2"
  || die "\nERROR: Cannot detect hisat2 index file: $hisat2_index!\n$usage";
my @fqs;
#check the fastq files;
foreach my $glob_str (@ARGV) {
	push @fqs, grep {-s $_} glob $glob_str;
}
die $usage unless @fqs;

#check the running environment
can_run('hisat2') or die 'hisat2 is not installed!';
can_run('samtools') or die 'samtools is not installed!';

#mkdir 
$out_folder =~ s/[\/|\|]+$//;
mkdir $out_folder unless -d $out_folder;
my $log_folder = $out_folder."/logs";
mkdir $log_folder unless -d $log_folder;
my $summary_folder = $out_folder."/summary";
mkdir $summary_folder unless -d $summary_folder;
my $detail_folder = $out_folder."/details";
mkdir $detail_folder unless -d $detail_folder;


##########################################################################################
#Getting the mapping infomation
##########################################################################################
my %fastqs = ();
my @suffx  = qw (_1.fq.gz _2.fq.gz _1.fastq.gz _2.fastq.gz _1.fq _2.fq _1.fastq _2.fastq
  _R1.fq.gz _R2.fq.gz _R1.fastq.gz _R2.fastq.gz _R1.fq _R2.fq _R1.fastq _R2.fastq
  .1.fq.gz .2.fq.gz .1.fastq.gz .2.fastq.gz .1.fq .2.fq .1.fastq .2.fastq
  .R1.fq.gz .R2.fq.gz .R1.fastq.gz .R2.fastq.gz .R1.fq .R2.fq .R1.fastq .R2.fastq
  .fq.gz .fq .fastq.gz .fastq);
if ( $reads eq 'pe' ) {
    foreach my $fq ( sort @fqs ) {
        if ( $fq =~ /(.*)[\.|_|\.R|_R]([1|2])\.(fq\.gz|fq|fastq\.gz|fastq)/ ) {
            my $basename = basename( $fq, @suffx );
            my $mate = $2;
            $fastqs{$basename}{$mate} = $fq;
        }
        else {
            warn "Not supported format:$fq\n";
        }
    }
}
else {
    foreach my $fq ( sort @fqs ) {
        if ( $fq =~ /(.*)\.(fq\.gz|fq|fastq\.gz|fastq)/ ) {
            my $basename = basename( $fq, @suffx );
            $fastqs{$basename} = $fq;
        }
        else {
            warn "Not supported format:$fq\n";
        }
    }
}

#######################################################################
#construct general hisat2 papameter
#######################################################################
my $hisat2_param = '';
$hisat2_param .= "$strandness " if defined $strandness;
$hisat2_param .= "--threads $progress ";

if ($args_file) {
	open my $args_fh, "<", $args_file or die "cannot open $args_file\n";
	while (<$args_fh>) {
		chomp;
		$hisat2_param .= "$_ ";		
	}
}

my %fq_parameters = ();
foreach my $basename ( sort keys %fastqs ) {
	my $fq_in_param          = "";
    if ( $reads eq "pe" ) {
        $fq_in_param  .= " -1 $fastqs{$basename}{1} -2 $fastqs{$basename}{2} ";
    }
    else {
        $fq_in_param  .= " -U $fastqs{$basename} ";
    }
	$fq_parameters{$basename} = $fq_in_param;
}

#######################################################################
#mkdir log files and output folders
#######################################################################
my $hisat2_log_folder = "$log_folder";
my $hisat2_summary_folder = "$summary_folder";
my $hisat2_out_folder = "$detail_folder";

warn "Staring hisat2 with $progress threads ...\n";	

foreach my $basename ( sort keys %fq_parameters ) {
	warn "$basename\n";	
	#mkdir
	mkdir $hisat2_out_folder unless -d $hisat2_out_folder;
	mkdir $hisat2_log_folder unless -d $hisat2_log_folder;
	my $sample_hisat2_out_folder = "$hisat2_out_folder/$basename";
	mkdir $sample_hisat2_out_folder unless -d $sample_hisat2_out_folder;

	my $hisat2_command = "hisat2 -x $hisat2_index $fq_parameters{$basename} $hisat2_param 2> $hisat2_summary_folder/$basename.summary.txt | samtools view -F 4 - -b | samtools sort -o $sample_hisat2_out_folder/$basename.bam > $hisat2_log_folder/$basename.stdout.log 2> $hisat2_log_folder/$basename.stderr.log";
	my $index_command = "samtools index $sample_hisat2_out_folder/$basename.bam >> $hisat2_log_folder/$basename.stdout.log 2>> $hisat2_log_folder/$basename.stderr.log";
	$hisat2_command =~ s/\s+/ /g;
	warn "\t$hisat2_command\n";
	run_command($hisat2_command);
	run_command($index_command);
}

##################################################################
################reformat the files
##################################################################
my $bam_folder = "$out_folder/Bams";
mkdir $bam_folder unless -e $bam_folder;

foreach my $basename ( sort keys %fq_parameters ) {
	move("$hisat2_out_folder/$basename/$basename.bam", "$bam_folder/$basename.bam");
	move("$hisat2_out_folder/$basename/$basename.bam.bai", "$bam_folder/$basename.bai");
}
	
sub make_pass_dir {
	my $lable = shift;
	
}

sub run_command {
	my $command = shift;
	my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
	  run( command => $command, verbose => 0 );
	if ($success) {
		warn "\tDone!\n";
	}
	else {
		warn "Something went wrong:\n@$stderr_buf";
	}		
}
