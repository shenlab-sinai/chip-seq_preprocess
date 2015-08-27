#!/usr/bin/perl

use strict;

my $logdir = $ARGV[0] or die "require log_path\n"; 
my $fastqcdir = $ARGV[1] or die "require fastqc_path\n"; 
my $alignmentdir = $ARGV[2] or die "require alignment_path\n"; 
my $rmdupdir = $ARGV[3] or die "require rmdup_path\n"; 
my $tdfdir = $ARGV[4] or die "require tdf_path\n"; 
my $phantompeakdir = $ARGV[5] or die "require phantompeak_path\n"; 
my $diffrepeatdir = $ARGV[6] or die "require diffrepeat_path\n"; 

system("mkdir $logdir");
system("mkdir $fastqcdir");
system("mkdir $alignmentdir");
system("mkdir $rmdupdir");
system("mkdir $tdfdir");
system("mkdir $phantompeakdir");
system("mkdir $diffrepeatdir");


