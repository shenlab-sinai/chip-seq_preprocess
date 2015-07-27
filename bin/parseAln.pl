#!/usr/bin/perl

use strict;

my $infile = $ARGV[0] or die "require input bam (or sam) filename\n"; 
my $mapcut = $ARGV[1] or die "require mapq cutoff\n";
my $diffcut = $ARGV[2] or die "require AS (alignment score) - XS (next best alignment score) difference cutoff\n";

print STDERR "Splitting $infile with MAPQ >= $mapcut and AS-XS >= $diffcut.\n\n";


my $bamfile;
my $samfile;
if ($infile =~ /\.bam$/gm) {
   $bamfile = $infile;
   ($samfile = $bamfile) =~ s/\.bam$/\.sam/gm;
   print STDERR "Converting input bam to sam...";
   system("samtools view -h $bamfile > $samfile");
   print STDERR "Done!\n\n";
}elsif ($infile =~ /\.sam$/gm) {
   $samfile = $infile;
   ($bamfile = $samfile) =~ s/\.sam$/\.bam/gm;
   print STDERR "Converting input sam to bam...";
   system("samtools view -Sb $samfile > $bamfile");
   print STDERR "Done!\n\n";
}else {
   die "Unrecognized input file, require .sam or .bam files";
}   
(my $uniqmapped = $samfile) =~ s/\.sam$/\.uniqmapped\.sam/gm;
(my $multimapped = $samfile) =~ s/\.sam$/\.multimapped\.sam/gm;
(my $unmapped = $samfile) =~ s/\.sam$/\.unmapped\.sam/gm;
(my $uniqmappedbam = $uniqmapped) =~ s/\.sam$/\.bam/gm;
(my $multimappedbam = $multimapped) =~ s/\.sam$/\.bam/gm;
(my $unmappedbam = $unmapped) =~ s/\.sam$/\.bam/gm;

print STDERR "Splitting samfiles...";
unless (open(UNIQ, ">$uniqmapped")){
   print STDERR "Cannot open file \"$uniqmapped\" to write to!!\n\n";
   exit;
}
unless (open(MULTI, ">$multimapped")){
   print STDERR "Cannot open file \"$multimapped\" to write to!!\n\n";
   exit;
}
unless (open(UN, ">$unmapped")){
   print STDERR "Cannot open file \"$unmapped\" to write to!!\n\n";
   exit;
}

my $countline=0;
my $counthead=0;
my $countuniq=0;
my $countmulti=0;
my $countunmap=0;
open (INPUT,"$samfile") or die "Cannot open input file 1!!! \n";
while (my $line=<INPUT>) {
   chomp $line;
   my $line2=$line;
   my $line3=$line;
   my $line4=$line;
   my $line5=$line;
   my $line6=$line;
   $countline++;
   print STDERR "\t$countline" if ($countline%10000==0);
   my $header='NA';
   my $flag='NA';
   my $MAPQ='NA';
   my $AS='NA';
   my $XS='NA';
   my $DIFF='NA';

   $header='T' if ($line2 =~ /^@/gm);
   $MAPQ = $1 if ($line3 =~ /^\S+\s+\d+\s+\S+\s+\d+\s+(\d+)/gm);
   $flag = $1 if ($line4 =~ /^\S+\s+(\S+)/gm);
   $AS = $1 if ($line5 =~ /AS:i:(\S+)/gm);
   $XS = $1 if ($line6 =~ /XS:i:(\S+)/gm);
   $DIFF = $AS-$XS if (($AS ne 'NA')&&($XS ne 'NA'));

   #print STDERR "$header\t$flag\t$MAPQ\t$AS\t$XS\t$DIFF\n";

   if (($mapcut eq "default")||($diffcut eq "default")) {
      if ($header eq 'T') {
         $counthead++;
         print UNIQ "$line\n";
         print MULTI "$line\n";
         print UN "$line\n";
      }elsif ($flag==4) {
         $countunmap++;   
         print UN "$line\n";
      }elsif (($AS ne 'NA')&&($XS eq 'NA')) {
         $countuniq++;
         print UNIQ "$line\n";
      }elsif (($AS ne 'NA')&&($XS ne 'NA')) {
         $countmulti++;
         print MULTI "$line\n";
      }else {
         die "unknown combination: $header\t$flag\t$AS\t$XS\t$DIFF\nfrom line: $line\n";
      }
   }else {
      if ($header eq 'T') {
         $counthead++;
         print UNIQ "$line\n";
         print MULTI "$line\n";
         print UN "$line\n";
      }elsif ($flag==4) {
         $countunmap++;   
         print UN "$line\n";
      }elsif (($AS ne 'NA')&&($XS eq 'NA')) {
         if ($MAPQ >= $mapcut) {
            $countuniq++;
            print UNIQ "$line\n";
         }else {
            $countmulti++;
            print MULTI "$line\n";
         }
      }elsif (($AS ne 'NA')&&($XS ne 'NA')) {
         if (($MAPQ >= $mapcut)&&($DIFF >= $diffcut)) {
            $countuniq++;
            print UNIQ "$line\n";
         }else {
            $countmulti++;
            print MULTI "$line\n";
         }
      }else {
         die "unknown combination: $header\t$flag\t$AS\t$XS\t$DIFF\nfrom line: $line\n";
      }
   }

   #last if ($countline==50);
}
print STDERR "\tDone!\n";
print STDERR "$countline lines... $counthead headers, $countuniq unique, $countmulti multi, $countunmap unmapped.\n\n";

print STDERR "Converting sam to bam...";
system("samtools view -Sb $uniqmapped > $uniqmappedbam");
system("samtools view -Sb $multimapped > $multimappedbam");
system("samtools view -Sb $unmapped > $unmappedbam");
print STDERR "Done!\n";

print STDERR "Cleaning up part 1...";
system("rm $uniqmapped $multimapped $unmapped");
print STDERR "Done!\n\n";

print STDERR "Sorting bam...";
my @tosortfiles=($bamfile,$unmappedbam,$uniqmappedbam,$multimappedbam);
foreach my $file (@tosortfiles) {
   system("samtools sort -m 5G $file $file.sorted"); 
   system("mv $file.sorted.bam $file");
   system("samtools index $file");
}
print STDERR "Done!\n\n";
print STDERR "Sorted Bam files $uniqmappedbam, $multimappedbam, $unmappedbam created...\n\n";
print STDERR "parseAln.pl Completed!\n\n";

