#!/usr/bin/perl

#script to extract reads from the forward and the reverse strand
#needs samtools 
#works on the SAM Flag description:
#Bit    Description
#0x1    template having multiple segments in sequencing - single end or paired end
#0x40   the first segment in the template - first/forward read
#0x80   the last segment in the template - second/reverse read

#first argument is the BAM file location
#second argument is the strand rule - syntax borrowed from RSeQC
#third argument is the output folder where the new bam files will be created


#here the description from RSeQC:
#1) pair-end or single-end
#2) if experiment is strand-specific, how reads were stranded.
# * For pair-end RNA-seq, there are two different ways to strand reads:
#  i) 1++,1--,2+-,2-+
#     read1 mapped to '+' strand indicates parental gene on '+' strand
#     read1 mapped to '-' strand indicates parental gene on '-' strand
#     read2 mapped to '+' strand indicates parental gene on '-' strand
#     read2 mapped to '-' strand indicates parental gene on '+' strand
#  ii) 1+-,1-+,2++,2--
#     read1 mapped to '+' strand indicates parental gene on '-' strand
#     read1 mapped to '-' strand indicates parental gene on '+' strand
#     read2 mapped to '+' strand indicates parental gene on '+' strand
#     read2 mapped to '-' strand indicates parental gene on '-' strand
# * For single-end RNA-seq, there are two different ways to strand reads:
#  i) ++,--
#     read mapped to '+' strand indicates parental gene on '+' strand
#     read mapped to '-' strand indicates parental gene on '-' strand
#  ii) +-,-+
#     read mapped to '+' strand indicates parental gene on '-' strand
#     read mapped to '-' strand indicates parental gene on '+' strand		
#
# if unsure use infer_experiment.py from RsEQC to determine strand rule

use File::Basename;

my @infos;
my @sam_infos;
my @char;
my %change;

#determine strand rule
@infos = split (/\,/, $ARGV[1]);
if ($#infos == 1)
{
    #single end sequencing
    
    #what happens to a read mapping to the forward strand
    @char = split (//, $infos[0]);
    $change{$char[0]} = $char[1];
    
    #what happens to a read mapping to the reverse strand
    @char = split (//, $infos[1]);
    $change{$char[0]} = $char[1];
   
}
else
{
    #paired end sequencing
    
    #first in pair
    #what happens to a read mapping to the forward strand
    @char = split (//, $infos[0]);
    $strand = $char[0].$char[1];
    $change{$strand} = $char[2];
    
    #what happens to a read mapping to the reverse strand
    @char = split (//, $infos[1]);
    $strand = $char[0].$char[1];
    $change{$strand} = $char[2];

    #second in pair
    #what happens to a read mapping to the forward strand
    @char = split (//, $infos[2]);
    $strand = $char[0].$char[1];
    $change{$strand} = $char[2];

    #what happens to a read mapping to the reverse strand
    @char = split (//, $infos[3]);
    $strand = $char[0].$char[1];
    $change{$strand} = $char[2];
   
}


#writing BAM header
#this is done once here not to test each line for being a header later

open BAMREAD, "samtools view -H $ARGV[0] |";
$fwd_bam = $ARGV[2]."/".basename($ARGV[0],".bam")."_fwd.bam";
$rev_bam = $ARGV[2]."/".basename($ARGV[0],".bam")."_rev.bam";

open FWDWRITE, "| samtools view -bS -> $fwd_bam";
open REVWRITE, "| samtools view -bS -> $rev_bam";

while (<BAMREAD>)
{
    print FWDWRITE $_;
    print REVWRITE $_;
}
close (BAMREAD);

#now processing through the BAM file
open BAMREAD, "samtools view $ARGV[0] |";

while (<BAMREAD>)
{
    @sam_infos = split (/\t/, $_);
 
    #determining strand of alignment
    if ($sam_infos[1] & 0x10)
    {
        #reverse strand
        $strand = "-";
    }
    else
    {
        $strand = "+";
    }

    
    if ($sam_infos[1] & 0x01)
    {
        #paired end sequencing - add "1" or "2" to the strand
        if ($sam_infos[1] & 0x40)
        {
            #first read
            $strand = "1".$strand;
        }
        else
        {
            #second read
            $strand = "2".$strand;
        }
    }
    
    #writing to correct pipe

    if ($change{$strand} =~ /\+/)
    {
        print FWDWRITE $_; 
    }
    else
    {
        print REVWRITE $_;
    }

}
close (BAMREAD);
close (FWDWRITE);
close (REVWRITE);