#!/usr/bin/perl -w 
use strict;
if( @ARGV != 4 ) {
    print "Usage: perl join_ID.pl 1st-file 2nd-file common-id-in-1st-file common-id-in-2nd-file\n";
    exit 0;
}
my $fh1=shift @ARGV;
my $fh2=shift @ARGV;
my $id1=shift @ARGV;
my $id2=shift @ARGV;
open FH2,"<$fh2";

my $rec2;
my $temp=0;
while (<FH2>)
{
 chomp($_);
 my @rec2=split(/[\s]+/,$_);
  $rec2[$id2-1] =~ s# ##g;
  $rec2->{$rec2[$id2-1]}{$temp++}=$_;
}
close FH2;
open FH1,"<$fh1";
while (<FH1>)
{
  chomp($_);
  my @rec1=split(/[\s]+/,$_);
  $rec1[$id1-1] =~ s# ##g;
  if (exists ($rec2->{$rec1[$id1-1]}))
  {
     my @value = values (%{$rec2->{$rec1[$id1-1]}});
     for (my $j =0 ;$j <= $#value; $j++)
     {
      print $_."\t".$value[$j]."\n";  
     }
  }  
}
close FH1;
