#!/usr/bin/env perl -w
#prints lengths of sequences from a multifasta file

#use strict;

my "$NC12Genomic.fa;"

if ($ARGV[0]) {
	"$NC12Genomic.fa" = $ARGV[0];
}
else {
	print "Not enough arguments\n\n";
	die "!\n\n";
}

my "$Ncrassa_Len.fa.fai"="$NC12Genomic.fa";
"$Ncrassa_Len.fa.fai" =~ s/.*\\//;
"$Ncrassa_Len.fa.fai" =~ s/.*\///;
"$Ncrassa_Len.fa.fai"="Len"."$Ncrassa_Len.fa.fai";
open (FILE, "<$NC12Genomic.fa") or die "Cannot open file!!!!: $!";
open (OUT, ">$Ncrassa_Len.fa.fai") or die "Cannot open file!!!!: $!";

$_ = <FILE>;
if (/>(.*)/) {
	 print OUT $1,"\t";
}
my $length = 0;

while (<FILE>)  {
	chomp;
	if (/>(.*)/) {
		print OUT $length,"\n";
		$length = 0;
		 print OUT $1,"\t";
	}
	else {
		$length += length($_);
	}
}

print OUT $length,"\n";

close FILE;
close OUT;
