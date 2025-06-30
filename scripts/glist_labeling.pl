#!/usr/bin/perl
use warnings;
use diagnostics;
if(@ARGV < 2){
	die "*.pl <gtdb_r214_all_ac2species.tsv> <L3K1*_gtdb_r214.glist> <[^]wanted species name>" ;
}
open $st,$ARGV[0] || die "can't open $ARGV[0]:$!";
while(<$st>){
	chomp;
	($ac,$sp)=(split '\t')[0,1];
	$sp =~ s/\s/_/g;
	$hash{$ac} = $sp;
	$ac =~ tr/AF/FA/ ;
	$hash{$ac} = $sp;
}
close $st;

#please remove size column if present
open $glist,$ARGV[1]  || die "can't open $ARGV[1]:$!";
$i = 1 ;
while(<$glist>){
	chomp;
	if(/(GC[AF]_\d+\.\d+)/){
		$ac = $1 ;
	}
	else{
		die "ac not matched! : $_";
	}

	if (!exists $hash{$ac} ){ # filter non gtdb
		print "0\tNULL\t$ac\n"; 
	}
	elsif (exists $seen{$ac} ){ # dedup
		print  "0\tNULL\t$ac\t$hash{$ac}\n";	
	}
	else {
		if (defined $ARGV[2]) {
			$sp = $ARGV[2];
    # Check if $ARGV[2] starts with '!' and remove it if present
				
    	 $negate = ($sp =~ s/^\^//);

    	# Determine whether to print "0 NULL" or the index based on condition
    	 $output = ( $negate >0 xor ($hash{$ac} eq $sp )) ? $i: "0\tNULL" ;

    	# Print the formatted output
    	print  "$output\t$ac\t$hash{$ac}\n";
		}
		else{
			print  $i,"\t",$ac,"\t",$hash{$ac},"\n";
		}
		$seen{$ac} = 1;
   	$i++;
	}
	
}

close $glist;


