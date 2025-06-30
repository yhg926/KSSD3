#!/usr/bin/perl
use strict;
use diagnostics;
use warnings;
if(@ARGV < 2){
	die '*.pl <trianle.mtx> <"trianlge"/"matrix"> <Start with Num of sample if set>';
}

my $shape = 0;
$shape = 1 if $ARGV[1] =~ /^m/ ;

die "$ARGV[0] is already trianlge" if @ARGV  == 2 and $shape == 0 ;

open(my $fh, '<', $ARGV[0]) or die "Can't open '$ARGV[0]': $!";
my @lines = <$fh>;
close $ARGV[0];

if($shape == 0 and @ARGV > 2){
	print $#lines + 1,"\n";
	#print $lines[3];
	for ( my $l=0 ; $l<@lines; $l++) {
		print $lines[$l];
	}
	exit(0);
}
#$num_smp = $#lines + 1;	



my @mtx ;
my @names;
for( my $i = 0; $i < @lines; $i++){
	chomp $lines[$i];
	my	@value = split /\t+/,$lines[$i];
	die "$ARGV[0]: line $i is empty\n" if @value < 1;
	push @names, shift @value; 	
	for(my $n = 0 ; $n < @value; $n++){
		$mtx[$i][$n] = $mtx[$n][$i] = $value[$n];
	}
	$mtx[$i][$i] = "0.000000" ;
}
my $num_smp = scalar @names;
print $num_smp,"\n" if @ARGV > 2 ;

for (my $i =0 ;$i< $num_smp;$i++){
	print $names[$i];
	for (my $j =0 ;$j < $num_smp;$j++ ){
		print "\t",$mtx[$i][$j];
	}
	print "\n";
}











