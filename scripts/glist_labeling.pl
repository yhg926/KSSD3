use warnings;
use diagnostics;
if(@ARGV < 2){
	die "*.pl <gtdb_r214_all_ac2species.tsv> <L3K1*_gtdb_r214.glist> <wanted species name>" ;
}
open $st,$ARGV[0] || die "can't open $ARGV[0]:$!";
while(<$st>){
	chomp;
	($ac,$sp)=(split '\t')[0,1];
	$sp =~ s/\s/_/g;
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
		die "$ac not matched!";
	}

	if (!exists $hash{$ac} ){ # filter non gtdb
		print "0\tNULL\n"; 
	}
	elsif (exists $seen{$ac} ){ # dedup
		print "0\tNULL\t$hash{$ac}\n";	
	}
	elsif (defined $ARGV[2]){
		if($hash{$ac} eq $ARGV[2]){
			print $i,"\t",$ac,"\t",$ARGV[2],"\n"; 
		}
		else{
			 print "0\tNULL","\t",$ac,"\t",$hash{$ac},"\n";
		}
		$i++;
	}
	else{
		print $i,"\t",$ac,"\t",$hash{$ac},"\n";
		$seen{$ac} = 1;
    $i++;
	}
}
close $glist;
