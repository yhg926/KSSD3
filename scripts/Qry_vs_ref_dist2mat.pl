use warnings;
use diagnostics;
if(@ARGV < 2){
	die "*.pl <distance.out> <dist col_num> <skip_head_n>";
}
open $fh,$ARGV[0] || die "can't open $ARGV[0]:$!";
$skipn = (!defined $ARGV[2] )? 0 : ($ARGV[2] =~ /^[1-9]\d*$/ ?  $ARGV[2] : 0) ;

#print $ARGV[2],"\t", $skipn,"\n"; exit(1);
while ( $skipn--)
{
	<$fh>;
} 

$qry_ind = $ref_ind = 0;
while(<$fh>){
	chomp;
	($qry,$ref,$dist) = (split /\t+/)[0,1,$ARGV[1]-1];
	
		 
  if (!exists $qryhash{$qry}){
		push @qry_names,$qry;
		$qryhash{$qry} = $qry_ind;
		$qry_ind++;
	}
	if(!exists $refhash{$ref}){
		push @ref_names,$ref;
		$refhash{$ref} = $ref_ind;
		$ref_ind++;
	}
 	$dist_arr[$qryhash{$qry}][$refhash{$ref}] = $dist;
}
close $fh;

for ($ref_ind=0; $ref_ind < @ref_names;$ref_ind++) {
	print "\t",$ref_names[$ref_ind];
}
print "\n";


for($qry_ind=0; $qry_ind < @qry_names;$qry_ind++){
	print $qry_names[$qry_ind];
	for ($ref_ind = 0; $ref_ind < @ref_names; $ref_ind++){
		print "\t",$dist_arr[$qry_ind][$ref_ind];
	}
	print "\n";
}




	

