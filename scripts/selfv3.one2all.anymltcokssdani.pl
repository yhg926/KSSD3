#!/usr/bin/perl
use warnings;
use diagnostics;
use Digest::MurmurHash3 qw( murmur32 );

if(@ARGV < 2){
	die "*.pl <query.fa> <ref1.fa> <ref2.fa> ...";
}

$UnitCL =  9; #20 is best
$UnitOL = 6; #16 is best
$UnitL = $UnitCL + $UnitOL ;
$UnitN = 2;  # 1 is best
$CTX = $UnitCL * $UnitN;
$TL = $UnitL*$UnitN + $UnitOL ;
#@nuqcount = (0) x @ARGV ;

$/ = '>';
for ($n=0 ; $n<@ARGV;$n++) {
	open ($f, $ARGV[$n] )|| die "can't open $ARGV[$n]:$! ";
	<$f>;

	while($seq=<$f>){
		chomp $seq;
		$seq =~ s/\n|\r//g; 
		next if length $seq < $TL ;

		for($i = 0; $i< (length $seq) - $TL +1 ; $i++){
			$tuple = substr($seq,$i,$TL);
			$ctx = $obj ="";
			for ($p = 0; $p < $TL; $p+=$UnitL ) {
				$obj.= substr($tuple,$p,$UnitOL) ;
				$ctx.= substr($tuple,$p+$UnitOL,$UnitCL) if $p+$UnitOL <  $TL;
			}
			$crctx = $ctx;
			$crctx =~ tr/ACGTacgt/TGCATGCA/ ; 
			$crctx = reverse $crctx;
	
			if ($ctx lt $crctx) {
				$unictx = $ctx;
			}
			elsif ($ctx gt $crctx) {
				$unictx = $crctx;
				$obj =~  tr /ACGTacgt/TGCATGCA/ ; 
				$obj = reverse $obj; 
			}
			else { next;
			}
			next if murmur32($unictx) % 256 != 1; 
		
		if ($n == 0) {	
			if(exists $hash[$n]->{$unictx}) {
				if ($hash[$n]->{$unictx} eq $obj ){					
					next;
				}
				elsif ($hash[$n]->{$unictx} ne 0 ){ 										
					$hash[$n]->{$unictx} = 0 ;
				}
			}
			else{
				$hash[$n]->{$unictx} = $obj; 	
			}
		}
	 else{
			$hash[$n]->{$unictx}->{$obj}++; 
		}
		}	

	}
 	close $f;
}

	$n = 0;
	$num_n = 0;
	$num_n_self_diff = 0;
	foreach	$ele (keys %{$hash[$n]}){
		  if ($hash[$n]->{$ele} ne 0){
				 $num_n++;
			}else{
				 $num_n_self_diff++;	
				delete $hash[$n]->{$ele};
			}
	}

	for ($j = 1; $j < @ARGV ; $j++){
		$num_j = 0 ;
		$num_ctx = 0;
		$diff_obj  = 0;
		$diff_obj_UnitN = 0;

		foreach $ele (keys %{$hash[$j]}){
			$num_j++;
			if( exists $hash[$n]->{$ele}  ){
				$num_ctx++;
#new 20250621
				@objs = keys %{$hash[$j]->{$ele}};	
				next if @objs > 1;	
				$obj_ref = $hash[$n]->{$ele};
				$obj_qry = $objs[0];	

				if ($obj_qry ne $obj_ref){
					$diff_obj++;
					for ($p = 0; $p < $UnitOL*$UnitN + 1 ; $p+=$UnitOL ){
						$diff_obj_UnitN++  if substr($obj_ref,$p,$UnitOL) ne substr($obj_qry,$p,$UnitOL);
					}
				#if ( !exists $hash[$j]->{$ele}->{$hash[$n]->{$ele}} ){
				#	$diff_obj++;
				#}
			}
		}
	}
		print $ARGV[$n],"\t", $ARGV[$j],"\t",$num_n,"\t",$num_n_self_diff,"\t",$num_j,"\t",$num_ctx,"\t",$diff_obj,"\t",$diff_obj_UnitN,"\t";
		if($num_ctx == 0) {
			print "0\t0\tnan\tnan\n";
		}
		else {
			print  $num_ctx/$num_n,"\t", $num_ctx/$num_j,"\t", $diff_obj / $num_ctx,"\t", (1-$diff_obj / ($num_ctx))**(1/($TL-$CTX)),"\n" ; 
	 }
	}



#my $string = "Hello, World!";
#my $hash = murmur32($string);

#print "MurmurHash3 (x86_32) hash: $hash\n";





