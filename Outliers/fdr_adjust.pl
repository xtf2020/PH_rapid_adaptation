#!/usr/bin/perl
$p_file = $ARGV[0];
$p_colum1 = $ARGV[1];
$out = $ARGV[2];
$p_colum = $p_colum1-1;
open(IN,"<$p_file");
open(STDOUT, ">$out");
#readline IN;
while(<IN>){
	chomp;
	s/^\s+//;
	s/\s+$//;
	s/\s+/\t/g;
	$old_p_pos = $.;
	$old_p_line = $_;
	$old_p{$old_p_pos} = $old_p_line;
	@old_p_list = split (/\t/,$old_p_line);
	$old{$old_p_pos} = $old_p_list[$p_colum];
}
@rank_p = sort {$old{$a} <=> $old{$b} or $a <=> $b} keys %old;
$i = 1;
foreach (@rank_p) {
	$rank_p = $_;
	# print "$rank_p\t$old{$rank_p}\n";
	$p_p{$i} = $rank_p;
	$np_op{$rank_p} = $i;
	$i++;
}
$lenth = $i-1;
$n = $lenth;
while($n >= 1){
	$new_p_value = ($old{$p_p{$n}} * $lenth / $n);
	$new{$n} = $new_p_value;
	unless ($n == $lenth){
		$max = $n+1;
		unless($new{$n} <= $new{$max}){
			$new{$n} = $new{$max};
		}
	}
	$n = $n-1;
}
$o = 1;
while ($o <= $lenth){
	print "$old_p{$o}\t$new{$np_op{$o}}\n";
	$o++;
}