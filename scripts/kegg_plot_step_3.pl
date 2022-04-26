#!/usr/bin/perl

open IN, "$ARGV[0]";

our @grabbed;
my $level2_count;
my $level2_name;

while (<IN>) {
	chomp;
        if (/^Metabolism$|^Genetic Information Processing|^Environmental Information Processing|^Cellular Processes|^Organismal Systems|^Human Diseases/) { 
	$level1_name = $_;

 }  elsif ($_ =~ /(^[A-Z])/)   {
	$level2_name = $_;                              
	push @grabbed, $_;		
	while (<IN>) {                               #这里使用的嵌套循环，第二个内循环会接着往下读取
	last if (/^[A-Z]/);
	if (/\((\d+)\)/) {
	$level2_count += $1;
}	
	
	push @grabbed, $_;
}
 
	 print "$level1_name"."\t"."$level2_name"."\t"."$level2_count\n"  ; 
	

	$level2_count = 0;
	
        #seek (IN, -length($line), 1);

}
}



