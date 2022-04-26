#!/usr/bin/perl

open IN, "$ARGV[0]";
open OUT1, ">tmp_Metabolism";
open OUT2, ">tmp_Genetic_Information_Processing";
open OUT3, ">tmp_Environmental_Information_Processing";
open OUT4, ">tmp_Cellular_Processes";
open OUT5, ">tmp_Organismal_Systems";
open OUT6, ">tmp_Human_Diseases";

#prepare different Level1 file

while (<IN>) {
	if (/Metabolism/../Genetic Information Processing/) {
	print OUT1 unless /Genetic Information Processing/;
} 
}
	seek (IN, 0, 0);

while (<IN>) {
	if (/Genetic Information Processing/../Environmental Information Processing/) {
	print OUT2 unless /Environmental Information Processing/;
}
}

	seek (IN, 0, 0);

while (<IN>) {
	if (/Environmental Information Processing/../Cellular Processes/) {
	print OUT3 unless /Cellular Processes/;
} 
}
	seek (IN, 0, 0);

while (<IN>) {
	if (/Cellular Processes/../Organismal Systems/) {
	print OUT4 unless /Organismal Systems/;
} 
}
	seek (IN, 0, 0);

while (<IN>) {
	if (/Organismal Systems/../Human Diseases/) {
	print OUT5 unless /Human Diseases/;
} 
}
	seek (IN, 0, 0);

while (<IN>) {
	if (/Human Diseases/../^$/) {
	print OUT6 unless /^$/;
} 
}
	seek (IN, 0, 0);

