#!/usr/bin/perl


open IN1, "<tmp_Metabolism";
open IN2, "<tmp_Genetic_Information_Processing";
open IN3, "<tmp_Environmental_Information_Processing";
open IN4, "<tmp_Cellular_Processes";
open IN5, "<tmp_Organismal_Systems";
open IN6, "<tmp_Human_Diseases";

open OUT7, ">tmp2_Metabolism";
open OUT8, ">tmp2_Genetic_Information_Processing";
open OUT9, ">tmp2_Environmental_Information_Processing";
open OUT10, ">tmp2_Cellular_Processes";
open OUT11, ">tmp2_Organismal_Systems";
open OUT12, ">tmp2_Human_Diseases";

while (<IN1>) {
        if (/^(\w+)/../^(\w+)/) {
                print OUT7 ;
                }
                }
while (<IN2>) {
        if (/^(\w+)/../^(\w+)/) {
                print OUT8 ;
                }
                }


while (<IN3>) {
        if (/^(\w+)/../^(\w+)/) {
                print OUT9 ;
                }
                }


while (<IN4>) {
        if (/^(\w+)/../^(\w+)/) {
                print OUT10 ;
                }
                }

while (<IN5>) {
        if (/^(\w+)/../^(\w+)/) {
                print OUT11 ;
                }
                }

while (<IN6>) {
        if (/^(\w+)/../^(\w+)/) {
                print OUT12 ;
                }
                }

