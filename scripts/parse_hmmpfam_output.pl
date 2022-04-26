#!/usr/bin/env perl

=head1 NAME

    parse_hmmpfam_output.pl

=head1 SYNOPSIS
 
    parse_hmmpfam_output.pl hmmpfam protein_fasta evalue_cutoff total_num_hmms cegma_dir
        where hmmpfam is the hmmpfam output file,
              protein_fasta is the fasta file of proteins for your species,
              evalue_cutoff is the evalue cutoff to use,
              total_num_hmms is the number of HMMs searched against (usually 458 for CEGMA KOGs), xxx
              cegma_dir is the directory where CEGMA was installed.

=head1 DESCRIPTION

    This script takes a hmmpfam output file, and counts the number of queries that have hits with evalue <= evalue_cutoff,
    and number of HMMs hit with evalue <= evalue_cutoff.

=head1 VERSION
  
    Perl script last edited 24-Jul-2013.

=head1 CONTACT

    alc@sanger.ac.uk (Avril Coghlan)

=cut

# 
# Perl script parse_hmmpfam_output.pl
# Written by Avril Coghlan (alc@sanger.ac.uk)
# 24-Jul-13.
# Last edited 24-Jul-2013.
# SCRIPT SYNOPSIS: parse_hmmpfam_output.pl: given a hmmpfam output file, counts the number of queries with hits of evalue <= evalue_cutoff, and number of HMMs hit with evalue <= evalue_cutoff.
#
#------------------------------------------------------------------#

# CHECK IF THERE ARE THE CORRECT NUMBER OF COMMAND-LINE ARGUMENTS:

use strict;
use warnings;

# xxx
BEGIN { 
    unshift (@INC, '/nfs/users/nfs_a/alc/Documents/git/helminth_scripts/modules'); 
}

use HelminthGenomeAnalysis::AvrilHMMUtils; 
use HelminthGenomeAnalysis::AvrilFastaUtils; 
use Carp::Assert; # HAS THE assert() FUNCTION 
use Math::Round;

my $num_args               = $#ARGV + 1;
if ($num_args != 5)
{
    print "Usage of parse_hmmpfam_output.pl\n\n";
    print "perl parse_hmmpfam_output.pl <hmmpfam> <protein_fasta> <evalue_cutoff> <total_num_hmms> <cegma_dir>\n";
    print "where <hmmpfam> is the hmmpfam output file,\n";
    print "      <protein_fasta> is the fasta file of proteins for your species,\n";
    print "      <evalue_cutoff> is the e-value cutoff to use,\n";
    print "      <total_number_hmms> is the total number of HMMs searched against (usually 458 for CEGMA KOGs),\n"; # xxx don't need?
    print "      <cegma_dir> is the directory where CEGMA was installed\n";
    print "For example, >perl parse_hmmpfam_output.pl ANCCEY_cegma.output ANCCEY.protein.fa 1e-05 458\n";
    print "/nfs/users/nfs_m/mz3/bin/cegma_v2/\n";
    exit;
}

# FIND THE PATH TO THE INPUT HMMPFAM FILE:                     

my $hmmpfam                = $ARGV[0];

# FIND THE FASTA FILE OF PROTEINS FOR YOUR SPECIES:

my $protein_fasta          = $ARGV[1];

# FIND THE EVALUE CUTOFF TO USE:

my $evalue_cutoff          = $ARGV[2];

# FIND THE TOTAL NUMBER OF HMMs SEARCHED AGAINST:

my $total_num_hmms         = $ARGV[3];

# FIND THE DIRECTORY WHERE CEGMA WAS INSTALLED:

my $cegma_dir              = $ARGV[4];

#------------------------------------------------------------------#

# RUN THE MAIN PART OF THE CODE:

&run_main_program($hmmpfam,$protein_fasta,$evalue_cutoff,$total_num_hmms,$cegma_dir);

print STDERR "FINISHED.\n";

#------------------------------------------------------------------#

# RUN THE MAIN PART OF THE CODE:

sub run_main_program
{
   my $hmmpfam             = $_[0]; # THE INPUT HMMPFAM FILE
   my $protein_fasta       = $_[1]; # THE FASTA FILE OF PROTEIN SEQUENCES
   my $evalue_cutoff       = $_[2]; # THE EVALUE CUTOFF TO USE
   my $total_num_hmms      = $_[3]; # TOTAL NUMBER OF HMMs SEARCHED AGAINST xxx need?
   my $cegma_dir           = $_[4]; # DIRECTORY WHERE CEGMA WAS INSTALLED
   my $protein_fasta_obj;           # OBJECT FOR THE INPUT FASTA FILE OF PROTEIN SEQUENCES
   my $num_hmmpfam_queries;         # NUMBER OF QUERIES SEEN IN $hmmpfam
   my $num_protein_queries;         # NUMBER OF PROTEIN SEQUENCES IN $protein_fasta 
   my $num_signif_matches;          # NUMBER OF QUERIES WITH SIGNIFICANT hmmpfam HITS 
   my $num_hmms_hit;                # NUMBER OF HMMs THAT ARE HIT WITH SIGNIFICANT MATCHES 
   my $pc_hmms_hit;                 # PERCENT OF $total_num_hmms THAT ARE HIT WITH SIGNIFICANT MATCHES
   my $CEGMA;                       # HASH TABLE OF SINGLE GENE CEGMA FAMILIES 
   my $num_conserved_hmms;          # NUMBER OF SINGLE GENE CEGMA FAMILIES  
 
   # READ IN THE FILE OF 248 SINGLE GENE CEGMA FAMILIES:
   $CEGMA                  = HelminthGenomeAnalysis::AvrilHMMUtils::read_cegma_fams($cegma_dir);
   $num_conserved_hmms     = keys (%{$CEGMA});
   assert($num_conserved_hmms == 248);

   # READ IN THE SEQUENCES IN THE INPUT FASTA FILE OF PROTEIN SEQUENCES:
   $protein_fasta_obj      = HelminthGenomeAnalysis::AvrilFastaUtils->new(fasta_file => $protein_fasta); 

   # COUNT THE NUMBER OF SEQUENCES IN THE INPUT FASTA FILE OF PROTEIN SEQUENCES:
   $num_protein_queries    = HelminthGenomeAnalysis::AvrilFastaUtils::count_number_seqs($protein_fasta_obj);  
 
   # COUNT THE NUMBER OF QUERIES IN THE INPUT HMMPFAM FILE:
   $num_hmmpfam_queries    = HelminthGenomeAnalysis::AvrilHMMUtils::count_hmmpfam_queries($hmmpfam,$evalue_cutoff);

   # CHECK $num_protein_queries IS EQUAL TO $num_hmmpfam_queries:
   if ($num_protein_queries != $num_hmmpfam_queries)
   {
      print STDERR "ERROR: run_main_program: have $num_protein_queries in $protein_fasta, but $num_hmmpfam_queries queries in $hmmpfam, so hmmpfam may not have run to completion!\n";
   }
  
   # COUNT THE NUMBER OF QUERIES THAT HAVE SIGNIFICANT MATCHES IN THE $hmmpfam FILE, AND NUMBER OF HMMs THAT ARE HIT WITH SIGNIFICANT MATCHES:
   ($num_signif_matches,$num_hmms_hit) = HelminthGenomeAnalysis::AvrilHMMUtils::count_signif_matches($hmmpfam,$evalue_cutoff,$CEGMA);
   $pc_hmms_hit             = $num_hmms_hit*100/$num_conserved_hmms; 
   $pc_hmms_hit             = nearest(0.1, $pc_hmms_hit); # round to 0.1 precision
   print "$hmmpfam: $num_signif_matches query sequences have significant matches (with e-value <= $evalue_cutoff), and $num_hmms_hit ($pc_hmms_hit %) HMMs were hit with e-value <= $evalue_cutoff\n";
}

#------------------------------------------------------------------#


