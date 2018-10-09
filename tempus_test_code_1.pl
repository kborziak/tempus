#!/usr/bin/perl -w
use strict;

##
#
#  Perl script to parse vcf files and query ExAC API for variant information
#  Written by Kirill Borziak 10/8/2018
#
##

# Import libraries
use JSON::Parse ':all';
use LWP::UserAgent;
use HTTP::Request;
my $ua = new LWP::UserAgent;
$ua->agent("Perl API Client/1.0");


# Variant consequence from most severe to least
my $consequences = {
  '3_prime_UTR_variant' => 10, 
  '5_prime_UTR_variant' => 11, 
  'initiator_codon_variant' => 7, 
  'intron_variant' => 9, 
  'missense_variant' => 6, 
  'non_coding_transcript_exon_variant' => 8, 
  'splice_acceptor_variant' => 3, 
  'splice_donor_variant' => 4, 
  'splice_region_variant' => 2, 
  'stop_gained' => 1, 
  'stop_lost' => 5, 
  'stop_retained_variant' => 12, 
  'synonymous_variant' => 13, 
  '' => 14, 
};


# Location of challenge data and output file
my $input_loc = "/home/kborziak/Documents/job search/tempus/Challenge_data\ (1).vcf";
my $output_loc = "/home/kborziak/Documents/job search/tempus/challenge_output";

# Open input file
open (my $input_file, "<", $input_loc);
# Open output file
open (my $output_file, ">", $output_loc);
print $output_file "chromosome\tposition\treference\tvariant\ttype_of_variation\tdepth\tvariant_reads\tpercentage_variant_reads\tExAC_allele_freq\tadditional_info\n";
# Read input file
while (my $line = <$input_file>) {
  chomp $line;
  # If not header
  if ($line !~ /^#/) {
    # Parse variant
    my @split = split (/\t/, $line);
    # Variant location, reference and alternate bp
    my $chr = $split[0];
    my $pos = $split[1];
    my $ref = $split[3];
    my $alt = $split[4];
    # Variant depth and alternate reads
    $split[7] =~ /DP=(\d+)/;
    my $depth = $1;
    $split[7] =~ /AO=(\d+)/;
    my $alt_reads = $1;
    # Percentage variant reads
    my $percentage = $alt_reads / $depth;
    $percentage = sprintf("%.4f", $percentage);
    # ExAC variables
    my $cons = '';
    my $allele_freq = 0;
    my $info = '';
    
    # Request ExAC API
    # Create URL
    my $url = "http://exac.hms.harvard.edu/rest/variant/$chr-$pos-$ref-$alt";
    # Fetch the ExAC API data from the URL
    my $request = HTTP::Request->new("GET" => $url);
    my $response = $ua->request($request);
    # Parse API output with JSON
    my $json_obj = parse_json ($response->content);
    
    # Get variant consequences
    if ($json_obj->{"consequence"}) {
      foreach my $key (keys $json_obj->{"consequence"}) {
        # Select most severe consequence
        if ($consequences->{$key} < $consequences->{$cons}) {
          $cons = $key;
        }
      }
      # Get allele freq
      $allele_freq = $json_obj->{"variant"}->{"allele_freq"};
    }
    # Else show variant not in ExAC
    else {
      $split[7] =~ /TYPE=([a-z]+)/;
      $cons = $1;
      $info = 'not in ExAC';
    }
    
    # Print output
    print $output_file "$chr\t$pos\t$ref\t$alt\t$cons\t$depth\t$alt_reads\t$percentage\t$allele_freq\t$info\n";
  }
}
close $input_file;

