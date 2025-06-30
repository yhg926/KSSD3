#!/usr/bin/perl
use strict;
use warnings;
use JSON;

# Check if the input file is provided
my $filename = $ARGV[0] or die "Usage: perl extract_edge_num.pl <file.jplace>\n";

# Read the JSON content from the file
open my $fh, '<', $filename or die "Could not open '$filename' $!\n";
my $json_text = do { local $/; <$fh> };
close $fh;

# Parse JSON
my $data = decode_json($json_text);

# Check if placements and fields are defined
die "No placements found in the file\n" unless $data->{placements};
die "No fields specified in the file\n" unless $data->{fields};

# Find the index of "edge_num" in the fields array
my ($edge_num_index) = grep { $data->{fields}[$_] eq 'edge_num' } 0..$#{$data->{fields}};
die "'edge_num' not found in fields\n" unless defined $edge_num_index;

# Extract "edge_num" for each query
for my $placement (@{$data->{placements}}) {
    my $query_name = $placement->{n}[0];  # Get the query name
    my @edge_nums = map { $_->[$edge_num_index] } @{$placement->{p}};
    print "Query: $query_name, Edge Numbers: @edge_nums\n";
}

