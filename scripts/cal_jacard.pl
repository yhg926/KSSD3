#!/usr/bin/perl
use strict;
use warnings;

# Function to read a genome sequence from a FASTA file
sub read_fasta {
    my ($file) = @_;
    open(my $fh, "<", $file) or die "Cannot open file $file: $!\n";
    my $seq = "";
    while (my $line = <$fh>) {
        chomp $line;
        next if $line =~ /^>/;  # Skip header lines
        $seq .= $line;
    }
    close $fh;
    return $seq;
}

# Function to compute the reverse complement of a DNA sequence
sub rev_comp {
    my ($seq) = @_;
    $seq = reverse $seq;
    $seq =~ tr/ACGTacgt/TGCAtgca/;
    return $seq;
}

# Function to extract canonical k-mers from a genome
sub extract_kmers {
    my ($seq, $k) = @_;
    my %kmers;
    for (my $i = 0; $i <= length($seq) - $k; $i++) {
        my $kmer = substr($seq, $i, $k);
        my $rc_kmer = rev_comp($kmer);
        my $canonical = ($kmer lt $rc_kmer) ? $kmer : $rc_kmer;  # Canonical form
        $kmers{$canonical} = 1;
    }
    return %kmers;
}

# Function to compute Jaccard index and shared k-mer count
sub jaccard_index {
    my ($set1, $set2) = @_;
    my $intersection = 0;
    my $size1 = scalar(keys %$set1);
    my $size2 = scalar(keys %$set2);

    foreach my $kmer (keys %$set1) {
        if (exists $set2->{$kmer}) {
            $intersection++;
        }
    }
    my $union = $size1 + $size2 - $intersection;
    my $jaccard = ($union == 0) ? 0 : $intersection / $union;

    return ($intersection, $size1, $size2, $jaccard);
}

# Main execution
if (@ARGV < 2) {
    die "Usage: $0 <k> <genome1.fasta> <genome2.fasta> [genome3.fasta ...]\n";
}

my $k = shift @ARGV;
die "k must be an integer greater than 0\n" if $k !~ /^\d+$/ || $k <= 0;

print "Processing k=$k with multiple genomes: @ARGV\n";

# Read and process all genome files
my %genomes;
foreach my $file (@ARGV) {
    print "Reading $file...\n";
    my $seq = read_fasta($file);
    $genomes{$file} = { kmers => { extract_kmers($seq, $k) } };
}

# Compute pairwise Jaccard index
my @files = keys %genomes;
for my $i (0 .. $#files - 1) {
    for my $j ($i + 1 .. $#files) {
        my ($genome1, $genome2) = ($files[$i], $files[$j]);
        my ($shared_kmers, $size1, $size2, $jaccard) = jaccard_index($genomes{$genome1}{kmers}, $genomes{$genome2}{kmers});
        printf "%s %s %d %d %d %.6f\n", $genome1, $genome2, $shared_kmers, $size1, $size2, $jaccard;
    }
}

