#!/usr/bin/env perl
# Mike Covington
# created: 2015-01-23
#
# Description:
#
use strict;
use warnings;
use Log::Reproducible;
use autodie;
use feature 'say';
use File::Path 'make_path';
use Getopt::Long;
use List::Util 'sum';
use POSIX;

my ( $id, $help ) ;
my $length_3_UTR = 500;
my $out_dir      = '.';

my $options = GetOptions(
    "id=s"           => \$id,
    "length_3_UTR=i" => \$length_3_UTR,
    "out_dir=s"      => \$out_dir,
    "help"           => \$help,
);

my @bam_files = @ARGV;

validate_options( $id, \@bam_files, $help );

my $coverage = get_coverage(@bam_files);
my $seq_lengths = get_seq_lengths( $coverage, $bam_files[0] );

my ($rel_cds_cov, $rel_utr_cov,    $abs_cds_cov,
    $abs_utr_cov, $scaled_cds_cov, $total_cov
) = build_coverage_profiles( $coverage, $seq_lengths, $length_3_UTR );

make_path $out_dir;
write_coverage_profile( $rel_cds_cov,    $total_cov, "$id.cds-rel.cov" );
write_coverage_profile( $rel_utr_cov,    $total_cov, "$id.utr-rel.cov" );
write_coverage_profile( $abs_cds_cov,    $total_cov, "$id.cds-abs.cov" );
write_coverage_profile( $abs_utr_cov,    $total_cov, "$id.utr-abs.cov" );
write_coverage_profile( $scaled_cds_cov, $total_cov, "$id.cds-scaled.cov" );

exit;

sub build_coverage_profiles {
    my ( $coverage, $seq_lengths, $length_3_UTR ) = @_;

    my @seq_length_values = values %$seq_lengths;
    my $mean_cds_length
        = ( sum(@seq_length_values) / scalar @seq_length_values )
        - $length_3_UTR;

    my ( %rel_cds_cov, %rel_utr_cov );
    my ( %abs_cds_cov, %abs_utr_cov );
    my %scaled_cds_cov;
    my $total_cov;

    for my $seq_id ( sort keys %$coverage ) {
        my $length  = $$seq_lengths{$seq_id};
        my $cds_end = $length - $length_3_UTR;

        for my $pos ( keys %{ $$coverage{$seq_id} } ) {
            my $depth = $$coverage{$seq_id}{$pos};

            $total_cov += $depth;

            if ( $pos <= $cds_end ) {

                my $cds_length = $length - $length_3_UTR;

                # last pos of cds == 0
                my $abs_pos = $pos - $cds_length;
                $abs_cds_cov{$abs_pos} += $depth;

                # scale cds position from -99 to 0
                my $rel_pos = ceil( 100 * $abs_pos / $cds_length );
                $rel_cds_cov{$rel_pos} += $depth;

                # scale cds position from -(mean cds length - 1) to 0
                my $scaled_pos
                    = ceil( $mean_cds_length * $abs_pos / $cds_length );
                $scaled_cds_cov{$scaled_pos} += $depth;
            }
            else {
                my $utr_pos = $pos - $cds_end;

                my $rel_pos = ceil( 100 * $utr_pos / $length_3_UTR );
                $rel_utr_cov{$rel_pos} += $depth;

                $abs_utr_cov{$utr_pos} += $depth;
            }
        }
    }

    return \%rel_cds_cov, \%rel_utr_cov, \%abs_cds_cov, \%abs_utr_cov,
        \%scaled_cds_cov, $total_cov;
}

sub do_or_die {
    my ( $help, $errors ) = @_;

    if ($help) {
        die usage();
    }
    elsif (@$errors) {
        my $error_string = join "\n", map {"ERROR: $_"} @$errors;
        die usage(), $error_string, "\n\n";
    }
}

sub get_coverage {
    my @bam_files = @_;
    my %coverage;

    for my $bam (@bam_files) {
        say "Processing $bam";
        open my $depth_fh, "-|", "samtools depth $bam";
        while (<$depth_fh>) {
            chomp;
            my ( $seq_id, $pos, $depth ) = split /\t/, $_;
            $coverage{$seq_id}{$pos} += $depth;
        }
    }

    return \%coverage;
}

sub get_seq_lengths {
    my ( $coverage, $bam ) = @_;
    my %seq_lengths;

    open my $header_fh, "-|", "samtools view -H $bam";
    while (<$header_fh>) {
        next unless /^\@SQ/;

        my ( $seq_id, $length ) = /\@SQ\tSN:(.+)\tLN:(\d+)/;

        if ( exists $$coverage{$seq_id} ) {
            $seq_lengths{$seq_id} = $length;
        }
    }

    return \%seq_lengths;
}

sub usage {
    return <<EOF;

Usage: $0 [options] --id <Sample ID> <bam file(s)>

Options:
  -i, --id              Sample ID (required)
  -l, --length_3_UTR    Length of 3'-UTR following CDS in reference [500]
  -o, --out_dir         Output directory [.]
  -h, --help            Display this usage information

EOF
}


sub validate_options {
    my ( $id, $bam_files, $help ) = @_;
    my @errors;

    push @errors, "Must specify Sample ID with '--id'" unless defined $id;
    push @errors, "Must specify at least one BAM file"
        if scalar @$bam_files == 0;

    do_or_die( $help, \@errors );
}

sub write_coverage_profile {
    my ( $coverage_profile, $total_cov, $cov_out_file ) = @_;

    open my $cov_out_fh, ">", $cov_out_file;
    for my $pos ( sort { $a <=> $b } keys %$coverage_profile ) {
        my $percent_coverage = 100 * $$coverage_profile{$pos} / $total_cov;
        say $cov_out_fh join "\t", $pos, $percent_coverage;
    }
    close $cov_out_fh;
}
