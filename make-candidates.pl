#! /usr/bin/env perl

use warnings;
use strict;

#use Data::Printer;

use Array::Circular;
use Bio::SeqFeature::Primer;

my ($fasta, $gff) = @ARGV;
use Excel::Writer::XLSX;

my %stats;

my %ori_check = ( '+' => undef, '-' => undef );

# read in the genome
my $genome = load_fasta($fasta);
warn "processed genome\n";

my $cds = load_gff($gff);
warn "processed gff\n";

my %output;

my @sorted_cds;

for my $chr ( qw( CP030788.1 CP030789.1 CP030790.1 ) ) {

  my $prev_start = 0;

  my %open;

  resolve_overlaps($cds->{$chr});

  for my $feat ( @{$cds->{$chr}} ) {
  
    my $cds = $feat->{cds};
    my $start = $feat->{t_start};
    my $stop = $feat->{t_stop};
    
    $output{$cds} = $feat;
    push @sorted_cds, $cds;

# order of parameters:
#   chromosome circular array
#   position on genome of truncated beginning or end of gene
#   offset from gene where primer starts
#   is the primer revcomp?
#   which direction to we go if Tm is too low?

    # p 1a (on + strand)
    $output{$cds}{p1a} = primer_seq($genome->{$chr}, $start, -3000, '+', '+');

    # p 1b (on - strand)
    $output{$cds}{p1b} = primer_seq($genome->{$chr}, $start, -1, '-', '+');

    # p 2a (on + strand)
    $output{$cds}{p2a} = primer_seq($genome->{$chr}, $stop, 0, '+', '-');

    # p 2b (on - strand)
    $output{$cds}{p2b} = primer_seq($genome->{$chr}, $stop, 3000, '-', '-');
    
    # v (on - strand)
    $output{$cds}{v} = primer_seq($genome->{$chr}, $stop, 250, '-', '-');
  }
}

warn "processed CDS\n";

my $wb = Excel::Writer::XLSX->new("Primer-Design-variable-length-and-temp-primers.xlsx");

my $f_number = $wb->add_format();
$f_number->set_num_format( '#,##0' );

my $d_number = $wb->add_format();
$d_number->set_num_format( '#0.0' );

my $p_number = $wb->add_format();
$p_number->set_num_format('.00%' );

my $f_center = $wb->add_format( align => 'center' );

my %font = ( 
  -font => 'Lucida Console',
  );

my $fw_font = $wb->add_format(%font);

my $ws = $wb->add_worksheet();

my ($row,$col) = (0,10);

$ws->merge_range($row, $col, 0, $col+3, 'Primer 1a', $f_center);
$col += 5;

$ws->merge_range($row, $col, 0, $col+3, 'Primer 1b', $f_center);
$col += 5;

$ws->merge_range($row, $col, 0, $col+3, 'Primer 2a', $f_center);
$col += 5;

$ws->merge_range($row, $col, 0, $col+3, 'Primer 2b', $f_center);
$col += 5;

$ws->merge_range($row, $col, 0, $col+3, 'Validation Primer', $f_center);
$col += 5;
 
$row++;
$col=0;

$ws->write($row,$col++, $_) for ( 'CDS', 'Gene', 'Chr', 'Start', 'Trunc Start', 'Stop', 'Trunc Stop', 'CDS Length', 'Del Length', 'Del %' );

for ( 1 .. 5 ) {
  $ws->write($row,$col++, $_) for qw( Primer Length Tm Start Shift );
}

# set some column widths
$ws->set_column(0,2, 12);

$ws->set_column(10,10,40);
$ws->set_column(15,15,40);
$ws->set_column(20,20,40);
$ws->set_column(25,25,40);
$ws->set_column(30,30,40);


for my $cds ( @sorted_cds ) {
  my $data = $output{$cds};

  $row++;
  $col=0;

  $data->{t_start} or warn "weird t_start for $cds!";
  $data->{start} or warn "weird start for $cds!";

  $data->{t_stop} or warn "weird t_stop for $cds!";
  $data->{stop} or warn "weird stop for $cds!";

  # calculate CDS length
  my $cds_length = $data->{stop} - $data->{start} + 1;

  # calculate how much was deleted
  my $del_length = ($data->{p2a}{start}-1) - ($data->{t_start} + $data->{p1b}{shift}) + 1;

  $data->{t_start} = '' if $data->{t_start} == $data->{start};
  $data->{t_stop} = '' if $data->{t_stop} == $data->{stop};

  # start
  $ws->write($row,$col++, $_) for ( $cds, $data->{gene}, $data->{chr} );

  $ws->write($row,$col++, $_, $f_number) for map { $data->{$_} } qw( start t_start stop t_stop );

  $ws->write($row,$col++, $cds_length, $f_number);
  $ws->write($row,$col++, $del_length, $f_number);
  $ws->write($row,$col++, $del_length/$cds_length, $p_number);

  for ( qw( p1a p1b p2a p2b v ) ) {
 
    $ws->write($row,$col++, $data->{$_}{seq}, $fw_font);
    $ws->write($row,$col++, $data->{$_}{length}, $f_number);
    $ws->write($row,$col++, $data->{$_}{tm}, $d_number);
    $ws->write($row,$col++, $data->{$_}{start}, $f_number);
    $ws->write($row,$col++, $data->{$_}{shift});
  }

}

$wb->close();

#p %stats;

############################################################################################
#
# hashref primer_seq( $chr, $index, $offset, $ori, $shift_dir);
#
# find a valid primer meetings the requirements and return the data for it
#  $chr = Chromosome circular array
#  $index = array index of first base in primer
#  $offset = how far to move from the index
#  $ori = '+' / '-' direction to grow primer
#  $shift = '+' / '-' direction to shift primer
#
############################################################################################
sub primer_seq {

  my ($chr, $index, $offset,  $ori, $shift_dir) = @_;

  # move the index
  $chr->index($index);

  if ( $offset > 0 ) {
    $chr->next($offset);
  } elsif ( $offset < 0 ) {
    $chr->previous(abs($offset));
  }

  my $p_index = $chr->index();

  my %results;

  my $shift = 0;

  my $found_primer = call_position($chr, $p_index, $ori, \%results );

  while ( ! $found_primer ) {

    # keep track of shifting
    $shift++;

    if ( $shift_dir eq '+' ) {
      $chr->next();
    } else {
      $chr->previous();
    }

    # reset pointer
    $p_index = $chr->index();

    $found_primer = call_position($chr, $p_index, $ori, \%results);
  }

  $results{shift} = $shift;

  return \%results;
}

############################################################################################
#
# bool call_position( $chr, $loc, $ori, $results );
#
# Try to call a primer from the supplied position with the supplied
# orientation. If successful, information on the barcode is loaded
# into $results.
#
#  $chr = Chromosome circular array
#  $loc = array index of first base in primer
#  $ori = '+' / '-' direction to grow primer
#  $results = hash reference to store results if we find a good primer
#
#
############################################################################################
sub call_position {

  my ( $chr, $loc, $ori, $results ) = @_;

  for my $length ( qw( 25 24 26 23 27 22 28 21 29 20 30 ) ) {
    
    # set the index
    $chr->index($loc);

    # let's be optimistic!
    my $is_good = 1;

    # start by getting the primer sequence
    my $seq = call_primer( $chr, $ori, $length );
    $chr->index($loc);

    $seq eq '' and die "weirdness at $ori - $length - $loc\n";

    my $primer = Bio::SeqFeature::Primer->new( -seq => $seq );
  
    $stats{tm_check}++;

    # skip homopolymers of 5 or more
    $seq =~ /AAAAA/ and $is_good = 0;
    $seq =~ /TTTTT/ and $is_good = 0;
    $seq =~ /CCCCC/ and $is_good = 0;
    $seq =~ /GGGGG/ and $is_good = 0;

    # skip primers ending in A or T
    $seq =~ /[AT]$/ and $is_good = 0;

    $primer->Tm < 58 and $is_good = 0;
    $primer->Tm > 63 and $is_good = 0;

    if ( $is_good ) {

      $results->{tm} = $primer->Tm;
      $results->{seq} = $seq;
      $results->{length} = $length;
      
      # need to shift start 1 to the right to move between 0-index and 1-index
      $chr->next();

      # for negative primers do shenanigans to find genomic start
      if ( $ori eq '-' ) {
	$chr->previous($length-1);
      }
      
      # mark start and reseek
      $results->{start} = $chr->index();
      $chr->index($loc);
      return 1;
    }    
  }
  
  # if we fall out of the acceptable lengths, return false
  return 0;
  
}

############################################################################################
#
# string call_primer( $chr, $ori, $length );
#
# Call the primer from the given position with the given length and orientation.
#
############################################################################################
sub call_primer {

  my ($chr, $ori, $length) = @_;

  exists $ori_check{$ori} or die "$ori is a bad orientation";

  if ($ori eq '-') {
    $chr->previous($length-1);
  }
    
  my $string = $chr->current();
  
  for ( 1 .. $length-1 ) {
    $string .= $chr->next();
  }

  $string = revcomp($string) if $ori eq '-';

  return $string;
}

# comb through a chromosome looking for overlapping genes. If so keep track of how truncated they need to be

sub resolve_overlaps {

  my $chr = shift;

  my %open;
  my $prev_start = 0;

  for my $feat ( @$chr ) {

    my $start = $feat->{start};
    my $stop = $feat->{stop};
    my $cds = $feat->{cds};

    # sanity check, start position always increments
    $start > $prev_start or die "$cds is bad news - $start == $prev_start\n";

    # look for overlap with 5' end of all previous CDS
    while ( my ($o_cds, $o_feat) = each %open ) {

      # normal, non-overlapping
      if ( $o_feat->{stop} < $start ) {
	delete $open{$o_cds};
       
      # otherwise we have an overlap
      } else {

	# first we trim from the old feature
	my $new_o_trunc = $start-1;
	if ( $new_o_trunc < $o_feat->{t_stop} ) {
	  $o_feat->{t_stop} = $new_o_trunc;
	}

	# then we trim from the new feature
	my $new_trunc = $o_feat->{stop}+1;
	if ( $new_trunc > $feat->{t_start} ) {
	  $feat->{t_start} = $new_trunc;
	}
      }
    }

    # save this one as open
    $open{$feat->{cds}} = $feat;

  }
}

sub revcomp { # operates on all elements passed in
  my $revcomp = reverse(shift);
  $revcomp =~ tr/ACGTacgt/TGCAtgca/;
  return $revcomp;
}

sub load_gff {

  my $gff = shift;

  my %features;

  my %gene_name;

  open my $fh, '<', $gff or die "unable to open $gff: $!";

  while (<$fh>) {
    next if /^#/;

    chomp();
    next unless $_;
    my @f = split /\t/;

    if ( $f[2] eq 'gene' ) {
      my %p = map { split('=', $_) } split(';', $f[8]);
      $gene_name{$p{ID}} = $p{Name};
    }

    next unless $f[2] eq 'CDS';

    my %p = map { split('=', $_) } split(';', $f[8]);

    if ( ! exists $gene_name{$p{Parent}} ) {
      warn "haven't seen gene $p{Parent} for $p{ID}\n";
      next;
    }

    push @{$features{$f[0]}}, { start => $f[3], t_start => $f[3], stop => $f[4], t_stop => $f[4], cds => $p{Name}, gene => $gene_name{$p{Parent}}, chr => $f[0] };
  }

  return \%features;
}


sub load_fasta {
  
  my $fasta = shift;
  
  open my $fh, '<', $fasta or die "can't open $fasta: $!";

  
  my %tmp;
  my %genome;

  my $mol;

  while (<$fh>) {
    
    if ( />/ ) {
      
      /(CP\d+\.\d)/ or die "can't find chromosome name in $_";

      $mol = $1;

      $tmp{$mol} = [];
      
    } else {
      chomp();
      push @{$tmp{$mol}}, split('',$_);
    }
  }

  while ( my ($name, $seq) = each %tmp ) {

    my $a = Array::Circular->new(@$seq);
    $genome{$name} = $a;
  }

  return \%genome;
}
