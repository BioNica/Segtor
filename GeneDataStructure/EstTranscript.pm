=head1 NAME

   GeneDataStructure::EstTranscript.pm


=head1 Disclaimer

This software and all the data therein were developed by the Bioinformatics and Computational Biology Lab (Laboratorio de Bioinformatica e Biologia Computational (LBBC)) of the Brazilian National Cancer Institute (Instituto Nacional de Cancer (INCA)) in Rio de Janeiro, RJ, Brazil. The INCA is part of the Health Ministry (Ministerio da Saude) of the Federal Government of Brazil. This software was developed by the staff of the LBBC as part of their official duties. According to federal regulations, this software is therefore public domain.

This software is provided as is and the INCA does make any guarantee regarding the accuracy/quality/reliability of the results. The INCA cannot be held accountable for any hardware or software damage that might be caused through the use or misuse of this software.

The LBBC encourages researchers using our software to acknowledge us. The software can be distributed freely and modified given that any derivative bears an acknowledgment that it was derived from our original software.

=head1 Author

Gabriel Renaud grenaud@inca.gov.br



=cut

package EstTranscript;

use strict;
use warnings;

use Data::Dumper;
use Transcript;
our @ISA = "Transcript";

=item new

   Constructor for the EstTranscript object. This object has Transcript.pm
   as parent class. This object is used to hold data for an EST.

   The only parameter that is required is the 'databaseline' which is the
   line from the all_est table from UCSC.

=cut
sub new{
  my ($class,%arg)=(@_);

  if( !(exists $arg{databaseline}) ){
    die "Enter the 'databaseline' parameter\n";
  }

  my $self;

  my $chromosome;
  my $strand;
  my $id;
  my $transcStart;
  my $transcEnd;

  my $codingStart;
  my $codingEnd;

  my $exonsStart;
  my $exonsSize;

  my $querySize;
  my $queryStart;
  my $queryEnd;


  #1    2       3       4       5       6       7       8       9      10       11              12     13       14      15      16              17      18     19       20                      21                      22
  #585	330	13	0	0	2	2	3	4	-	AA663731	346	1	346	chr1	247249719	2802	3149	6	47,26,73,27,165,5,	0,47,74,147,174,340,	2802,2850,2876,2951,2979,3144,


  #                             1       2       3       4       5       6       7       8       9      10         11      12      13      14      15      16      17      18       19        20        21          22
  if ( $arg{databaseline} =~ /^(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([\+\-])\s+(\w+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([\.\w]+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([,\d]+)\s+([,\d]+)\s+([,\d]+)/ ) {

    $chromosome   = $15;
    $strand       = $10;
    $id           = $11;
    $transcStart  = $17;
    $transcEnd    = $18;
    $exonsStart   = $22;
    $exonsSize    = $20;
    $querySize    = $12;
    $queryStart   = $13;
    $queryEnd     = $14;

  } else {
    die "Database line ".$arg{databaseline}." did not parse\n";
  }

  $self = $class->SUPER::new(chromosome   => $chromosome,
			     strand       => $strand,
			     id           => $id,
			     transcStart  => $transcStart,
			     transcEnd    => $transcEnd,
			     exonsStart   => $exonsStart,
			     exonsSize    => $exonsSize,
			     querySize    => $querySize,
			     queryStart   => $queryStart,
			     queryEnd     => $queryEnd,
			     baseDir      => $arg{baseDir});

  return $self;
}




1;
