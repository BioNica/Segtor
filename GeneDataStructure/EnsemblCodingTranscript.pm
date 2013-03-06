=head1 NAME

   GeneDataStructure::EnsemblCodingTranscript.pm


=head1 Disclaimer

This software and all the data therein were developed by the Bioinformatics and Computational Biology Lab (Laboratorio de Bioinformatica e Biologia Computational (LBBC)) of the Brazilian National Cancer Institute (Instituto Nacional de Cancer (INCA)) in Rio de Janeiro, RJ, Brazil. The INCA is part of the Health Ministry (Ministerio da Saude) of the Federal Government of Brazil. This software was developed by the staff of the LBBC as part of their official duties. According to federal regulations, this software is therefore public domain.

This software is provided as is and the INCA does make any guarantee regarding the accuracy/quality/reliability of the results. The INCA cannot be held accountable for any hardware or software damage that might be caused through the use or misuse of this software.

The LBBC encourages researchers using our software to acknowledge us. The software can be distributed freely and modified given that any derivative bears an acknowledgment that it was derived from our original software.

=head1 Author

Gabriel Renaud grenaud@inca.gov.br



=cut

package EnsemblCodingTranscript;

use strict;
use warnings;

use Data::Dumper;
use Range;
use CodingTranscript;
our @ISA = "CodingTranscript";


=item new

   Constructor for the EnsemblCodingTranscript object. This object inherits from
   CodingTranscript.pm and is used to store information pertaining to an Ensembl Transcript.
   This object can be constructed by providing a line from the ensGene.txt table from UCSC.

=cut
sub new{
  my ($class,%arg)=(@_);

  if( !(exists $arg{ensemblline}) ){
    die "Enter the 'ensemblline' parameter\n";
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
  my $exonsEnd;

  my $framesFromFile;

  #1       2               3       4       5       6       7       8       9       10                              11                              12      13                  14      15      16
  #585     ENST00000404059 chr1    +       1872    3533    3533    3533    6       1872,2041,2475,2837,3083,3315,  1920,2090,2560,2915,3237,3533,  0       ENSG00000219789     none    none    -1,-1,-1,-1,-1,-1o,
  #                             1       2       3        4         5       6       7       8       9       10         11        12       13     14      15     16
  if ( $arg{ensemblline} =~ /^(\d+)\s+(\S+)\s+(\w+)\s+([\+\-])\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([,\d]+)\s+([,\d]+)\s+(\d+)\s+(.+)\s+(\w+)\s+(\w+)\s+([\-,\d]+)$/ ) {

    $chromosome  = $3;
    $strand      = $4;
    $id          = $2;

    $transcStart = $5;
    $transcEnd   = $6;

    $codingStart = $7;
    $codingEnd   = $8;

    $exonsStart  = $10;
    $exonsEnd    = $11;

    $framesFromFile = $16;


  } else {
    die "Database line ".$arg{ensemblline}." did not parse\n";
  }



  $self = $class->SUPER::new(chromosome     => $chromosome,
			     strand         => $strand,
			     id             => $id,
			     transcStart    => $transcStart,
			     transcEnd      => $transcEnd,
			     codingStart    => $codingStart,
			     codingEnd      => $codingEnd,
			     exonsStart     => $exonsStart,
			     exonsEnd       => $exonsEnd,
			     baseDir        => $arg{baseDir},
			     seqExtractor   => $arg{seqExtractor},
			     framesFromFile => $framesFromFile);




  return $self;
}




1;
