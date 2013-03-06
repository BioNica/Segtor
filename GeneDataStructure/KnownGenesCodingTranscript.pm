=head1 NAME

   GeneDataStructure::KnownGenesCodingTranscript.pm


=head1 Disclaimer

This software and all the data therein were developed by the Bioinformatics and Computational Biology Lab (Laboratorio de Bioinformatica e Biologia Computational (LBBC)) of the Brazilian National Cancer Institute (Instituto Nacional de Cancer (INCA)) in Rio de Janeiro, RJ, Brazil. The INCA is part of the Health Ministry (Ministerio da Saude) of the Federal Government of Brazil. This software was developed by the staff of the LBBC as part of their official duties. According to federal regulations, this software is therefore public domain.

This software is provided as is and the INCA does make any guarantee regarding the accuracy/quality/reliability of the results. The INCA cannot be held accountable for any hardware or software damage that might be caused through the use or misuse of this software.

The LBBC encourages researchers using our software to acknowledge us. The software can be distributed freely and modified given that any derivative bears an acknowledgment that it was derived from our original software.

=head1 Author

Gabriel Renaud grenaud@inca.gov.br



=cut

package KnownGenesCodingTranscript;

use strict;
use warnings;

use Data::Dumper;
use Range;
use CodingTranscript;
our @ISA = "CodingTranscript";

=item new

   Constructor for the KnownGenesCodingTranscript object. This object inherits from
   CodingTranscript.pm and is used to store information pertaining to an KnownGenes Transcript.
   This object can be constructed by providing a line from the knownGene.txt table from UCSC.

=cut
sub new{
  my ($class,%arg)=(@_);

  if ( !(exists $arg{knowngeneline}) ) {
    die "Enter the 'knowngeneline' parameter\n";
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

  my $proteinID;

  #1               2       3       4       5       6       7       8       9               10                      12
  #uc001aaa.2      chr1    +       1115    4121    1115    1115    3       1115,2475,3083, 2090,2584,4121,         uc001aaa.2

  #                                1           2        3        4       5       6       7       8       9           10              11          12
  if ( $arg{knowngeneline} =~ /^(\w+\.\d+)\s+(\w+)\s+([\+\-])\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([,\d]+)\s+([,\d]+)\s+(\w+\-?\d+?)?\s+(\w+\.\d+)$/) {
    $chromosome  = $2;
    $strand      = $3;
    $id          = $1;

    $transcStart = $4;
    $transcEnd   = $5;


    $codingStart = $6;
    $codingEnd   = $7;

    $exonsStart  = $9;
    $exonsEnd    = $10;
    if($11){
      $proteinID=$11;
    }
    if($id ne $12){
      die "KnownGenesCodingTranscript.pm Parsing problem\n";
    }
  } else {
    die "Database line ".$arg{knowngeneline}." did not parse\n";
  }


  $self = $class->SUPER::new(chromosome   => $chromosome,
			     strand       => $strand,
			     id           => $id,
			     transcStart  => $transcStart,
			     transcEnd    => $transcEnd,
			     codingStart  => $codingStart,
			     codingEnd    => $codingEnd,
			     exonsStart   => $exonsStart,
			     exonsEnd     => $exonsEnd,
			     baseDir      => $arg{baseDir},
			     seqExtractor => $arg{seqExtractor} );

  if($proteinID){
    $self->{_proteinID}=$proteinID;
  }else{
    $self->{_proteinID}=undef;
  }


  return $self;
}




1;
