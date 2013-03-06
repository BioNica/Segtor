=head1 NAME

   GeneDataStructure::RefSeqCodingTranscript.pm


=head1 Disclaimer

This software and all the data therein were developed by the Bioinformatics and Computational Biology Lab (Laboratorio de Bioinformatica e Biologia Computational (LBBC)) of the Brazilian National Cancer Institute (Instituto Nacional de Cancer (INCA)) in Rio de Janeiro, RJ, Brazil. The INCA is part of the Health Ministry (Ministerio da Saude) of the Federal Government of Brazil. This software was developed by the staff of the LBBC as part of their official duties. According to federal regulations, this software is therefore public domain.

This software is provided as is and the INCA does make any guarantee regarding the accuracy/quality/reliability of the results. The INCA cannot be held accountable for any hardware or software damage that might be caused through the use or misuse of this software.

The LBBC encourages researchers using our software to acknowledge us. The software can be distributed freely and modified given that any derivative bears an acknowledgment that it was derived from our original software.

=head1 Author

Gabriel Renaud grenaud@inca.gov.br



=cut

package RefSeqCodingTranscript;

use strict;
use warnings;

use Data::Dumper;
use Range;
use CodingTranscript;
our @ISA = "CodingTranscript";

=item new

   Constructor for the RefSeqCodingTranscript object. This object inherits from
   CodingTranscript.pm and is used to store information pertaining to an RefGene Transcript.
   This object can be constructed by providing a line from the refGeneline.txt and
   the corresponding line from refSeqAliline.txt table from UCSC.

=cut

sub new{
  my ($class,%arg)=(@_);

  if( !(exists $arg{refGeneline}) ){
    die "RefSeqCodingTranscript.pm Enter the 'refGeneline' parameter\n";
  }

  if( !(exists $arg{refSeqAliline}) ){
    die "RefSeqCodingTranscript.pm Enter the 'refSeqAliline' parameter\n";
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
  my $querySize;
  my $queryStart;
  my $queryEnd;

  my $framesFromFile;

  #585     NM_207180       chrE64_random   -       113563  120203  113805  120203  7       113563,115464,116307,116970,117758,118919,120172,   113840,115794,116637,117300,118088,118946,120203,       4909    PIT 54  cmpl    cmpl    1,1,1,1,1,1,0,

  #1       2               3       4       5               6               7               8               9       10                                         11                                       12      13              14      15      16
  #1643    NM_016459       chr5    -       138751155       138753504       138751352       138753444       4       138751155,138751608,138752048,138753267,   138751509,138751719,138752173,138753504, 0       MGC29506        cmpl    cmpl    2,2,0,0,

  #                             1       2        3           4         5       6       7       8       9       10         11        12      13         14      15       16

  if ( $arg{refGeneline} =~ /^(\d+)\s+(\w+)\s+([\.\w]+)\s+([\+\-])\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([,\d]+)\s+([,\d]+)\s+(\d+)\s+([ \S]+)\s+(\w+)\s+(\w+)\s+([\-,012]+)$/ ) {

    $chromosome  = $3;
    $strand      = $4;
    $id          = $2;

    $transcStart = $5;
    $transcEnd   = $6;

    $codingStart = $7;
    $codingEnd   = $8;

    $exonsStart  = $10;
    $exonsEnd    = $11;

    $framesFromFile=$16;

  } else {
    die "RefSeqCodingTranscript.pm Database refgene line ".$arg{refGeneline}." did not parse\n";
  }

  #1       2       3       4       5       6       7       8       9      10       11              12     13       14      15      16              17         18               19      20                      21              22
  #1643    796     0       31      0       0       0       3       1522    -       NM_016459       845     0       827     chr5    180857866       138751155  138753504        4       354,111,125,237,        18,372,483,608, 138751155,138751608,138752048,138753267,

  #                               1       2       3       4       5       6       7       8       9       10        11      12      13      14      15      16      17      18       19      20         21          22
  if ( $arg{refSeqAliline} =~ /^(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([\+\-])\s+(\w+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([\.\w]+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([,\d]+)\s+([,\d]+)\s+([,\d]+)$/ ) {

    if($id ne $11){
      die "RefSeqCodingTranscript.pm The refSeq ID $id found in refGene does not match the one found in refSeqAli: $11\n";
    }

    if($chromosome ne $15){
      die "RefSeqCodingTranscript.pm The chromosome $chromosome found in refGene does not match the one found in refSeqAli: $15\n";
    }

    if($strand ne $10){
      die "RefSeqCodingTranscript.pm The strand $strand found in refGene does not match the one found in refSeqAli: $10\n";
    }

    if($transcStart != $17){
      die "RefSeqCodingTranscript.pm The transcription start $transcStart found in refGene does not match the one found in refSeqAli: $17\n";
    }

    if($transcEnd != $18){
      die "RefSeqCodingTranscript.pm The transcription end $transcEnd found in refGene does not match the one found in refSeqAli: $18\n";
    }

    if($exonsStart ne $22){
      #warn "RefSeqCodingTranscript.pm For ID = $id, the exons start $exonsStart found in refGene does not match the one found in refSeqAli: $22\n";
      #warn "RefSeqCodingTranscript.pm for ID = $id, the exons start found in refGene does not match the one found in refSeqAli\n";
    }

    $querySize  = $12;
    $queryStart = $13;
    $queryEnd   = $14;

  }else{
    die "RefSeqCodingTranscript.pm Database refseqali line ".$arg{refSeqAliline}." did not parse\n";
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
			     querySize      => $querySize,
			     queryStart     => $queryStart,
			     queryEnd       => $queryEnd,
			     baseDir        => $arg{baseDir},
			     seqExtractor   => $arg{seqExtractor},
			     framesFromFile => $framesFromFile);




  return $self;
}




1;
