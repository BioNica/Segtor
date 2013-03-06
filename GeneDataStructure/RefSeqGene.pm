=head1 NAME

   GeneDataStructure::RefSeqGene.pm


=head1 Disclaimer

This software and all the data therein were developed by the Bioinformatics and Computational Biology Lab (Laboratorio de Bioinformatica e Biologia Computational (LBBC)) of the Brazilian National Cancer Institute (Instituto Nacional de Cancer (INCA)) in Rio de Janeiro, RJ, Brazil. The INCA is part of the Health Ministry (Ministerio da Saude) of the Federal Government of Brazil. This software was developed by the staff of the LBBC as part of their official duties. According to federal regulations, this software is therefore public domain.

This software is provided as is and the INCA does make any guarantee regarding the accuracy/quality/reliability of the results. The INCA cannot be held accountable for any hardware or software damage that might be caused through the use or misuse of this software.

The LBBC encourages researchers using our software to acknowledge us. The software can be distributed freely and modified given that any derivative bears an acknowledgment that it was derived from our original software.

=head1 Author

Gabriel Renaud grenaud@inca.gov.br



=cut

package RefSeqGene;


use strict;
use warnings;

use Data::Dumper;
use Gene;
our @ISA = "Gene";



=item determineClassName

  Subroutine to determine the class
  name to make sure the transcripts are
  instances of EnsemblCodingTranscript

=cut
sub determineClassName{
  my ($class)=@_;
  my $stringToParse=Dumper($class);
  if($stringToParse =~ /,\s\'(\w+)\'\s\);$/){
    return $1;
  }else{
    die "determineClassName() Big problem\n";
  }
}


=item new

   Constructor for RefSeqGene.pm. This object is used to
   store a single RefGene gene. See CreateRefSeqGenes.pm
   as to how to build the data structure.

   This modules needs the following additional arguments:
   'genesymbol' and 'locuslinkid'


=cut
sub new{
  my ($class,%arg)=(@_);


  foreach my $transcript (@{$arg{transcriptArray}}){
    if(determineClassName($transcript) ne "RefSeqCodingTranscript"){
      die "RefSeqGene.pm, the transcriptArray must be of the type RefSeqCodingTranscript\n";
    }
  }
  my $self = $class->SUPER::new(%arg);

  if( !(exists $arg{genesymbol}) ){
    die "RefSeqGene.pm : The 'genesymbol' parameter is mandatory\n";
  }

  if( !(exists $arg{locuslinkid}) ){
    die "RefSeqGene.pm : The 'locuslinkid' parameter is mandatory\n";
  }


  $self->{_genesymbol}=$arg{genesymbol};
  $self->{_locuslinkid}=$arg{locuslinkid};

  return $self;
}

=item getGeneSymbol

  Subroutine to return the
  gene symbol

=cut
sub getGeneSymbol{
  my($self)=@_;
  return $self->{_genesymbol};
}


=item getLocusLinkID

  Subroutine to return the
  locus link id

=cut
sub getLocusLinkID{
  my($self)=@_;
  return $self->{_locuslinkid};
}


=item print

  Subroutine to return the
  a string representation of the
  RefSeqGene.pm object

=cut
sub print{
  my($self)=@_;
  my $stringToReturn="";
  $stringToReturn.="Gene ".$self->{_genesymbol}." locuslinkID = ".$self->{_locuslinkid}."\n";
  $stringToReturn.=$self->{_chromosome}.":".$self->{_strand}.":".$self->{_farthest5PrimeEnd}."-".$self->{_farthest3PrimeEnd}."\n";
  $stringToReturn.=$self->{_numberOfTranscript}." variants :\n";

  for(my $transcriptIndex=0;$transcriptIndex<=$#{$self->{_transcriptsArray}};$transcriptIndex++){
    $stringToReturn.=$self->{_transcriptsArray}->[$transcriptIndex]->print();
  }
  return $stringToReturn;
}



1;
