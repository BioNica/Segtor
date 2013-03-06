=head1 NAME

   GeneDataStructure::CreateGeneStructure.pm


=head1 Disclaimer

This software and all the data therein were developed by the Bioinformatics and Computational Biology Lab (Laboratorio de Bioinformatica e Biologia Computational (LBBC)) of the Brazilian National Cancer Institute (Instituto Nacional de Cancer (INCA)) in Rio de Janeiro, RJ, Brazil. The INCA is part of the Health Ministry (Ministerio da Saude) of the Federal Government of Brazil. This software was developed by the staff of the LBBC as part of their official duties. According to federal regulations, this software is therefore public domain.

This software is provided as is and the INCA does make any guarantee regarding the accuracy/quality/reliability of the results. The INCA cannot be held accountable for any hardware or software damage that might be caused through the use or misuse of this software.

The LBBC encourages researchers using our software to acknowledge us. The software can be distributed freely and modified given that any derivative bears an acknowledgment that it was derived from our original software.

=head1 Author

Gabriel Renaud grenaud@inca.gov.br



=cut

package CreateGeneStructure;

use strict;
use warnings;

use Data::Dumper;



use Cwd qw(abs_path);
my $path;

BEGIN {
  my @pathtemp = split("/", abs_path(__FILE__) );
  delete($pathtemp[$#pathtemp]);
  delete($pathtemp[$#pathtemp]);
  $path=join("/",@pathtemp);
}

use lib $path."/ChromosomeIndexing/";
use ChromosomeIndexing;



=item getArrayOfGenes

  Subroutine to return the array
  of genes

=cut
sub getArrayOfGenes{
  my ($self)=@_;
  return $self->{_arrayOfGenes};
}

=item getArrayOfTranscript

  Subroutine to return the array
  of transcripts

=cut
sub getArrayOfTranscript{
  my ($self)=@_;
  return $self->{_arrayOfTranscripts};
}






=item internalGetArrayOfGenesForGivenChrIndex

   Subroutine to return an array of RefSeqGene.pm
   for a given chromosome index. Do not use externally, this is
   meant as a private subroutine

=cut
sub internalGetArrayOfGenesForGivenChrIndex{
  my ($self,$index)=@_;
  return @{$self->{_arrayOfGenesPerChr}->[$index]};
}

=item getArrayOfGenesForGivenChr

   Subroutine to return an array of RefSeqGene.pm
   for a given chromosome. Please see :
   http://10.46.8.53/wiki/index.php5/ChromosomeIndexing
   for more info on accepted chromosome names.

=cut
sub getArrayOfGenesForGivenChr{
  my ($self,$chr)=@_;
  my $index=ChromosomeIndexing::chr2index($self->{_species}, $chr);
  return internalGetArrayOfGenesForGivenChrIndex($self,$index);
}

=item getArrayOfGenesForGivenIndex

   Subroutine to return an array of RefSeqGene.pm
   for a given chromosome index. Please see :
   http://10.46.8.53/wiki/index.php5/ChromosomeIndexing
   for more info on chromosome indices.

=cut
sub getArrayOfGenesForGivenIndex{
  my ($self,$index)=@_;
  if(!ChromosomeIndexing::chrIndexEXISTS($self->{_species}, $index)){
    die "CreateRefSeqGene.pm : The index $index does not exists\n";
  }
  return internalGetArrayOfGenesForGivenChrIndex($self,$index);
}



=item internalGetArrayOfTranscriptsForGivenChrIndex

   Subroutine to return an array of RefSeqCodingTranscript.pm
   for a given chromosome index. Do not use externally, this is
   meant as a private subroutine

=cut
sub internalGetArrayOfTranscriptsForGivenChrIndex{
  my ($self,$index)=@_;
  if($self->{_arrayOfTranscriptsPerChr}->[$index]){
    return @{$self->{_arrayOfTranscriptsPerChr}->[$index]};
  }else{
    my @arrayTemp=();
    return @arrayTemp;
  }
}

=item getArrayOfTranscriptsForGivenChr

   Subroutine to return an array of RefSeqCodingTranscript.pm
   for a given chromosome. Please see :
   http://10.46.8.53/wiki/index.php5/ChromosomeIndexing
   for more info on accepted chromosome names.

=cut
sub getArrayOfTranscriptsForGivenChr{
  my ($self,$chr)=@_;
  my $index=ChromosomeIndexing::chr2index($self->{_species}, $chr);
  return internalGetArrayOfTranscriptsForGivenChrIndex($self,$index);
}

=item internalGetArrayOfTranscriptsForGivenChrIndex

   Subroutine to return an array of RefSeqCodingTranscript.pm
   for a given chromosome index. Please see :
   http://10.46.8.53/wiki/index.php5/ChromosomeIndexing
   for more info on chromosome indices.

=cut
sub getArrayOfTranscriptsForGivenIndex{
  my ($self,$index)=@_;
  if(!ChromosomeIndexing::chrIndexEXISTS($self->{_species}, $index)){
    die "CreateRefSeqGene.pm : The index $index does not exists\n";
  }
  return internalGetArrayOfTranscriptsForGivenChrIndex($self,$index);
}





1;
