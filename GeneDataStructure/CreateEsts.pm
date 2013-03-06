=head1 NAME

   GeneDataStructure::CreateEsts.pm


=head1 Disclaimer

This software and all the data therein were developed by the Bioinformatics and Computational Biology Lab (Laboratorio de Bioinformatica e Biologia Computational (LBBC)) of the Brazilian National Cancer Institute (Instituto Nacional de Cancer (INCA)) in Rio de Janeiro, RJ, Brazil. The INCA is part of the Health Ministry (Ministerio da Saude) of the Federal Government of Brazil. This software was developed by the staff of the LBBC as part of their official duties. According to federal regulations, this software is therefore public domain.

This software is provided as is and the INCA does make any guarantee regarding the accuracy/quality/reliability of the results. The INCA cannot be held accountable for any hardware or software damage that might be caused through the use or misuse of this software.

The LBBC encourages researchers using our software to acknowledge us. The software can be distributed freely and modified given that any derivative bears an acknowledgment that it was derived from our original software.

=head1 Author

Gabriel Renaud grenaud@inca.gov.br



=cut

package CreateEsts;


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

use EstTranscript;
use lib $path."/ChromosomeIndexing/";
use ChromosomeIndexing;


use CreateGeneStructure;
our @ISA = "CreateGeneStructure";




=item new

   Constructor to build the array of EstTranscript.pm
   using the species name. For more information on accepted
   species name see:
   http://10.46.8.53/wiki/index.php5/GenomeBuildNaming

   The modules stores the transcripts in 2 arrays:
   _arrayOfTranscripts         : all of the EstTranscript.pm
   _arrayOfTranscriptsPerChr   : The EstTranscript.pm on a given chromosome

=cut

sub new{
  my ($class,%arg)=(@_);
  my ($self) =bless{_species                   => $arg{species},
		    _baseDir                   => $arg{baseDir}}, $class;

  if($self->{_baseDir}){
    ChromosomeIndexing::setBaseDir($self->{_baseDir});
  }

  if(! ChromosomeIndexing::indexExistsForSpecies($self->{_species})){
    die "The species ".$self->{_species}." was not indexed\n";
  }

  if(!(-e ChromosomeIndexing::getBaseDir().$self->{_species}."/database/all_est.txt")){
    die "The Est file ".ChromosomeIndexing::getBaseDir().$self->{_species}."/database/all_ests.txt does not exists\n";
  }

  $self->{_arrayOfTranscripts}=[];
  $self->{_arrayOfTranscriptsPerChr}=[];
  my $file=ChromosomeIndexing::getBaseDir().$self->{_species}."/database/all_est.txt";

  open(INFO, $file) or die ("Can't open $file ");
  while (my $line = <INFO>) {
    my $trans;
    if($self->{_baseDir}){
      $trans=EstTranscript->new(databaseline    => $line,
				baseDir         => $self->{_baseDir}.$self->{_species}."/chromosomes/");
    }else{
      $trans=EstTranscript->new(databaseline   => $line );
    }
    #print Dumper($trans);
    push(@{$self->{_arrayOfTranscripts}},$trans);
    push(@{$self->{_arrayOfTranscriptsPerChr}->[ChromosomeIndexing::chr2index($self->{_species}, $trans->getChromosome()) ]}  ,$trans);
  }
  close(INFO);

  return $self;
}




#
#
#=item internalGetArrayOfTranscriptsForGivenChrIndex
#
#   Subroutine to return an array of EstTranscript.pm
#   for a given chromosome index. Do not use externally, this is
#   meant as a private subroutine
#
#=cut
#sub internalGetArrayOfTranscriptsForGivenChrIndex{
#  my ($self,$index)=@_;
#  return @{$self->{_arrayOfTranscriptsPerChr}->[$index]};
#}
#
#=item getArrayOfTranscriptsForGivenChr
#
#   Subroutine to return an array of EstTranscript.pm
#   for a given chromosome. Please see :
#   http://10.46.8.53/wiki/index.php5/ChromosomeIndexing
#   for more info on accepted chromosome names.
#
#=cut
#sub getArrayOfTranscriptsForGivenChr{
#  my ($self,$chr)=@_;
#  my $index=ChromosomeIndexing::chr2index($self->{_species}, $chr);
#  return internalGetArrayOfTranscriptsForGivenChrIndex($index);
#}
#
#=item internalGetArrayOfTranscriptsForGivenChrIndex
#
#   Subroutine to return an array of EstTranscript.pm
#   for a given chromosome index. Please see :
#   http://10.46.8.53/wiki/index.php5/ChromosomeIndexing
#   for more info on chromosome indices.
#
#=cut
#sub getArrayOfTranscriptsForGivenIndex{
#  my ($self,$index)=@_;
#  if(!ChromosomeIndexing::chrIndexEXISTS($self->{_species}, $index)){
#    die "CreateEsts.pm : The index $index does not exists\n";
#  }
#  return internalGetArrayOfTranscriptsForGivenChrIndex($index);
#}
#
#
#
1;


