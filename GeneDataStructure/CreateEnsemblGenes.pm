=head1 NAME

   GeneDataStructure::CreateEnsemblGenes.pm


=head1 Disclaimer

This software and all the data therein were developed by the Bioinformatics and Computational Biology Lab (Laboratorio de Bioinformatica e Biologia Computational (LBBC)) of the Brazilian National Cancer Institute (Instituto Nacional de Cancer (INCA)) in Rio de Janeiro, RJ, Brazil. The INCA is part of the Health Ministry (Ministerio da Saude) of the Federal Government of Brazil. This software was developed by the staff of the LBBC as part of their official duties. According to federal regulations, this software is therefore public domain.

This software is provided as is and the INCA does make any guarantee regarding the accuracy/quality/reliability of the results. The INCA cannot be held accountable for any hardware or software damage that might be caused through the use or misuse of this software.

The LBBC encourages researchers using our software to acknowledge us. The software can be distributed freely and modified given that any derivative bears an acknowledgment that it was derived from our original software.

=head1 Author

Gabriel Renaud grenaud@inca.gov.br



=cut


package CreateEnsemblGenes;


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


use EnsemblCodingTranscript;
use EnsemblGene;

use lib $path."/ChromosomeIndexing/";
use ChromosomeIndexing;


use CreateGeneStructure;
our @ISA = "CreateGeneStructure";



=item new

   Constructor to build the array of EnsemblCodingTranscript.pm
   and the array of EnsemblGenes.pm
   using the species name. For more information on accepted
   species name see:
   http://10.46.8.53/wiki/index.php5/GenomeBuildNaming

   The modules stores the transcripts in 2 arrays:
   _arrayOfTranscripts         : all of the EnsemblCodingTranscript.pm
   _arrayOfTranscriptsPerChr   : The EnsemblCodingTranscript.pm on a given chromosome

   The modules stores the genes in 2 arrays:
   _arrayOfGenes         : all of the EnsemblGenes.pm
   _arrayOfGenesPerChr   : The EnsemblGenes.pm on a given chromosome

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


  if(!(-e ChromosomeIndexing::getBaseDir().$self->{_species}."/database/ensGene.txt")){
    die "The Ensembl file ".ChromosomeIndexing::getBaseDir().$self->{_species}."/database/ensGene.txt does not exists\n";
  }

  $self->{_arrayOfTranscripts}=[];
  $self->{_arrayOfGenes}=[];
  $self->{_arrayOfGenesPerChr}=[];
  $self->{_arrayOfTranscriptsPerChr}=[];


  foreach my $indexChr  (ChromosomeIndexing::listIndexChr($self->{_species})){
    $self->{_arrayOfGenesPerChr}->[$indexChr]=[];
    $self->{_arrayOfTranscriptsPerChr}->[$indexChr]=[];
  }

  my %hashEnsGene2ArrayOfLines;
  my $file=ChromosomeIndexing::getBaseDir().$self->{_species}."/database/ensGene.txt";
  open(INFO, $file) or die ("Can't open $file ");
  while(my $line = <INFO>){
    my @arrayTemp=split("\t",$line);
    my $key=$arrayTemp[2] ."@==@".$arrayTemp[3] ."@==@". $arrayTemp[12];
    if(exists $hashEnsGene2ArrayOfLines{$key}){
      push(@{$hashEnsGene2ArrayOfLines{$key}}, $line);
    }else{
      $hashEnsGene2ArrayOfLines{$key}=[$line];
    }
  }
  close(INFO);


  foreach my $geneid (keys %hashEnsGene2ArrayOfLines){
    my $arrayOfTranscript=[];
    foreach my $lineToAdd (@{$hashEnsGene2ArrayOfLines{$geneid}}){
      my $trans;
      if($self->{_baseDir}){
	$trans=EnsemblCodingTranscript->new(ensemblline   => $lineToAdd,
					    baseDir       => $self->{_baseDir}.$self->{_species}."/chromosomes/");
      }else{
	$trans=EnsemblCodingTranscript->new(ensemblline   => $lineToAdd);
      }

      push(@{$self->{_arrayOfTranscripts}},$trans);
      push(@{$self->{_arrayOfTranscriptsPerChr}->[ChromosomeIndexing::chr2index($self->{_species}, $trans->getChromosome()) ]}  ,$trans);
      push(@{$arrayOfTranscript},$trans);
    }
    my $geneSymbolToSend;
    if($geneid =~ /^(.*)@==@(.*)@==@(.*)/){
      $geneSymbolToSend=$3;
    }else{
      die "CreateEnsemblGenes.pm, wrong key $geneid for hash ".Dumper($hashEnsGene2ArrayOfLines{$geneid});
    }

    my $gene=EnsemblGene->new(id                => $geneSymbolToSend,
			      transcriptArray   => $arrayOfTranscript);
    push(@{$self->{_arrayOfGenes}},$gene);
    push(@{$self->{_arrayOfGenesPerChr}->[ChromosomeIndexing::chr2index($self->{_species}, $gene->getChromosome()) ]}  ,$gene);
  }

  return $self;
}










1;

