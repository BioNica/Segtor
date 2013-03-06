=head1 NAME

   GeneDataStructure::CreateRefSeqGenes.pm


=head1 Disclaimer

This software and all the data therein were developed by the Bioinformatics and Computational Biology Lab (Laboratorio de Bioinformatica e Biologia Computational (LBBC)) of the Brazilian National Cancer Institute (Instituto Nacional de Cancer (INCA)) in Rio de Janeiro, RJ, Brazil. The INCA is part of the Health Ministry (Ministerio da Saude) of the Federal Government of Brazil. This software was developed by the staff of the LBBC as part of their official duties. According to federal regulations, this software is therefore public domain.

This software is provided as is and the INCA does make any guarantee regarding the accuracy/quality/reliability of the results. The INCA cannot be held accountable for any hardware or software damage that might be caused through the use or misuse of this software.

The LBBC encourages researchers using our software to acknowledge us. The software can be distributed freely and modified given that any derivative bears an acknowledgment that it was derived from our original software.

=head1 Author

Gabriel Renaud grenaud@inca.gov.br



=cut

package CreateRefSeqGenes;


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

use RefSeqCodingTranscript;
use RefSeqGene;
use lib $path."/ChromosomeIndexing/";

use ChromosomeIndexing;
use CreateGeneStructure;
our @ISA = "CreateGeneStructure";

my $skipNonNms=1;



=item clusterTranscripts

    This subroutine is used to cluster the transcripts
    that overlap each other on the same strand and same chromosome
    The subroutine returns an array of array where each bin contains
    an array of transcripts that overlap each other

=cut
sub clusterTranscripts{
  my ($arrayOfTranscriptsParam)=@_;


  my $index1=0;
  my $index2=$index1+1;
  my $keepLooping=1;
  my $transcripts=[];
  my $initialArraySize=$#{$arrayOfTranscriptsParam}+1;
  my $initialNumberOfElements=0;

  #Copying the number of elements
  for(my $index1=0;$index1<=$#{$arrayOfTranscriptsParam};$index1++){
    $transcripts->[$index1]  =$arrayOfTranscriptsParam->[$index1];
    $initialNumberOfElements+=($#{$arrayOfTranscriptsParam->[$index1]}+1);
  }


  while($keepLooping){ #we will keep until no overlaps are found
    my $transcriptsNextLoop=[];
    my $overlapFoundForCurrentIteration=0;

    #For each element of the array, we will check if it overlaps
    #another record in the array.
    for(my $index1=0;$index1<=$#{$transcripts};$index1++){
      my $hasHadAnOverlap=0;
      my $indexOverlap=0;

      #Checking for the remaining elements starting at element at index1+1
      for(my $index2=($index1+1);$index2<=$#{$transcripts};$index2++){

	#test if transcripts in index1 and transcripts in index2 overlap
	#Transcript.pm checks to see if the transcripts are the same chromosomes
	foreach my $transcript1 (@{ $transcripts->[$index1]  }) {
	  foreach my $transcript2 (@{ $transcripts->[$index2]  }) {
	    if ( $transcript1->overlapsWithStrand($transcript2)  ) {
	      $hasHadAnOverlap=1;
	    }
	  }
	}

	# we found an overlap, we break out of the loop and store
	# the index of the element where the overlap was found
	if($hasHadAnOverlap){
	  $indexOverlap=$index2;
	  last;
	}
      }#end loop index2

      if($hasHadAnOverlap){#an overlap was found between index1 and index2,
	#we copy the rest of the elements
	for(my $index2=($index1);$index2<=$#{$transcripts};$index2++){
	  if($index2<$indexOverlap){ #we store the elements that occur prior to the overlapping one
	    $transcriptsNextLoop->[$index2]=$transcripts->[$index2];
	  }elsif($index2>$indexOverlap){ #we store the elements that occur after to the overlapping one at an the same index minus 1
	    $transcriptsNextLoop->[$index2-1]=$transcripts->[$index2];
	  }elsif($index2 == $indexOverlap){ #we store the element that overlapped with index1 with index1 for the next iteration;
	    push(@{$transcriptsNextLoop->[$index1]},@{$transcripts->[$index2]});
	  }else{
	    die "Big trouble\n";
	  }
	}#end copying for
	$overlapFoundForCurrentIteration=1;

	last;
      }else{#If not overlap was found, we store the element as is.
	$transcriptsNextLoop->[$index1]=$transcripts->[$index1];
      }


    }#end loop index1

    if($overlapFoundForCurrentIteration){ #if overlap was found
      $transcripts=$transcriptsNextLoop;
    }else{
      $keepLooping=0;
    }

  }#end while(keepLooping)

  my $finalArraySize=$#{$transcripts}+1;
  my $finalNumberOfElements=0;

  #Copying the number of elements
  for(my $index1=0;$index1<=$#{$transcripts};$index1++){
    $finalNumberOfElements+=($#{$transcripts->[$index1]}+1);
  }

  if($initialNumberOfElements != $finalNumberOfElements){
    die "CreateRefSeqGenes.pm clusterTranscripts failed\n";
  }

  return $transcripts;
}


=item new

   Constructor to build the array of RefSeqCodingTranscript.pm
   and the array of RefSeqGenes.pm
   using the species name. For more information on accepted
   species name see:
   http://10.46.8.53/wiki/index.php5/GenomeBuildNaming

   The modules stores the transcripts in 2 arrays:
   _arrayOfTranscripts         : all of the RefSeqCodingTranscript.pm
   _arrayOfTranscriptsPerChr   : The RefSeqCodingTranscript.pm on a given chromosome

   The modules stores the genes in 2 arrays:
   _arrayOfGenes         : all of the RefSeqGenes.pm
   _arrayOfGenesPerChr   : The RefSeqGenes.pm on a given chromosome

=cut
sub new{
  my ($class,%arg)=(@_);
  my ($self) =bless{_species                   => $arg{species},
		    _baseDir                   => $arg{baseDir} }, $class;

  if($self->{_baseDir}){
    ChromosomeIndexing::setBaseDir($self->{_baseDir});
  }

  if(! ChromosomeIndexing::indexExistsForSpecies($self->{_species})){
    print "The species ".$self->{_species}." was not indexed\n";
    die;
  }


  if(!(-e ChromosomeIndexing::getBaseDir().$self->{_species}."/database/refGene.txt")){
    print "The RefSeq file ".ChromosomeIndexing::getBaseDir().$self->{_species}."/database/refGene.txt does not exists\n";
    die;
  }
  if(!(-e ChromosomeIndexing::getBaseDir().$self->{_species}."/database/refFlat.txt")){
    print "The RefSeq file ".ChromosomeIndexing::getBaseDir().$self->{_species}."/database/refFlat.txt does not exists\n";
    die;
  }
  if(!(-e ChromosomeIndexing::getBaseDir().$self->{_species}."/database/refLink.txt")){
    print "The RefSeq file ".ChromosomeIndexing::getBaseDir().$self->{_species}."/database/refLink.txt does not exists\n";
    die;
  }
  if(!(-e ChromosomeIndexing::getBaseDir().$self->{_species}."/database/refSeqAli.txt")){
    print "The RefSeq file ".ChromosomeIndexing::getBaseDir().$self->{_species}."/database/refSeqAli.txt does not exists\n";
    die;
  }



  $self->{_arrayOfTranscripts}=[];
  $self->{_arrayOfGenes}=[];
  $self->{_arrayOfGenesPerChr}=[];
  $self->{_arrayOfTranscriptsPerChr}=[];

  foreach my $indexChr  (ChromosomeIndexing::listIndexChr($self->{_species})){
    $self->{_arrayOfGenesPerChr}->[$indexChr]=[];
    $self->{_arrayOfTranscriptsPerChr}->[$indexChr]=[];
  }



  my $file;

  my %hashNM2LocusLink;
  $file=ChromosomeIndexing::getBaseDir().$self->{_species}."/database/refLink.txt";
  open(INFO, $file) or die ("Can't open $file ");
  while(my $line = <INFO>){
    my @arrayTemp=split("\t",$line);
    if($skipNonNms && $arrayTemp[2] !~ /^NM/){
      next;
    }


    if(exists $hashNM2LocusLink{$arrayTemp[2]}){
      print "CreateRefSeqGenes.pm Big problem #1, id ".$arrayTemp[2]." found twice\n";
      die;
    }else{
      $hashNM2LocusLink{$arrayTemp[2]}= $arrayTemp[6] ;
    }
  }
  close(INFO);



  my %hashNM2GeneSymbol;
  $file=ChromosomeIndexing::getBaseDir().$self->{_species}."/database/refFlat.txt";
  open(INFO, $file) or die ("Can't open $file ");
  while(my $line = <INFO>){
    my @arrayTemp=split("\t",$line);

    if($skipNonNms && $arrayTemp[1] !~ /^NM/){
      next;
    }

    if(exists $hashNM2GeneSymbol{$arrayTemp[1]}){
      if($hashNM2GeneSymbol{$arrayTemp[1]} ne $arrayTemp[0] ){
	print "CreateRefSeqGenes.pm Big problem #2, id ".$arrayTemp[1]." found twice with different symbols\n";
	die;
      }
    }else{
      $hashNM2GeneSymbol{$arrayTemp[1]}= $arrayTemp[0] ;
    }
  }
  close(INFO);




  my %hashLL2RefGeneLine;
  my %hashLL2RefSeqAliLine;
  my %hashLL2Key;

  my %hashNM2RefGeneLine;
  $file=ChromosomeIndexing::getBaseDir().$self->{_species}."/database/refGene.txt";

  open(INFO, $file) or die ("Can't open $file ");
  while(my $line = <INFO>){

    my @arrayTemp=split("\t",$line);

    if($skipNonNms && $arrayTemp[1] !~ /^NM/){
      next;
    }

    if(!(exists $hashNM2LocusLink{$arrayTemp[1]})){
      print $arrayTemp[1]." not found\n";
      die;
    }

    my $locusLink=$hashNM2LocusLink{$arrayTemp[1]};
    if(exists $hashLL2RefGeneLine{$locusLink}){
      push( @{$hashLL2RefGeneLine{$locusLink}} ,$line );
    }else{
      $hashLL2RefGeneLine{$locusLink}  = [$line];
    }



    #       bin           NM                chr               start              end
    my $key=$arrayTemp[0].$arrayTemp[1]."#".$arrayTemp[2]."#".$arrayTemp[4]."-".$arrayTemp[5];

    if(exists $hashNM2RefGeneLine{$key}){
      warn "CreateRefSeqGenes.pm Big problem #3, id ".$key." found twice\n";
    }else{
      $hashNM2RefGeneLine{$key}= $line;

      if (exists $hashLL2Key{$locusLink}) {
	push( @{$hashLL2Key{$locusLink}} ,$key );
      } else {
	$hashLL2Key{$locusLink}  = [$key];
      }
    }


  }
  close(INFO);



  my %hashNM2RefSeqAliLine;
  $file=ChromosomeIndexing::getBaseDir().$self->{_species}."/database/refSeqAli.txt";
  open(INFO, $file) or die ("Can't open $file ");
  while(my $line = <INFO>){
    my @arrayTemp=split("\t",$line);

    if($skipNonNms && $arrayTemp[10] !~ /^NM/){
      next;
    }

    if(!(exists $hashNM2LocusLink{$arrayTemp[10]})){
      print $arrayTemp[10]." not found\n";
      die;
    }
    my $locusLink=$hashNM2LocusLink{$arrayTemp[10]};
    if(exists $hashLL2RefSeqAliLine{$locusLink}){
      push( @{$hashLL2RefSeqAliLine{$locusLink}} ,$line );
    }else{
      $hashLL2RefSeqAliLine{$locusLink}  = [$line];
    }

    #       bin           NM                chr                 start              end
    my $key=$arrayTemp[0].$arrayTemp[10]."#".$arrayTemp[14]."#".$arrayTemp[16]."-".$arrayTemp[17];

    chomp($key);

    if(exists $hashNM2RefSeqAliLine{ $key }){
      warn "CreateRefSeqGenes.pm Big problem #4, id ".$key." found twice\n";
    }else{
      $hashNM2RefSeqAliLine{ $key }= $line;
    }
  }
  close(INFO);




  foreach my $locusLink (keys %hashLL2Key) {


    my $arrayOfTranscripts;

    foreach my $key (@{$hashLL2Key{$locusLink}}) {
      if (! (exists $hashNM2RefGeneLine{$key}) ) {
	warn "CreateRefSeqGenes.pm Big problem #6 id ".$key." not found, skipping\n";
	next;
      }
      if (! (exists $hashNM2RefSeqAliLine{$key}) ) {
	warn "CreateRefSeqGenes.pm Big problem #7 id ".$key." not found, skipping\n";
	next;
      }
      my $lineRefGene   = $hashNM2RefGeneLine{$key};
      my $lineRefSeqAli = $hashNM2RefSeqAliLine{$key};

      my $trans;
      if($self->{_baseDir}){
	$trans=RefSeqCodingTranscript->new(refGeneline   => $lineRefGene,
					   refSeqAliline => $lineRefSeqAli,
					   baseDir       => $self->{_baseDir}.$self->{_species}."/chromosomes/");
      }else{
	$trans=RefSeqCodingTranscript->new(refGeneline   => $lineRefGene,
					   refSeqAliline => $lineRefSeqAli);
      }

      $trans->setGeneSymbol( $hashNM2GeneSymbol{ $trans->getId() } );

      push(@{$self->{_arrayOfTranscripts}},$trans);
      push(@{$self->{_arrayOfTranscriptsPerChr}->[ChromosomeIndexing::chr2index($self->{_species}, $trans->getChromosome()) ]}  ,$trans);

      push(@{$arrayOfTranscripts},[$trans]);
    }


    #This is done to cluster the transcripts that overlap each other
    #on the same chromosome
    $arrayOfTranscripts=clusterTranscripts($arrayOfTranscripts);



    my $geneCounterForID=0;
    foreach my $arrayOfClusters (@{$arrayOfTranscripts}){
      my $geneSymbol;
      my $arrayOfClustersToInsert=[];

      if(exists $hashNM2GeneSymbol{$arrayOfClusters->[0]->getId()}){
	$geneSymbol=$hashNM2GeneSymbol{$arrayOfClusters->[0]->getId()};
      }else{
	print "Key ".$arrayOfClusters->[0]->getId()." not found in hash hashNM2GeneSymbol\n";
	die;
      }
      push(@{$arrayOfClustersToInsert},$arrayOfClusters->[0]);



      for(my $index=1;$index<=$#{$arrayOfClusters};$index++){
	if(!(exists $hashNM2GeneSymbol{$arrayOfClusters->[$index]->getId()})){
	  print "Gene symbol was not found for ".$arrayOfClusters->[$index]->getId()."\n";
	  die;
	}


	if($geneSymbol ne $hashNM2GeneSymbol{$arrayOfClusters->[$index]->getId()}){
	  warn "Gene symbol $geneSymbol found in ".$arrayOfClusters->[0]->getId()." is different from ".$arrayOfClusters->[$index]->getId()." that was found for locus link id = $locusLink\n";
	}
	push(@{$arrayOfClustersToInsert},$arrayOfClusters->[$index]);

      }

      my $gene=RefSeqGene->new(id                => $geneSymbol, #$geneCounterForID,
			       transcriptArray   => $arrayOfClustersToInsert,
			       locuslinkid       => $locusLink,
			       genesymbol        => $geneSymbol);

      push(@{$self->{_arrayOfGenes}},$gene);
      push(@{$self->{_arrayOfGenesPerChr}->[ChromosomeIndexing::chr2index($self->{_species}, $gene->getChromosome()) ]}  ,$gene);

      $geneCounterForID++;

    }



  }#for each locus link


  return $self;
}







1;

