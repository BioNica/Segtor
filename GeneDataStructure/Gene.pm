=head1 NAME

   GeneDataStructure::Gene.pm


=head1 Disclaimer

This software and all the data therein were developed by the Bioinformatics and Computational Biology Lab (Laboratorio de Bioinformatica e Biologia Computational (LBBC)) of the Brazilian National Cancer Institute (Instituto Nacional de Cancer (INCA)) in Rio de Janeiro, RJ, Brazil. The INCA is part of the Health Ministry (Ministerio da Saude) of the Federal Government of Brazil. This software was developed by the staff of the LBBC as part of their official duties. According to federal regulations, this software is therefore public domain.

This software is provided as is and the INCA does make any guarantee regarding the accuracy/quality/reliability of the results. The INCA cannot be held accountable for any hardware or software damage that might be caused through the use or misuse of this software.

The LBBC encourages researchers using our software to acknowledge us. The software can be distributed freely and modified given that any derivative bears an acknowledgment that it was derived from our original software.

=head1 Author

Gabriel Renaud grenaud@inca.gov.br



=cut

package Gene;

use strict;
use warnings;

use Data::Dumper;
use Range;


use Cwd qw(abs_path);
my $path;

BEGIN {
  my @pathtemp = split("/", abs_path(__FILE__) );
  delete($pathtemp[$#pathtemp]);
  delete($pathtemp[$#pathtemp]);
  $path=join("/",@pathtemp);
}


use lib $path."/SegmentTree/";
use SegmentTreeRange;
use SegmentTree;

 use vars qw($geneCounter);


=item new

   Constructor for the Gene object. This module
   holds data pertaining to a gene and an array of the various
   transcripts associated to it.

   The object can be created using an ID and
   an array of transcripts as such :

   my $trans1=Transcript->new(chromosome   => "chr20",
   			      strand       => "+",
   			      id           => "NM_001011546",
   			      transcStart  => 17498598,
   			      transcEnd    => 17536652,
   			      exonsStart   => "17498598,17520576,17529382,17533199,17535681,",
   			      exonsEnd     => "17498856,17520708,17529690,17533276,17536652,");

   my $trans2=Transcript->new(chromosome   => "chr20",
   			      strand       => "+",
   			      id           => "NM_006870",
   			      transcStart  => 17498598,
   			      transcEnd    => 17536652,
   			      exonsStart   => "17498598,17529382,17533199,17535681,",
   			      exonsEnd     => "17498856,17529690,17533276,17536652,");


   my $gene=Gene->new(id                => "DSTN",
       		      transcriptArray   => [$trans1,$trans2]);

=cut
sub new{
  my ($class,%arg)=(@_);
  my ($self) =bless{_id                   => $arg{id},
		    _transcriptsArray     => $arg{transcriptArray}}, $class;

  for(my $transcriptIndex=0;$transcriptIndex<=$#{$self->{_transcriptsArray}};$transcriptIndex++){
    $self->{_transcriptsArray}->[$transcriptIndex]->setGene($self);
  }

  $self->{_geneCounter}=$geneCounter++;

  my $chromosome = $self->{_transcriptsArray}->[0]->getChromosome();
  my $strand     = $self->{_transcriptsArray}->[0]->getStrand();

  my $longestTranscript        = $self->{_transcriptsArray}->[0];
  my $longestTranscriptLength  = $self->{_transcriptsArray}->[0]->getTranscriptionRange()->getLength();
  my $shortestTranscript       = $self->{_transcriptsArray}->[0];
  my $shortestTranscriptLength = $self->{_transcriptsArray}->[0]->getTranscriptionRange()->getLength();

  my $farthest5PrimeEnd   = $self->{_transcriptsArray}->[0]->get5Prime();
  my $closest5PrimeEnd    = $self->{_transcriptsArray}->[0]->get5Prime();
  my $farthest3PrimeEnd   = $self->{_transcriptsArray}->[0]->get3Prime();
  my $closest3PrimeEnd    = $self->{_transcriptsArray}->[0]->get3Prime();

  for(my $transcriptIndex=1;$transcriptIndex<=$#{$self->{_transcriptsArray}};$transcriptIndex++){
    if($chromosome ne $self->{_transcriptsArray}->[$transcriptIndex]->getChromosome()){
      die "\n".Dumper($self)."\nGene.pm The chromosomes among the transcripts in 'transcriptArray' has to be the same\n";
    }
    if($strand ne $self->{_transcriptsArray}->[$transcriptIndex]->getStrand()){
      die  "\n".Dumper($self)."Gene.pm The strand among the transcripts in 'transcriptArray' has to be the same\n";
    }

    if($longestTranscriptLength < $self->{_transcriptsArray}->[$transcriptIndex]->getTranscriptionRange()->getLength()){
      $longestTranscriptLength  = $self->{_transcriptsArray}->[$transcriptIndex]->getTranscriptionRange()->getLength();
      $longestTranscript        = $self->{_transcriptsArray}->[$transcriptIndex];
    }

    if($shortestTranscriptLength > $self->{_transcriptsArray}->[$transcriptIndex]->getTranscriptionRange()->getLength()){
      $shortestTranscriptLength  = $self->{_transcriptsArray}->[$transcriptIndex]->getTranscriptionRange()->getLength();
      $shortestTranscript        = $self->{_transcriptsArray}->[$transcriptIndex];
    }


    if($strand eq "+"){
      if( $farthest5PrimeEnd   > $self->{_transcriptsArray}->[$transcriptIndex]->get5Prime() ){
	$farthest5PrimeEnd     = $self->{_transcriptsArray}->[$transcriptIndex]->get5Prime();
      }
      if( $closest5PrimeEnd    < $self->{_transcriptsArray}->[$transcriptIndex]->get5Prime() ){
	$closest5PrimeEnd      = $self->{_transcriptsArray}->[$transcriptIndex]->get5Prime();
      }

      if( $farthest3PrimeEnd   < $self->{_transcriptsArray}->[$transcriptIndex]->get3Prime() ){
	$farthest3PrimeEnd     = $self->{_transcriptsArray}->[$transcriptIndex]->get3Prime();
      }
      if( $closest3PrimeEnd    > $self->{_transcriptsArray}->[$transcriptIndex]->get3Prime() ){
	$closest3PrimeEnd      = $self->{_transcriptsArray}->[$transcriptIndex]->get3Prime();
      }
    }elsif($strand eq "-"){
      if( $farthest5PrimeEnd   < $self->{_transcriptsArray}->[$transcriptIndex]->get5Prime() ){
	$farthest5PrimeEnd     = $self->{_transcriptsArray}->[$transcriptIndex]->get5Prime();
      }
      if( $closest5PrimeEnd    > $self->{_transcriptsArray}->[$transcriptIndex]->get5Prime() ){
	$closest5PrimeEnd      = $self->{_transcriptsArray}->[$transcriptIndex]->get5Prime();
      }

      if( $farthest3PrimeEnd   > $self->{_transcriptsArray}->[$transcriptIndex]->get3Prime() ){
	$farthest3PrimeEnd     = $self->{_transcriptsArray}->[$transcriptIndex]->get3Prime();
      }
      if( $closest3PrimeEnd    < $self->{_transcriptsArray}->[$transcriptIndex]->get3Prime() ){
	$closest3PrimeEnd      = $self->{_transcriptsArray}->[$transcriptIndex]->get3Prime();
      }
    }else{
      die "Big problem\n";
    }
  }



  $self->{_numberOfTranscript}=  ($#{$self->{_transcriptsArray}}+1);

  if( $self->{_numberOfTranscript} > 1 ){
    $self->{_altSplicedTranscript}=1;
  }else{
    $self->{_altSplicedTranscript}=0;
  }

  $self->{_chromosome}=$chromosome;
  $self->{_strand}    =$strand;


  $self->{_longestTranscript}     = $longestTranscript;
  $self->{_shortestTranscript}    = $shortestTranscript;

  $self->{_farthest5PrimeEnd}     = $farthest5PrimeEnd;
  $self->{_closest5PrimeEnd}      = $closest5PrimeEnd;
  $self->{_farthest3PrimeEnd}     = $farthest3PrimeEnd;
  $self->{_closest3PrimeEnd}      = $closest3PrimeEnd;




  return $self;
}



=item getGeneCounter

  Subroutine to return the internal
  gene counter.

=cut
sub getGeneCounter{
  my($self)=@_;
  return $self->{_geneCounter};
}



=item getArrayOfTranscripts

  Subroutine to return the array of Transcripts.pm
  associated with this object.

=cut
sub getArrayOfTranscripts{
  my($self)=@_;
  return $self->{_transcriptsArray};
}


=item getStrand

  Subroutine to return the strand for the given transcript

=cut
sub getStrand{
  my($self)=@_;
  return $self->{_strand};
}



=item getChromosome

  Subroutine to return the chromosome for the given transcript

=cut
sub getChromosome{
  my($self)=@_;
  return $self->{_chromosome};
}



=item getId

  Subroutine to return the id for the given transcript

=cut
sub getId{
  my($self)=@_;
  return $self->{_id};
}


=item getNumberOfTranscript

  Subroutine to return the number of
  transcript that was found

=cut
sub getNumberOfTranscript{
  my($self)=@_;
  return $self->{_numberOfTranscript};
}



=item hasAltSplicedTranscript

  Subroutine to determine if the
  gene has more than one transcript.
  Return 1 if it does, 0 if it has a
  single transcript

=cut
sub hasAltSplicedTranscript{
  my($self)=@_;
  if($self->{_numberOfTranscript} > 1){
    return 1;
  }else{
    return 0;
  }
}




=item getLongestTranscript

  Subroutine to return the longest transcript for the given transcript

=cut
sub getLongestTranscript{
  my($self)=@_;
  return $self->{_longestTranscript};
}



=item getShortestTranscript

  Subroutine to return the shortest transcript for the given transcript

=cut
sub getShortestTranscript{
  my($self)=@_;
  return $self->{_shortestTranscript};
}





=item getFarthest5PrimeEnd

  Subroutine to return the farthest 5'
  end (the one closer to the chromosomal
   5' end) for all transcripts

=cut
sub getFarthest5PrimeEnd{
  my($self)=@_;
  return $self->{_farthest5PrimeEnd};
}



=item getFarthest3PrimeEnd

  Subroutine to return the farthest 3'
  end (the one closer to the chromosomal
   3' end) for all transcripts

=cut
sub getFarthest3PrimeEnd{
  my($self)=@_;
  return $self->{_farthest3PrimeEnd};
}






=item getClosest5PrimeEnd

  Subroutine to return the closest 5'
  end (the one closer to the 3' end
  of the transcript) for all transcripts

=cut
sub getClosest5PrimeEnd{
  my($self)=@_;
  return $self->{_closest5PrimeEnd};
}



=item getClosest3PrimeEnd

  Subroutine to return the closest 3'
  end (the one closer to the 5' end
  of the transcript) for all transcripts

=cut
sub getClosest3PrimeEnd{
  my($self)=@_;
  return $self->{_closest3PrimeEnd};
}





=item print

  Subroutine to return the
  a string representation of the
  Gene.pm object

=cut
sub print{
  my($self)=@_;
  my $stringToReturn="";
  $stringToReturn.="Gene ".$self->{_id}."\n";
  $stringToReturn.=$self->{_chromosome}.":".$self->{_strand}.":".$self->{_farthest5PrimeEnd}."-".$self->{_farthest3PrimeEnd}."\n";
  $stringToReturn.=$self->{_numberOfTranscript}." variants :\n";

  for(my $transcriptIndex=0;$transcriptIndex<=$#{$self->{_transcriptsArray}};$transcriptIndex++){
    $stringToReturn.=$self->{_transcriptsArray}->[$transcriptIndex]->print();
  }
  return $stringToReturn;
}


=item clusterOverlappingRanges

  Subroutine an array of Ranges.pm
  and returns an array of Ranges that do not overlap.
  Arguments:
   arrayOfRangesTemp : The array of ranges to cluster


=cut
sub clusterOverlappingRanges{
  my($self,$arrayOfRangesTemp)=@_;
  my $arrayOfRanges=[];
  my $i=0;

  while($i<=$#{$arrayOfRangesTemp}){
    my $myCurrentRange=$arrayOfRangesTemp->[$i];

    my $j=($i+1);
    while($j<=$#{$arrayOfRangesTemp}){
      if($myCurrentRange->overlaps($arrayOfRangesTemp->[$j])){
	$myCurrentRange->setEnd( $arrayOfRangesTemp->[$j]->getEnd() )
      }else{
	last;
      }
      $j++;
    }
    $i=$j;
    push(@{$arrayOfRanges},$myCurrentRange);
  }
  return $arrayOfRanges;
}



=item returnAtlSpliced

  Subroutine to return the
  an array of alternatively spliced regions
  and an array of regions that are expressed in every
  Transcript.pm object associated with this Gene.pm object

=cut
sub returnAtlSpliced{
  my($self)=@_;
  #building an array of exons for each transcripts
  my $arrayOfExonsTree=[];
  my %hashOfCoordinates;
  $hashOfCoordinates{$self->getFarthest5PrimeEnd()}=0;
  if(!(exists $hashOfCoordinates{$self->getFarthest3PrimeEnd()}) ){
  }else{
    die "Farthest 3 prime end is the same as the 5 prime end\n";
  }

  for (my $indexTrans=0;$indexTrans<$self->{_numberOfTranscript};$indexTrans++) {
    my @arrayExons=@{$self->{_transcriptsArray}->[$indexTrans]->getExons()};
    foreach my $exonToAdd (@arrayExons) {
      my $segmentTreeRange=new SegmentTreeRange(start => $exonToAdd->getStart(),
						end   => $exonToAdd->getEnd()   );
      if(!(exists $hashOfCoordinates{$exonToAdd->getStart()} ) ){
	$hashOfCoordinates{$exonToAdd->getStart()} = 0;
      }

      if(!(exists $hashOfCoordinates{$exonToAdd->getEnd()} ) ){
	$hashOfCoordinates{$exonToAdd->getEnd()} = 0;
      }

      my $hash={transID    =>  $self->{_transcriptsArray}->[$indexTrans]->getId(),
		range      =>  $segmentTreeRange};
      push(@{$arrayOfExonsTree},$hash);
    }

  }
  my $segTree=new SegmentTree(arrayOfInputs => $arrayOfExonsTree);

  my $arrayOfAltRanges  =[];
  my $arrayOfExprRanges =[];

  my @arraySortedCoord=sort(keys %hashOfCoordinates);


  for(my $i=0;$i<$#arraySortedCoord;$i++){ #skip last one
    if( ( $arraySortedCoord[$i+1] - $arraySortedCoord[$i] ) >1 ) {
      my $arrayOfResults=$segTree->seekCoord( $arraySortedCoord[$i] + 1 );

      if ($#{$arrayOfResults} >= 0 ) {
	if ($#{$arrayOfResults} == ($self->{_numberOfTranscript}-1) ) {
	  my $rangeToPush=Range->new(start => $arraySortedCoord[$i],end   => $arraySortedCoord[$i+1]);
	  push(@{$arrayOfExprRanges},$rangeToPush);
	}else{
	  my $rangeToPush=Range->new(start => $arraySortedCoord[$i],end   => $arraySortedCoord[$i+1]);
	  push(@{$arrayOfAltRanges},$rangeToPush);
	}
      } else {
	#print "no hits ".$arraySortedCoord[$i]." - ".$arraySortedCoord[$i+1]."\n";
      }
    }
  }


  $arrayOfAltRanges  = $self->clusterOverlappingRanges($arrayOfAltRanges);
  $arrayOfExprRanges = $self->clusterOverlappingRanges($arrayOfExprRanges);


  return ($arrayOfAltRanges,$arrayOfExprRanges);
}







1;
