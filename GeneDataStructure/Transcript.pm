=head1 NAME

   GeneDataStructure::Transcript.pm


=head1 Disclaimer

This software and all the data therein were developed by the Bioinformatics and Computational Biology Lab (Laboratorio de Bioinformatica e Biologia Computational (LBBC)) of the Brazilian National Cancer Institute (Instituto Nacional de Cancer (INCA)) in Rio de Janeiro, RJ, Brazil. The INCA is part of the Health Ministry (Ministerio da Saude) of the Federal Government of Brazil. This software was developed by the staff of the LBBC as part of their official duties. According to federal regulations, this software is therefore public domain.

This software is provided as is and the INCA does make any guarantee regarding the accuracy/quality/reliability of the results. The INCA cannot be held accountable for any hardware or software damage that might be caused through the use or misuse of this software.

The LBBC encourages researchers using our software to acknowledge us. The software can be distributed freely and modified given that any derivative bears an acknowledgment that it was derived from our original software.

=head1 Author

Gabriel Renaud grenaud@inca.gov.br



=cut

package Transcript;

use strict;
use warnings;

use Data::Dumper;
use Range;

use vars qw($transcriptCounter);




my $limitForFlagging5pPartial=10; # Limit to use for partial5PrimeEnd() to return 1
my $limitForFlagging3pPartial=15; # Limit to use for partial3PrimeEnd() to return 1


=item new

   Constructor for the Transcript object. This object holds data for
   transcripts (refSeq transcript, ests), it should be used as a superclass
   and the user should only use the child classes.

   The object can be created using a list of exonStart and
   exonEnd as such:

   $trans=Transcript->new(chromosome   => "chr12",
	    	          strand       => "-",
		          id           => "NM_057088",
		          transcStart  => 51469735,
		          transcEnd    => 51476159,
		          exonsStart   => "51469735,51470888,51471256,51471741,51472289,51472761,51473213,51474161,51475448,",
		          exonsEnd     => "51470409,51470923,51471477,51471867,51472454,51472857,51473274,51474382,51476159,");

   or using a list of exonStart and respective size:

   $trans=Transcript->new(chromosome   => "chr1",
	  	          strand       => "+",
		          id           => "BQ428753",
		          transcStart  => 3280,
		          transcEnd    => 3838,
		          exonsStart   => "3280,3760,3788,3800,3811,3819,3826,3833,",
		          exonsSize    => "480,28,12,11,8,6,7,5,");


   The constructor requires the following parameters :

   MANDATORY PARAMETERS:

       chromosome    # chromosome name (must be the same as the one used by UCSC)
       strand        # strand, must be either '+' or '-'
       id            # Identifier for the transcript
       transcStart   # Genomic coordinate of the start of where the transcript aligns (must be 0-based like UCSC)
       transcEnd     # Genomic coordinate of the end of where the transcript aligns (must be 1-based like UCSC)

       Either:
         exonsStart  # Genomic coordinates of where the exons begin (must be 0-based like UCSC)
         exonsSize   # Length for each exon defined above respectively.
       or
         exonsStart  # Genomic coordinates of where the exons begin (must be 0-based like UCSC)
         exonsEnd    # Genomic coordinates of where the exons end   (must be 0-based like UCSC)

       exonsStart, exonsEnd and exonsSize must be in the format "number1,number2,...numberN,"

   OPTIONAL PARAMETERS:

       querySize     # Size of the initial query sequence that was used to map the transcript
       queryStart    # Start coordinate on the initial query sequence where the alignment begins (starts at 0)
       queryEnd      # End coordinate on the initial query sequence where the alignment ends
   You need to specify all 3 querySize, queryStart and queryEnd if you specify one
       species       # The species code is need if getSequence or getExonicSequence is called

=cut
sub new{
  my ($class,%arg)=(@_);

  if($arg{strand} ne "+" &&
     $arg{strand} ne "-" ){
    die "Transcript.pm The strand must either be '+' or '-'\n";
  }

  if($arg{transcStart} !~ /^\d+$/){
    die "Transcript.pm The transcStart must numerical\n";
  }

  if($arg{transcEnd} !~ /^\d+$/){
    die "Transcript.pm The transcEnd must numerical\n";
  }

  my ($self) =bless{_chromosome       => $arg{chromosome},
		    _strand           => $arg{strand},
		    _id               => $arg{id},
		    _exons            => [],
		    _introns          => [],
		   }, $class;

  $self->{_transcriptCounter}=$transcriptCounter++;

  if( (exists $arg{querySize})  ||
      (exists $arg{queryStart}) ||
      (exists $arg{queryEnd})   ){
    if( (exists $arg{querySize})  &&
	(exists $arg{queryStart}) &&
	(exists $arg{queryEnd})   ){

      if($arg{querySize} !~ /^\d+$/){
	die "Transcript.pm The querySize must numerical\n";
      }

      if($arg{queryStart} !~ /^\d+$/){
	die "Transcript.pm The queryStart must numerical\n";
      }

      if($arg{queryEnd} !~ /^\d+$/){
	die "Transcript.pm The queryEnd must numerical\n";
      }

      if($arg{querySize} < $arg{queryStart}){
	die "Transcript.pm The queryStart cannot be greater than querySize\n";
      }

      if($arg{querySize} < $arg{queryEnd}){
	die "Transcript.pm The queryEnd cannot be greater than querySize\n";
      }

      if($arg{queryStart} < 0){
	die "Transcript.pm The queryStart cannot be lesser than 0\n";
      }

      $self->{_queryParamDefined}  = 1;
      $self->{_querySize}  = $arg{querySize};
      $self->{_queryStart} = $arg{queryStart};
      $self->{_queryEnd}   = $arg{queryEnd};



    }else{
      die "Transcript.pm If either one of querySize, queryStart or queryEnd is specified, all three must be specified\n";
    }
  }else{
    $self->{_queryParamDefined}  = 0;
  }


  #Setting transcript range
  $self->{_transRange}=Range->new(start => ($arg{transcStart}+1),
				  end   => $arg{transcEnd});

  #Setting 5prime end 3prime end
  if( $self->{_strand} eq "+" ){
    $self->{_5primeEnd}=$self->{_transRange}->getStart();
    $self->{_3primeEnd}=$self->{_transRange}->getEnd();
  }else{
    $self->{_5primeEnd}=$self->{_transRange}->getEnd();
    $self->{_3primeEnd}=$self->{_transRange}->getStart();
  }





  #############################
  # Setting exons and introns #
  #############################
  if( (exists $arg{exonsStart}) &&
      (exists $arg{exonsSize}) ){
    my $exonStartString = $arg{exonsStart};
    my $exonSizeString  = $arg{exonsSize};
    if(substr($exonStartString,-1) eq ","){
      $exonStartString=substr($exonStartString,0,-1);
    }
    if(substr($exonSizeString,-1) eq ","){
      $exonSizeString=substr($exonSizeString,0,-1);
    }

    my @exonsSt=split(",",$exonStartString);
    my @exonsSi=split(",",$exonSizeString);

    if($#exonsSt != $#exonsSi){
      die "Transcript.pm The number of exons start must be the same as the number of exons size\n";
    }

    if( $self->{_strand} eq "+" ){
    }else{
      @exonsSt=reverse(@exonsSt);
      @exonsSi=reverse(@exonsSi);
    }


    for(my $indexEx=0;$indexEx<=$#exonsSt;$indexEx++){
      if($exonsSt[$indexEx] !~ /^\d+$/){
	die "Transcript.pm The exon start ".$exonsSt[$indexEx]." must be numerical\n";
      }
      if($exonsSi[$indexEx] !~ /^\d+$/){
	die "Transcript.pm The exon size ".$exonsSi[$indexEx]." must be numerical\n";
      }


      my $exon = Range->new(start => ($exonsSt[$indexEx]+1),
			    end   => ($exonsSt[$indexEx]+$exonsSi[$indexEx]));
      push(@{$self->{_exons}},$exon);

      if($indexEx == 0){
	$self->{_firstexon}=$exon;
      }else{

	my $intron;
	my $startIntron;
	my $endIntron;

	if( $self->{_strand} eq "+" ){
	  if( ($exonsSt[$indexEx-1]+$exonsSi[$indexEx-1]) == $exonsSt[$indexEx]){ #an intron of length 1
	    $startIntron  = ($exonsSt[$indexEx-1]+$exonsSi[$indexEx-1]+1);
	    $endIntron    = ($exonsSt[$indexEx-1]+$exonsSi[$indexEx-1]+1);
	  }else{
	    $startIntron  = ($exonsSt[$indexEx-1]+$exonsSi[$indexEx-1]+1);
	    $endIntron    = $exonsSt[$indexEx];
	  }
	}else{
	  if( ($exonsSt[$indexEx-1]) == ($exonsSt[$indexEx]+$exonsSi[$indexEx])){ #an intron of length 1
	    $startIntron  = ($exonsSt[$indexEx]+$exonsSi[$indexEx]+1);
	    $endIntron    = ($exonsSt[$indexEx]+$exonsSi[$indexEx]+1);
	  }else{
	    $startIntron  = ($exonsSt[$indexEx]+$exonsSi[$indexEx]+1);
	    $endIntron    = $exonsSt[$indexEx-1];
	  }
	}

	$intron = Range->new(start => $startIntron,
			     end   => $endIntron);

	push(@{$self->{_introns}},$intron);
	if($indexEx == 1 ){
	  $self->{_firstintron}=$intron;
	}
	if($indexEx == $#exonsSt){
	  $self->{_lastintron}=$intron;
	}
      }

      if($indexEx == $#exonsSt){
	$self->{_lastexon}=$exon;
      }
    }

  }elsif( (exists $arg{exonsStart}) &&
	  (exists $arg{exonsEnd}) ){
    my $exonStartString = $arg{exonsStart};
    my $exonEndString  = $arg{exonsEnd};
    if(substr($exonStartString,-1) eq ","){
      $exonStartString=substr($exonStartString,0,-1);
    }
    if(substr($exonEndString,-1) eq ","){
      $exonEndString=substr($exonEndString,0,-1);
    }

    my @exonsSt=split(",",$exonStartString);
    my @exonsEn=split(",",$exonEndString);


    if($#exonsSt != $#exonsEn){
      die "Transcript.pm The number of exons start must be the same as the number of exons end\n";
    }

    if( $self->{_strand} eq "+" ){
    }else{
      @exonsSt=reverse(@exonsSt);
      @exonsEn=reverse(@exonsEn);
    }



    for(my $indexEx=0;$indexEx<=$#exonsSt;$indexEx++){
      if($exonsSt[$indexEx] !~ /^\d+$/){
	die "Transcript.pm The exon start ".$exonsSt[$indexEx]." must be numerical\n";
      }
      if($exonsEn[$indexEx] !~ /^\d+$/){
	die "Transcript.pm The exon size ".$exonsEn[$indexEx]." must be numerical\n";
      }


      my $exon = Range->new(start => ($exonsSt[$indexEx]+1),
			    end   => ($exonsEn[$indexEx]));

      push(@{$self->{_exons}},$exon);

      if($indexEx == 0){
	$self->{_firstexon}=$exon;
      }else{
	my $intron;
	my $startIntron;
	my $endIntron;

	if( $self->{_strand} eq "+" ){
	  if($exonsEn[$indexEx-1] == $exonsSt[$indexEx]){
	    $startIntron =($exonsEn[$indexEx-1]+1);
	    $endIntron   =($exonsEn[$indexEx-1]+1);
	  }else{
	    $startIntron = ($exonsEn[$indexEx-1]+1);
	    $endIntron   =  $exonsSt[$indexEx];
	  }
	}else{
	  if($exonsEn[$indexEx] == $exonsSt[$indexEx-1]){
	    $startIntron = $exonsEn[$indexEx]+1;
	    $endIntron   = $exonsEn[$indexEx]+1;
	  }else{
	    $startIntron = $exonsEn[$indexEx]+1;
	    $endIntron   = $exonsSt[$indexEx-1];
	  }
	}


	$intron = Range->new(start => $startIntron,
			     end   => $endIntron);

	push(@{$self->{_introns}},$intron);
	if($indexEx == 1 ){
	  $self->{_firstintron}=$intron;
	}
	if($indexEx == $#exonsSt){
	  $self->{_lastintron}=$intron;
	}
      }

      if($indexEx == $#exonsSt){
	$self->{_lastexon}=$exon;
      }
    }


  }else{
    die "Transcript.pm Please specify either exonsStart and exonsSize or, exonsStart and exonsEnd\n";
  }


  if( $self->{_strand} eq "+" ){
    if($self->{_firstexon}->getStart() != $self->{_5primeEnd}){
      die "Transcript.pm The start of the first exon ".$self->{_firstexon}->getEnd()." does not match the 5' end ".$self->{_5primeEnd}."\n";
    }
    if($self->{_lastexon}->getEnd() != $self->{_3primeEnd}){
      die "Transcript.pm The end of the last exon ".$self->{_lastexon}->getEnd()." does not match the 3' end ".$self->{_3primeEnd}."\n";
    }
  }else{
    if($self->{_firstexon}->getEnd() != $self->{_5primeEnd}){
      die "Transcript.pm The end of the first exon ".$self->{_firstexon}->getEnd()." does not match the 5' end ".$self->{_5primeEnd}."\n";
    }
    if($self->{_lastexon}->getStart() != $self->{_3primeEnd}){
      die "Transcript.pm The start of the last exon ".$self->{_lastexon}->getEnd()." does not match the 3' end ".$self->{_3primeEnd}."\n";
    }
  }

  return $self;
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


=item getTranscriptionRange

  Subroutine to return the transcription range
  (as a Range object) for the given transcript

=cut
sub getTranscriptionRange{
  my($self)=@_;
  return $self->{_transRange};
}




=item getTranscriptCounter

  Subroutine to return the internal
  transcript counter.

=cut
sub getTranscriptCounter{
  my($self)=@_;
  return $self->{_transcriptCounter};
}


=item getGeneCounter

  Subroutine to return the internal
  counter for the parent gene if it is set.
  Returns transcriptCounter otherwise

=cut
sub getGeneCounter{
  my($self)=@_;

  if($self->{_parentGene}){
    return $self->getGene()->getGeneCounter();
  }else{
    return $self->{_transcriptCounter};
  }
}





=item getFirstExon

  Subroutine to return the first exon
  (as a Range object) for the given transcript

=cut
sub getFirstExon{
  my($self)=@_;
  return $self->{_firstexon};
}

=item getLastExon

  Subroutine to return the last exon
  (as a Range object) for the given transcript

=cut
sub getLastExon{
  my($self)=@_;
  return $self->{_lastexon};
}


=item getFirstIntron

  Subroutine to return the first intron
  (as a Range object) for the given transcript

=cut
sub getFirstIntron{
  my($self)=@_;
  return $self->{_firstintron};
}

=item getLastIntron

  Subroutine to return the last intron
  (as a Range object) for the given transcript

=cut
sub getLastIntron{
  my($self)=@_;
  return $self->{_lastintron};
}


=item get5Prime

  Subroutine to return the 5' end
  for the given transcript

=cut
sub get5Prime{
  my($self)=@_;
  return $self->{_5primeEnd};
}

=item get3Prime

  Subroutine to return the 3' end
  for the given transcript

=cut
sub get3Prime{
  my($self)=@_;
  return $self->{_3primeEnd};
}


=item getExons

  Subroutine to return the exons as
  an array of Range objects

=cut
sub getExons{
  my($self)=@_;
  return $self->{_exons};
}

=item getIntrons

  Subroutine to return the introns as
  an array of Range objects

=cut
sub getIntrons{
  my($self)=@_;
  return $self->{_introns};
}





=item getExonCount

  Subroutine to return the
  number of exons

=cut
sub getExonCount{
  my($self)=@_;
  return ($#{$self->{_exons}}+1);
}

=item getIntronsCount

  Subroutine to return the
  number of introns

=cut
sub getIntronCount{
  my($self)=@_;
  return ($#{$self->{_introns}}+1);
}





=item getQuerySize

  Subroutine to return the size of the
  initial query sequence that was
  used to map the transcript

=cut
sub getQuerySize{
  my($self)=@_;

  if(!($self->{_queryParamDefined})){
    die "Transcript.pm Cannot call getQuerySize() if the object was contructed without defining querySize";
  }
  return $self->{_querySize};
}


=item getQueryStart

  Subroutine to return the start
  coordinate of the initial query
  sequence that was used to map
  the transcript

=cut
sub getQueryStart{
  my($self)=@_;

  if(!($self->{_queryParamDefined})){
    die "Cannot call getQueryStart() if the object was contructed without defining queryStart";
  }
  return $self->{_queryStart};
}



=item getQueryEnd

  Subroutine to return the end
  coordinate of the initial query
  sequence that was used to map
  the transcript

=cut
sub getQueryEnd{
  my($self)=@_;

  if(!($self->{_queryParamDefined})){
    die "Cannot call getQueryEnd() if the object was contructed without defining queryEnd";
  }
  return $self->{_queryEnd};
}


=item partial5PrimeEnd

  Subroutine to determine if the 5' end of the
  initial query sequence that was used to map
  the transcript aligns partially
  to the genome. Returns 1 if the getQueryStart()
  is greater than limitForFlagging5pPartial

=cut
sub partial5PrimeEnd{
  my($self)=@_;

  if(!($self->{_queryParamDefined})){
    die "Transcript.pm Cannot call getQueryEnd() if the object was contructed without defining queryEnd";
  }

  if($self->{_queryStart} >= $limitForFlagging5pPartial){
    return 1;
  }else{
    return 0;
  }
}







=item partial3PrimeEnd

  Subroutine to determine if the 3' end of the
  initial query sequence that was used to map
  the transcript aligns partially
  to the genome. Returns 1 if the difference between
  getQuerySize() and getQueryEnd() is greater than
  limitForFlagging3pPartial

=cut
sub partial3PrimeEnd{
  my($self)=@_;

  if(!($self->{_queryParamDefined})){
    die "Transcript.pm Cannot call getQueryEnd() if the object was contructed without defining queryEnd";
  }

  if(  ($self->{_querySize}-$self->{_queryEnd}) >= $limitForFlagging3pPartial){
    return 1;
  }else{
    return 0;
  }
}



=item getGene

  Subroutine to return the
  Gene.pm object that contains
  the transcript if it was set.

=cut
sub getGene{
  my($self)=@_;

  if( !($self->{_parentGene}) ){
    print Dumper($self);
    die "Transcript.pm Cannot call getGene() if the Gene.pm object was not defined";
  }
  return $self->{_parentGene};
}



=item setGene

  Subroutine to set the
  Gene.pm object that contains
  the transcript.

=cut
sub setGene{
  my($self,$gene)=@_;
  if($self->{_parentGene}){
    die "Transcript.pm The gene was already set\n";
  }

  $self->{_parentGene}=$gene;
}




=item getGeneSymbol

  Subroutine to return the gene symbol
  associated with this transcript.

=cut
sub getGeneSymbol{
  my($self)=@_;
  if( $self->{_internalGeneSymbol} ){
    return $self->{_internalGeneSymbol};
  }else{
    return $self->{_parentGene}->getId();
  }
}



=item setGeneSymbol

  Subroutine to set the force the transcript to
  use that gene symbol instead of querying the parent gene
  of the transcript.

=cut
sub setGeneSymbol{
  my($self,$geneSymbol)=@_;
  if($self->{_internalGeneSymbol}){
    die "Transcript.pm The internalGeneSymbol was already set\n";
  }
  $self->{_internalGeneSymbol}=$geneSymbol;
}





=item overlaps

  Subroutine to determine if this transcripts
  overlaps another transcript. Does not take strand into account.

=cut
sub overlaps{
  my($self,$otherTranscript)=@_;

  if($self->getChromosome() ne $otherTranscript->getChromosome()){
    return 0;
  }

  return ($self->getTranscriptionRange()->overlaps( $otherTranscript->getTranscriptionRange() ));

}



=item overlapsWithStrand

  Subroutine to determine if this transcripts
  overlaps another transcript. Takes strand into account.

=cut
sub overlapsWithStrand{
  my($self,$otherTranscript)=@_;

  if($self->getChromosome() ne $otherTranscript->getChromosome()){
    return 0;
  }
  if($self->getStrand() ne $otherTranscript->getStrand()){
    return 0;
  }

  return ($self->getTranscriptionRange()->overlaps( $otherTranscript->getTranscriptionRange() ));

}



=item overlapsCoordinate

  Subroutine to determine if this transcripts
  overlaps a given coordinate.

=cut
sub overlapsCoordinate{
  my($self,$coordinate)=@_;
  return ($self->getTranscriptionRange()->contains( $coordinate ));

}




=item overlapsRange

  Subroutine to determine if this transcripts
  overlaps a given coordinate.

=cut
sub overlapsRange{
  my($self,$range)=@_;
  return ( $self->getTranscriptionRange()->overlaps($range) );
}







=item annotate

  Subroutine to determine the relative position
  of a coordinate with respect to this transcript
  given the chromosome and the range

=cut
sub annotate{
  my($self,$chr,$coordinate,$radius)=@_;

  if($chr ne $self->getChromosome() ){ #Not on the same chromosome
    return 0;
  }else{# On the same chromosome

    my $range = Range->new(start => $coordinate-$radius,
			   end   => $coordinate+$radius);
    if( ! ($self->overlapsRange($range)) ){ #if the range does not overlap the transcription range
      return 0;
    }else{ #if the range overlaps the transcription range
      return annotateCoordinate($self,$coordinate);
    }#end overlap check

  }#end chromosome check

  die "Transcription.pm Invalid state ".Dumper($self)." coordinate $coordinate\n";

}




=item annotateCoordinate

  Subroutine to determine the relative position
  of a coordinate with respect to this transcript
  given that we know that the chromosome is the same
  and the range is acceptable.
  Code used by the returned hash:
  key        value meaning
  'code'     1     exon
             2     intron
             3     upstream
             4     downstream

  'partial'  #     flag to determine if the 5' or 3'
             #     end is partial for up/downstream cases

  'index'    #     index of the exon or intron
                   (only used when code = 3 or 4)

  'distance5P'    #     distance to the 5 prime end
  'distance3P'    #     distance to the 3 prime end


=cut
sub annotateCoordinate{
  my($self,$coordinate)=@_;
  #print "coordinate $coordinate\n";
  if($coordinate !~ /\d+/){
    die "Error, Transcript.pm The coordinate must be numerical\n";
  }

  my $resultsToReturn={};

  if( $self->getTranscriptionRange()->contains($coordinate) ){ #if the coordinate is within the transcript
    $resultsToReturn->{distance5P}=abs($self->get5Prime() - $coordinate);
    $resultsToReturn->{distance3P}=abs($self->get3Prime() - $coordinate);

    for(my $indexExon=0;$indexExon<$self->getExonCount();$indexExon++){
      if($self->getExons()->[$indexExon]->contains($coordinate)){
	$resultsToReturn->{code}=1;
	$resultsToReturn->{index}=($indexExon+1);
	return $resultsToReturn;
      }
    }



    for(my $indexIntron=0;$indexIntron<$self->getIntronCount();$indexIntron++){
      if($self->getIntrons()->[$indexIntron]->contains($coordinate)){
	$resultsToReturn->{code}=2;
	$resultsToReturn->{index}=($indexIntron+1);
	return $resultsToReturn;
      }
    }

    die "Transcription.pm annotateCoordinate1 Invalid state ".Dumper($self)." coordinate $coordinate\n";
  }else{
    if($self->getStrand eq "+"){
      if($coordinate < $self->get5Prime()){     #upstream
	$resultsToReturn->{code}=3;
	$resultsToReturn->{distance5P}=$self->get5Prime() - $coordinate;
	$resultsToReturn->{distance3P}=$self->get3Prime() - $coordinate;


	if($self->{_queryParamDefined} && $self->partial5PrimeEnd()){
	  $resultsToReturn->{partial}=1;
	}else{
	  $resultsToReturn->{partial}=0;
	}

	return $resultsToReturn;
      }elsif($coordinate > $self->get3Prime()){ #downstream
	$resultsToReturn->{code}=4;
	$resultsToReturn->{distance5P}=$coordinate - $self->get5Prime();
	$resultsToReturn->{distance3P}=$coordinate - $self->get3Prime();

	if($self->{_queryParamDefined} && $self->partial3PrimeEnd()){
	  $resultsToReturn->{partial}=1;
	}else{
	  $resultsToReturn->{partial}=0;
	}

	return $resultsToReturn;
      }else{

	die "Transcription.pm Invalid state annotateCoordinate2 ".Dumper($self)." coordinate $coordinate\n";
      }

    }elsif($self->getStrand eq "-"){
      if($coordinate > $self->get5Prime()){     #upstream
	$resultsToReturn->{code}=3;
	$resultsToReturn->{distance5P}=$coordinate-$self->get5Prime();
	$resultsToReturn->{distance3P}=$coordinate-$self->get3Prime();

	if($self->{_queryParamDefined} && $self->partial5PrimeEnd()){
	  $resultsToReturn->{partial}=1;
	}else{
	  $resultsToReturn->{partial}=0;
	}

	return $resultsToReturn;
      }elsif($coordinate < $self->get3Prime()){ #downstream
	$resultsToReturn->{code}=4;
	$resultsToReturn->{distance5P}=$self->get5Prime()-$coordinate;
	$resultsToReturn->{distance3P}=$self->get3Prime()-$coordinate;

	if($self->{_queryParamDefined} && $self->partial3PrimeEnd()){
	  $resultsToReturn->{partial}=1;
	}else{
	  $resultsToReturn->{partial}=0;
	}

	return $resultsToReturn;
      }else{
	die "Transcription.pm Invalid state annotateCoordinate3 ".Dumper($self)." coordinate $coordinate\n";
      }
    }else{
      die "Transcription.pm Invalid strand\n";
    }
  }

  die "Transcription.pm Invalid state annotateCoordinate4 ".Dumper($self)." coordinate $coordinate\n";
}


=item annotateRange

  Subroutine to determine the relative position
  of a range with respect to this transcript
  given that we know that the chromosome is the same

  Code used by the returned hash:
  key        value meaning
  'code'     1     contained in transcript and totally contained in a single exon
             2     contained in transcript and overlaps exonic and intronic regions
             3     contained in transcript and totally contained in intronic region
             4     range completely contains the transcript
             5     overlaps 5' end and first exon only
             6     overlaps 5' end and exons, introns
             7     overlaps 3' end and last exon only
             8     overlaps 3' end and exons, introns
             9     upstream
            10     downstream

  'partial'  #     flag to determine if the 5' or 3'
             #     end is partial for up/downstream cases

  'indicesEx'        #     indices of the exons
                     (only used when code = 1,2,5,6,7,8)
  'indicesIntron'    #     indices of the introns
                     (only used when code = 2,3,6,8)


=cut
sub annotateRange{
  my($self,$coordinate1,$coordinate2)=@_;
  if ($coordinate1 > $coordinate2) {
    die "Transcript.pm annotateRange() The field coordinate1 $coordinate1 cannot be greater than coordinate2 $coordinate2\n";
  }
  my $rangeToAnnotate=new Range(start => $coordinate1,
				end   => $coordinate2);

  my $resultsToReturn={};
  $resultsToReturn->{indicesEx}     ="";
  $resultsToReturn->{indicesIntron} ="";

  if ( $self->getTranscriptionRange()->overlaps($rangeToAnnotate) ) { #if at least one coordinate is within the transcript

    if ( $self->getTranscriptionRange()->contains($coordinate1)  &&
	 $self->getTranscriptionRange()->contains($coordinate2) ) { #if both coordinates is within the transcript

      my @exonsToAdd;
      for (my $indexExon=0;$indexExon<$self->getExonCount();$indexExon++) {
	if ($self->getExons()->[$indexExon]->overlaps($rangeToAnnotate)) {
	  push(@exonsToAdd,$indexExon+1);
	}
      }

      my @intronsToAdd;
      for (my $indexIntron=0;$indexIntron<$self->getIntronCount();$indexIntron++) {
	if ($self->getIntrons()->[$indexIntron]->overlaps($rangeToAnnotate)) {
	  push(@intronsToAdd,$indexIntron+1);
	}
      }
      $resultsToReturn->{indicesEx}     = join(",",@exonsToAdd);
      $resultsToReturn->{indicesIntron} = join(",",@intronsToAdd);

      if (      ($#exonsToAdd == 0)  && ($#intronsToAdd == -1)  ) {
	$resultsToReturn->{code}=1; # 1     contained in transcript and totally contained in a single exon
      } elsif ( ($#exonsToAdd == -1) && ($#intronsToAdd == 0)  ) {
	$resultsToReturn->{code}=3; # 3     contained in transcript and totally contained in intronic region
      } else {
	$resultsToReturn->{code}=2; # 2     contained in transcript and overlaps exonic and intronic regions
      }

      return  $resultsToReturn;

    } elsif ( ($self->getStrand eq "+" && ( ($coordinate1 < $self->get5Prime()) &&  ($self->getTranscriptionRange()->contains($coordinate2)) ) ) ||
	      ($self->getStrand eq "-" && ( ($coordinate2 > $self->get5Prime()) &&  ($self->getTranscriptionRange()->contains($coordinate1)) ) ) ) {


      if ($self->{_queryParamDefined} && $self->partial5PrimeEnd()) {
	$resultsToReturn->{partial}=1;
      } else {
	$resultsToReturn->{partial}=0;
      }

      my @exonsToAdd;
      for (my $indexExon=0;$indexExon<$self->getExonCount();$indexExon++) {
	if ($self->getExons()->[$indexExon]->overlaps($rangeToAnnotate)) {
	  push(@exonsToAdd,$indexExon+1);
	}
      }

      my @intronsToAdd;
      for (my $indexIntron=0;$indexIntron<$self->getIntronCount();$indexIntron++) {
	if ($self->getIntrons()->[$indexIntron]->overlaps($rangeToAnnotate)) {
	  push(@intronsToAdd,$indexIntron+1);
	}
      }
      $resultsToReturn->{indicesEx}     = join(",",@exonsToAdd);
      $resultsToReturn->{indicesIntron} = join(",",@intronsToAdd);


      if ( ($#exonsToAdd == 0) && ($#intronsToAdd == -1) && ($exonsToAdd[0] == 1) ) {
	$resultsToReturn->{code}=5; # 5     overlaps 5' end and first exon only
      } else {
	$resultsToReturn->{code}=6; # 6     overlaps 5' end and exons, introns
      }

      return $resultsToReturn;

    } elsif ( ($self->getStrand eq "+" && ( ($self->getTranscriptionRange()->contains($coordinate1)) && ($coordinate2 > $self->get3Prime()) ) ) ||
	      ($self->getStrand eq "-" && ( ($self->getTranscriptionRange()->contains($coordinate2)) && ($coordinate1 < $self->get3Prime()) ) ) ) { #downstream

      if ($self->{_queryParamDefined} && $self->partial3PrimeEnd()) {
	$resultsToReturn->{partial}=1;
      } else {
	$resultsToReturn->{partial}=0;
      }


      my @exonsToAdd;
      for (my $indexExon=0;$indexExon<$self->getExonCount();$indexExon++) {
	if ($self->getExons()->[$indexExon]->overlaps($rangeToAnnotate)) {
	  push(@exonsToAdd,$indexExon+1);
	}
      }

      my @intronsToAdd;
      for (my $indexIntron=0;$indexIntron<$self->getIntronCount();$indexIntron++) {
	if ($self->getIntrons()->[$indexIntron]->overlaps($rangeToAnnotate)) {
	  push(@intronsToAdd,$indexIntron+1);
	}
      }
      $resultsToReturn->{indicesEx}     = join(",",@exonsToAdd);
      $resultsToReturn->{indicesIntron} = join(",",@intronsToAdd);


      if ( ($#exonsToAdd == 0) && ($#intronsToAdd == -1) && ($exonsToAdd[0] == $self->getExonCount()) ) {
	$resultsToReturn->{code}=7; # 7     overlaps 3' end and last exon only
      } else {
	$resultsToReturn->{code}=8; # 8     overlaps 3' end and exons, introns
      }

      return $resultsToReturn;

    } elsif ( $rangeToAnnotate->contains($self->get5Prime()) &&
	      $rangeToAnnotate->contains($self->get3Prime())  ) { #4     range completely contains the transcript
      $resultsToReturn->{code}=4;
      return $resultsToReturn;
    } else {
      die "Transcription.pm Invalid state annotateRange1 ".Dumper($self)." coordinate1 $coordinate1 coordinate2 $coordinate2\n";
    }

  } else {		 #both coordinates are outside the transcript.

    if ($self->getStrand eq "+") {
      if ( ($coordinate1 < $self->get5Prime()) &&
	   ($coordinate2 < $self->get5Prime()) ) { #9     upstream
	$resultsToReturn->{code}=9;

	if ($self->{_queryParamDefined} && $self->partial5PrimeEnd()) {
	  $resultsToReturn->{partial}=1;
	} else {
	  $resultsToReturn->{partial}=0;
	}

	return $resultsToReturn;
      } elsif (  ($coordinate1 > $self->get3Prime()) &&
		 ($coordinate2 > $self->get3Prime()) ) { #10     downstream
	$resultsToReturn->{code}=10;

	if ($self->{_queryParamDefined} && $self->partial3PrimeEnd()) {
	  $resultsToReturn->{partial}=1;
	} else {
	  $resultsToReturn->{partial}=0;
	}

	return $resultsToReturn;
      } else {
	  #die "Transcription.pm Invalid state\n";
	  die "Transcription.pm Invalid state annotateRange2 ".Dumper($self)." coordinate1 $coordinate1 coordinate2 $coordinate2\n";
      }
    } elsif ($self->getStrand eq "-") {
      if (  ($coordinate1 > $self->get5Prime()) &&
	    ($coordinate2 > $self->get5Prime()) ) { #9     upstream
	$resultsToReturn->{code}=9;

	if ($self->{_queryParamDefined} && $self->partial5PrimeEnd()) {
	  $resultsToReturn->{partial}=1;
	} else {
	  $resultsToReturn->{partial}=0;
	}

	return $resultsToReturn;
      } elsif ( ($coordinate1 < $self->get3Prime()) &&
		($coordinate2 < $self->get3Prime()) ) { #10     downstream
	$resultsToReturn->{code}=10;

	if ($self->{_queryParamDefined} && $self->partial3PrimeEnd()) {
	  $resultsToReturn->{partial}=1;
	} else {
	  $resultsToReturn->{partial}=0;
	}

	return $resultsToReturn;
      } else {
	die "Transcription.pm Invalid state annotateRange3 ".Dumper($self)." coordinate1 $coordinate1 coordinate2 $coordinate2\n";
      }
    } else {
      die "Transcription.pm Invalid strand\n";
    }
  }

}


=item print

  Subroutine to return the
  a string representation of the
  Transcript.pm object

=cut
sub print{
  my($self)=@_;
  return $self->{_id}."#".$self->{_chromosome}.":".$self->{_strand}.":".$self->{_transRange}->print()."\n";
}




1;
