=head1 NAME

   GeneDataStructure::CodingTranscript.pm


=head1 Disclaimer

This software and all the data therein were developed by the Bioinformatics and Computational Biology Lab (Laboratorio de Bioinformatica e Biologia Computational (LBBC)) of the Brazilian National Cancer Institute (Instituto Nacional de Cancer (INCA)) in Rio de Janeiro, RJ, Brazil. The INCA is part of the Health Ministry (Ministerio da Saude) of the Federal Government of Brazil. This software was developed by the staff of the LBBC as part of their official duties. According to federal regulations, this software is therefore public domain.

This software is provided as is and the INCA does make any guarantee regarding the accuracy/quality/reliability of the results. The INCA cannot be held accountable for any hardware or software damage that might be caused through the use or misuse of this software.

The LBBC encourages researchers using our software to acknowledge us. The software can be distributed freely and modified given that any derivative bears an acknowledgment that it was derived from our original software.

=head1 Author

Gabriel Renaud grenaud@inca.gov.br



=cut


package CodingTranscript;

use strict;
use warnings;
use Cwd;
use Cwd 'abs_path';
use Data::Dumper;
use File::Basename;
use Time::HiRes;


use Range;
use CodingExon;
use Transcript;
our @ISA = "Transcript";


my %complement = ('A' => 'T',
		  'C' => 'G',
		  'G' => 'C',
		  'T' => 'A');



=item new

   Constructor for the CodingTranscript object. This object holds data for
   a coding transcript. It inherits from Transcript.pm. It merely adds the
   coding range parameters.

   A CodingTranscript.pm object can be constructed as such:

   $trans=CodingTranscript->new(chromosome   => "chr19",
                                strand       => "+",
                                id           => "NM_020415",
                                transcStart  => 7639971,
                                transcEnd    => 7641340,
                                codingStart  => 7640212,
                                codingEnd    => 7641235,
                                exonsStart   => "7639971,7640202,7640706,7641104,",
                                exonsEnd     => "7640007,7640330,7640784,7641340,",
   Optional parameter:
                                baseDir      => "/directory/ucsc/data/"
                                seqExtractor => "/directory/software/seq-substring/seq-substring"
                                framesFromFile=> "-1,0,1,1,"
);

=cut

sub new{
  my ($class,%arg)=(@_);

  my $self = $class->SUPER::new(%arg);

  if($arg{codingStart} !~ /^\d+$/){
    die "CodingTranscript.pm The codingStart must numerical\n";
  }



  if($arg{baseDir}){
    $self->{_baseDir}=$arg{baseDir};
  }


  if($arg{seqExtractor}){
    $self->{_seqExtractor}=$arg{seqExtractor};
  }

  if($arg{codingEnd} !~ /^\d+$/){
    die "CodingTranscript.pm The codingEnd must numerical\n";
  }

  if ($arg{codingStart} == $arg{codingEnd}) {
    $self->{_isCoding}=0;
  } else {
    $self->{_isCoding}=1;

    $self->{_codingRange}=Range->new(start => ($arg{codingStart}+1),
				     end   =>  $arg{codingEnd});

    if ($self->{_transRange}->getStart() > $self->{_codingRange}->getStart()) {
      die "CodingTranscript.pm The transcription start ".$self->{_transRange}->getStart()." cannot be greater than the coding start ".$arg{codingStart}."\n";
    }


    if ($self->{_transRange}->getEnd()   <  $self->{_codingRange}->getStart()) {
      die "CodingTranscript.pm The transcription end ".$self->{_transRange}->getEnd()." cannot be lesser than the coding start ".$arg{codingStart}."\n";
    }


    if ($self->{_transRange}->getEnd()   <  $self->{_codingRange}->getEnd()) {
      die "CodingTranscript.pm The transcription end ".$self->{_transRange}->getEnd()." cannot be greater than the coding end ".$arg{codingEnd}."\n";
    }

    if ($self->{_transRange}->getStart() > $self->{_codingRange}->getEnd()) {
      die "CodingTranscript.pm The transcription start ".$self->{_transRange}->getStart()." cannot be greater than the coding end ".$arg{codingEnd}."\n";
    }



    if ( $self->{_strand} eq "+" ) {
      $self->{_5primeEndCoding}=$self->{_codingRange}->getStart();
      $self->{_3primeEndCoding}=$self->{_codingRange}->getEnd();
      if($arg{framesFromFile}){
	$self->{_framesFromFile}=[split(",",$arg{framesFromFile})];
      }
    } else {
      $self->{_5primeEndCoding}=$self->{_codingRange}->getEnd();
      $self->{_3primeEndCoding}=$self->{_codingRange}->getStart();
      if($arg{framesFromFile}){
	$self->{_framesFromFile}=[reverse(split(",",$arg{framesFromFile}))];
      }

    }
  }

  $self->{_isReliable}=1;


  return $self;
}



=item isReliable

  Subroutine to determine if the transcript
  contains oddities such as :
   * frames that do not match the ones on file
   * does not begin with a START codon
   * does not end with a STOP codon
   * The length of the coding seq
     is not a multiple of 3

=cut
sub isReliable{
  my($self)=@_;
  return $self->{_isReliable};
}



=item isCoding

  Subroutine to determine if the transcript
  is coding or not. Return 1 if coding, 0 otherwise.
  This subroutine will return 1 if the coding start
  is equal to the coding end.

=cut
sub isCoding{
  my($self)=@_;
  return $self->{_isCoding};
}


=item getCodingRange

  Subroutine to return the coding range
  (as a Range object) for the given transcript

=cut
sub getCodingRange{
  my($self)=@_;
  if(!$self->{_isCoding}){
    die "getCodingRange() The current transcript ".Dumper($self)."is not coding\n";
  }
  return $self->{_codingRange};
}




=item get5PrimeCoding

  Subroutine to return the 5' end
  of the coding part of the given transcript

=cut
sub get5PrimeCoding{
  my($self)=@_;
  if(!$self->{_isCoding}){
    die "get5PrimeCoding() The current transcript ".Dumper($self)." is not coding\n";
  }

  return $self->{_5primeEndCoding};
}

=item get3PrimeCoding

  Subroutine to return the 3' end
  of the coding part of the given transcript

=cut
sub get3PrimeCoding{
  my($self)=@_;
  if(!$self->{_isCoding}){
    die "get3PrimeCoding() The current transcript ".Dumper($self)." is not coding\n";
  }

  return $self->{_3primeEndCoding};
}





=item isCoordinateInCodingRange

  Subroutine to determine if a given coordinate is in
  the coding range for a given object.

=cut
sub isCoordinateInCodingRange{
  my($self,$coordinate)=@_;
  if(!$self->{_isCoding}){
    die "isCoordinateInCodingRange() The current transcript ".Dumper($self)."is not coding\n";
  }
  return $self->{_codingRange}->contains($coordinate);
}



=item isRangeInCodingRange

  Subroutine to determine if a given range object is in
  the coding range for a given object.

=cut
sub isRangeInCodingRange{
  my($self,$range)=@_;
  if(!$self->{_isCoding}){
    die "isCoordinateInCodingRange() The current transcript ".Dumper($self)."is not coding\n";
  }
  return $self->{_codingRange}->overlaps($range);
}



=item computeNextFrame

  Subroutine called by buildExonCoding to compute the frame for the next exon
  using the frame # of the current frame and the length of the (current exon mod 3).

=cut
sub computeNextFrame{
  my($currentFrame,$mod3length)=@_;
  my $nextFrame;


  if ($currentFrame == -1 || $currentFrame == 0 ) {

    if ($mod3length     == 0 ) {
      $nextFrame = 0;
    } elsif ($mod3length == 1 ) {
      $nextFrame = 1;
    } elsif ($mod3length == 2 ) {
      $nextFrame = 2;
    } else {
      die "Wrong mod $mod3length\n";
    }

  } elsif ($currentFrame == 1 ) {

    if ($mod3length     == 0 ) {
      $nextFrame = 1;
    } elsif ($mod3length == 1 ) {
      $nextFrame = 2;
    } elsif ($mod3length == 2 ) {
      $nextFrame = 0;
    } else {
      die "Wrong mod $mod3length\n";
    }

  } elsif ($currentFrame == 2 ) {


    if ($mod3length     == 0 ) {
      $nextFrame = 2;
    } elsif ($mod3length == 1 ) {
      $nextFrame = 0;
    } elsif ($mod3length == 2 ) {
      $nextFrame = 1;
    } else {
      die "Wrong mod $mod3length\n";
    }

  } else {
    die "Wrong frame $currentFrame\n";
  }

  return $nextFrame;
}


=item extractSequence

  Subroutine to use the sequence extractor program to extract
  the sequence of the coding exon from the chromosome file
  and return it. The internal variables _seqExtractor and
  _baseDir must be set. The chromosome file must contain only
  one
  Arguments:
   startChr  : The start coordinate on the chromosome
   endChr    : The end coordinate on the chromosome
   chrFile   : The name of the file containing the chromosome

=cut
sub extractSequence{
  my($self,$startChr,$endChr,$chrFile)=@_;

  my $command=$self->{_seqExtractor}." -s ".$startChr." -e ".$endChr." -f ".$self->{_baseDir}.$chrFile.".fa -d";
  my $sequenceToReturn=`$command`;

  if ( $? == -1 ){
    print "Error CodingTranscript.pm command failed: $!\n";
  }

  $sequenceToReturn=uc(join("",split("\n",$sequenceToReturn)));

  return $sequenceToReturn;

}





=item buildExonCoding

  Subroutine to set the genomic sequence for each coding exon
  that are part of the coding transcript. This subroutine must be
  called before calling annotateSNP(), annotateINSERT, annotateDELETION()
  The steps are:
   1) Set the $self->{_seqExtractor} variable
   2) Compute the frames for each exon, store them in $self->{_framesCalculated}
   3) If $self->{_framesFromFile} differ from $self->{_framesCalculated}, set the transcript as unreliable
   4) Extract the entire sequence, taking into account the strand
   5) Set each part of the genomic sequence to each CodingExon.pm object if they are within the protein coding region
   6) Put each CodingExon.pm into $self->{_codingExons}
   7) For each CodingExon.pm, append at the start or end of the coding exon the base pairs needed to complete the
      codons at the beginning or the end.


=cut
sub buildExonCoding{
  my($self)=@_;


  if (!$self->isCoding()) {
    warn "Cannot call buildExonCoding() on non-coding transcript ".$self->getId()."\n";
  } else {

    #Detecting the seq-substring program
    if (!($self->{_baseDir})) {
      die "Error in CodingTranscript The base directory must be set to call the buildExonCoding() subroutine\n";
    }

    if (!($self->{_seqExtractor})) {
      my $abs_path = dirname(abs_path(__FILE__));
      my @arrayOfStrings=split("/",$abs_path);

      delete $arrayOfStrings[$#arrayOfStrings];

      $self->{_seqExtractor}=join("/",@arrayOfStrings)."/seq-substring/seq-substring";

    }






























    #BEGIN DETECTING FRAMES



    my $currentFrame=-1;
    my $inCoding=0;
    $self->{_codingExons}=[];
    my $codingExonIndex=0;

    $self->{_framesCalculated}=[];

    for (my $exonToAddIndex=0 ; $exonToAddIndex<$self->getExonCount() ; $exonToAddIndex++) {
      my $exonToAdd=$self->getExons()->[$exonToAddIndex];


      if ($exonToAdd->contains($self->get5PrimeCoding()) && $exonToAdd->contains($self->get3PrimeCoding())) { #case where the coding st and end is contained in a single exon
	my $length=abs($self->get3PrimeCoding()-$self->get5PrimeCoding())+1;

	$currentFrame=0;
	push(@{$self->{_framesCalculated}},$currentFrame);

	$currentFrame=-1;
      } else {			#if not single coding exon

	if ( $exonToAdd->contains($self->get5PrimeCoding()) ) { #case where the coding st (but not end) is contained in the exon
	  $inCoding=1;
	  my $length;

	  $currentFrame=0;

	  if ($self->getStrand() eq "+") {
	    $length=abs($exonToAdd->getEnd()-$self->get5PrimeCoding())+1;
	  } else {
	    $length=abs($exonToAdd->getStart()-$self->get5PrimeCoding())+1;
	  }

	  push(@{$self->{_framesCalculated}},$currentFrame);

	  my $mod3length=$length%3;


	  $currentFrame=computeNextFrame($currentFrame,$mod3length);

	} elsif ( $exonToAdd->contains($self->get3PrimeCoding()) ) { #case where the coding end (but not st) is contained in the exon
	  my $length;

	  if ($self->getStrand() eq "+") {
	    $length=abs($self->get3PrimeCoding() - $exonToAdd->getStart() )+1;
	  } else {
	    $length=abs($self->get3PrimeCoding() - $exonToAdd->getEnd() )+1;
	  }

	  my $mod3length=$length%3;

	  push(@{$self->{_framesCalculated}},$currentFrame);

	  $inCoding=0;
	  $currentFrame=-1;
	} else {	  #if there no coding start or end in the exon

	  if ($inCoding) {
	    my $length;

	    if ($self->getStrand() eq "+") {
	      $length=abs( $exonToAdd->getEnd() - $exonToAdd->getStart() )+1;
	    } else {
	      $length=abs( $exonToAdd->getEnd() - $exonToAdd->getStart() )+1;
	    }
	    my $mod3length=$length%3;

	    push(@{$self->{_framesCalculated}},$currentFrame);
	    $currentFrame=computeNextFrame($currentFrame,$mod3length);


	  } else {

	    #check current frame
	    push(@{$self->{_framesCalculated}},$currentFrame);
	  }


	}	      #end if there no coding start or end in the exon

      }				#end if not single coding exon

    }				#end for

    my $frameProblem=0;
    my @framesWithProblems;


    for (my $exonToAddIndex=0 ; $exonToAddIndex<$self->getExonCount() ; $exonToAddIndex++) {
      if ($self->{_framesFromFile}) {
	if ($self->{_framesFromFile}->[$exonToAddIndex] != $self->{_framesCalculated}->[$exonToAddIndex]) {
	  $frameProblem=1;
	  push(@framesWithProblems, ($exonToAddIndex+1) );
	}
      }
    }

    if ($frameProblem) {
      if ($#framesWithProblems == 0) {
	warn "Error, CodingTranscript.pm The computed frame for exon #".$framesWithProblems[0]." in  ".$self->getId()." is not the same as the one on file\n";
      } else {
	warn "Error, CodingTranscript.pm The computed frames for exons #".join(",",@framesWithProblems)." in  ".$self->getId()." are not the same as the one on file\n";
      }
      $self->{_isReliable}=0;
    }

    #END DETECTING FRAMES











































    #BEGIN EXTRACTING THE EXONIC SEQUENCES
    $inCoding=0;

    #EXTRACT THE ENTIRE SEQUENCE
    my $sequenceTranscript=$self->extractSequence($self->get5Prime(),$self->get3Prime(),$self->getChromosome());

    for (my $exonToAddIndex=0 ; $exonToAddIndex<$self->getExonCount() ; $exonToAddIndex++) {
      my $exonToAdd=$self->getExons()->[$exonToAddIndex];


      if ($exonToAdd->contains($self->get5PrimeCoding()) && $exonToAdd->contains($self->get3PrimeCoding())) { #case where the coding st and end is contained in a single exon

	my $codingEx;
	my $sequenceExon;

	if ($self->getStrand() eq "+") {
	  $codingEx=CodingExon->new(start         =>  $self->get5PrimeCoding(),
				    end           =>  $self->get3PrimeCoding(),
				    frame         =>  $self->{_framesFromFile}->[$exonToAddIndex], #$currentFrame,
				    strand        =>  $self->getStrand(),
				    indexExon     =>  $codingExonIndex++);
	  $sequenceExon=substr($sequenceTranscript, ( $self->get5PrimeCoding()-$self->get5Prime() ), ($self->get3PrimeCoding()-$self->get5PrimeCoding()+1) );
	} else {
	  $codingEx=CodingExon->new(start         =>  $self->get3PrimeCoding(),
				    end           =>  $self->get5PrimeCoding(),
				    frame         =>  $self->{_framesFromFile}->[$exonToAddIndex], #$currentFrame,
				    strand        =>  $self->getStrand(),
				    indexExon     =>  $codingExonIndex++);
	  $sequenceExon=substr($sequenceTranscript, ( $self->get5Prime()-$self->get5PrimeCoding() ), ($self->get5PrimeCoding()-$self->get3PrimeCoding()+1) );
	}

	push(@{$self->{_codingExons}},$codingEx);


	$codingEx->setGenomicSequence($sequenceExon);
      } else {			#if not single coding exon

	if ( $exonToAdd->contains($self->get5PrimeCoding()) ) { #case where the coding st (but not end) is contained in the exon

	  my $sequenceExon;
	  my $codingEx;

	  if ($self->getStrand() eq "+") {
	    $sequenceExon=substr($sequenceTranscript, ( $self->get5PrimeCoding()-$self->get5Prime() ), ($exonToAdd->getEnd()-$self->get5PrimeCoding()+1) );

	    $codingEx=CodingExon->new(start         =>  $self->get5PrimeCoding(),
				      end           =>  $exonToAdd->getEnd(),
				      frame         =>  $self->{_framesFromFile}->[$exonToAddIndex], #$currentFrame,
				      strand        =>  $self->getStrand(),
				      indexExon     =>  $codingExonIndex++);
	  } else {
	    $sequenceExon=substr($sequenceTranscript, ( $self->get5Prime()-$self->get5PrimeCoding() ), ($self->get5PrimeCoding()-$exonToAdd->getStart()+1) );
	    $codingEx=CodingExon->new(start         =>  $exonToAdd->getStart(),
				      end           =>  $self->get5PrimeCoding(),
				      frame         =>  $self->{_framesFromFile}->[$exonToAddIndex], #$currentFrame,
				      strand        =>  $self->getStrand(),
				      indexExon     =>  $codingExonIndex++);
	  }
	  push(@{$self->{_codingExons}},$codingEx);
	  $inCoding=1;
	  #store current frame and compute next


	  $codingEx->setGenomicSequence($sequenceExon);



	} elsif ( $exonToAdd->contains($self->get3PrimeCoding()) ) { #case where the coding end (but not st) is contained in the exon
	  my $sequenceExon;
	  my $codingEx;


	  if ($self->getStrand() eq "+") {
	    $sequenceExon=substr($sequenceTranscript, ( $exonToAdd->getStart()-$self->get5Prime() ), ($self->get3PrimeCoding()-$exonToAdd->getStart()+1) );
	    $codingEx=CodingExon->new(start         =>  $exonToAdd->getStart(),
				      end           =>  $self->get3PrimeCoding(),
				      frame         =>  $self->{_framesFromFile}->[$exonToAddIndex], #$currentFrame,
				      strand        =>  $self->getStrand(),
				      indexExon     =>  $codingExonIndex);
	  } else {
	    $sequenceExon=substr($sequenceTranscript, ( $self->get5Prime()-$exonToAdd->getEnd() ), ($exonToAdd->getEnd()-$self->get3PrimeCoding()+1) );

	    $codingEx=CodingExon->new(start         =>  $self->get3PrimeCoding(),
				      end           =>  $exonToAdd->getEnd(),
				      frame         =>  $self->{_framesFromFile}->[$exonToAddIndex], #$currentFrame,
				      strand        =>  $self->getStrand(),
				      indexExon     =>  $codingExonIndex);

	  }
	  push(@{$self->{_codingExons}},$codingEx);
	  #print "sequenceExon ".$sequenceExon."\n";

	  $codingEx->setGenomicSequence($sequenceExon);


	  $inCoding=0;
	  $currentFrame=-1;
	} else {	  #if there no coding start or end in the exon

	  if ($inCoding) {
	    my $sequenceExon;
	    my $codingEx;



	    if ($self->getStrand() eq "+") {
	      $sequenceExon=substr($sequenceTranscript, ( $exonToAdd->getStart()-$self->get5Prime() ), ($exonToAdd->getEnd()-$exonToAdd->getStart()+1) );

	      $codingEx=CodingExon->new(start         =>  $exonToAdd->getStart(),
					end           =>  $exonToAdd->getEnd(),
					frame         =>  $self->{_framesFromFile}->[$exonToAddIndex], #$currentFrame,
					strand        =>  $self->getStrand(),
					indexExon     =>  $codingExonIndex++);
	    } else {
	      $sequenceExon=substr($sequenceTranscript, ( $self->get5Prime()-$exonToAdd->getEnd() ), ($exonToAdd->getEnd()-$exonToAdd->getStart()+1) );
	      $codingEx=CodingExon->new(start         =>  $exonToAdd->getStart(),
					end           =>  $exonToAdd->getEnd(),
					frame         =>  $self->{_framesFromFile}->[$exonToAddIndex], #$currentFrame,
					strand        =>  $self->getStrand(),
					indexExon     =>  $codingExonIndex++);
	    }
	    push(@{$self->{_codingExons}},$codingEx);



	    $codingEx->setGenomicSequence($sequenceExon);



	  } else {
	    my $codingEx;
	    if ($self->getStrand() eq "+") {
	      $codingEx=CodingExon->new(start         =>  $exonToAdd->getStart(),
					end           =>  $exonToAdd->getEnd(),
					frame         =>  $self->{_framesFromFile}->[$exonToAddIndex], #$currentFrame,
					strand        =>  $self->getStrand(),
					indexExon     =>  -1);
	    } else {
	      $codingEx=CodingExon->new(start         =>  $exonToAdd->getStart(),
					end           =>  $exonToAdd->getEnd(),
					frame         =>  $self->{_framesFromFile}->[$exonToAddIndex], #$currentFrame,
					strand        =>  $self->getStrand(),
					indexExon     =>  -1);
	    }
	    push(@{$self->{_codingExons}},$codingEx);



	  }


	}	      #end if there no coding start or end in the exon

      }				#end if not single coding exon

    }				#end for

    #END EXTRACTING THE EXONIC SEQUENCES

































    #BEGIN FIXING THE GENOMIC SEQUENCES
    my @exonsNoMult3;
    my $exonsWithNoMult3=0;

    for (my $exonToAddIndex=0 ; $exonToAddIndex<$self->getExonCount() ; $exonToAddIndex++) {

      if ($self->{_codingExons}->[$exonToAddIndex]->getIndexExon() == 0) { #first coding exon

	if ( ( $exonToAddIndex != ($self->getExonCount()-1) ) &&
	     $self->{_codingExons}->[$exonToAddIndex+1]->isCoding() ) {

	  if ($self->{_codingExons}->[$exonToAddIndex+1]->getFrame() == 1) { #add 2 at the end
	    $self->{_codingExons}->[$exonToAddIndex]->appendEndGenomicSequence( $self->{_codingExons}->[$exonToAddIndex+1]->getFirst2CharOfOrigSeq() );
	  }

	  if ($self->{_codingExons}->[$exonToAddIndex+1]->getFrame() == 2) { #add 1 at the end
	    $self->{_codingExons}->[$exonToAddIndex]->appendEndGenomicSequence( $self->{_codingExons}->[$exonToAddIndex+1]->getFirstCharOfOrigSeq() );
	  }

	}

      } elsif ($self->{_codingExons}->[$exonToAddIndex]->getIndexExon() == $codingExonIndex ) { #last coding exon

	if ( ( $exonToAddIndex > 0 ) &&
	     $self->{_codingExons}->[$exonToAddIndex-1]->isCoding() ) {

	  if ($self->{_codingExons}->[$exonToAddIndex]->getFrame() == 1) { #add 1 at the beginning
	    $self->{_codingExons}->[$exonToAddIndex]->appendStartGenomicSequence( $self->{_codingExons}->[$exonToAddIndex-1]->getLastCharOfOrigSeq() );
	  }

	  if ($self->{_codingExons}->[$exonToAddIndex]->getFrame() == 2) { #add 2 at the beginning
	    $self->{_codingExons}->[$exonToAddIndex]->appendStartGenomicSequence( $self->{_codingExons}->[$exonToAddIndex-1]->getLast2CharOfOrigSeq() );
	  }

	}


      } else {
	if ( $self->{_codingExons}->[$exonToAddIndex]->isCoding()  ) { #coding exon

	  if ( ( $exonToAddIndex != ($self->getExonCount()-1) ) &&
	       $self->{_codingExons}->[$exonToAddIndex+1]->isCoding() ) {

	    if ($self->{_codingExons}->[$exonToAddIndex+1]->getFrame() == 1) { #add 2 at the end
	      $self->{_codingExons}->[$exonToAddIndex]->appendEndGenomicSequence( $self->{_codingExons}->[$exonToAddIndex+1]->getFirst2CharOfOrigSeq() );
	    }

	    if ($self->{_codingExons}->[$exonToAddIndex+1]->getFrame() == 2) { #add 1 at the end
	      $self->{_codingExons}->[$exonToAddIndex]->appendEndGenomicSequence( $self->{_codingExons}->[$exonToAddIndex+1]->getFirstCharOfOrigSeq() );
	    }

	  }


	  if ( ( $exonToAddIndex > 0 ) &&
	       $self->{_codingExons}->[$exonToAddIndex-1]->isCoding() ) {

	    if ($self->{_codingExons}->[$exonToAddIndex]->getFrame() == 1) { #add 1 at the beginning
	      $self->{_codingExons}->[$exonToAddIndex]->appendStartGenomicSequence( $self->{_codingExons}->[$exonToAddIndex-1]->getLastCharOfOrigSeq() );
	    }

	    if ($self->{_codingExons}->[$exonToAddIndex]->getFrame() == 2) { #add 2 at the beginning
	      $self->{_codingExons}->[$exonToAddIndex]->appendStartGenomicSequence( $self->{_codingExons}->[$exonToAddIndex-1]->getLast2CharOfOrigSeq() );
	    }

	  }



	}
      }


      if ($self->{_codingExons}->[$exonToAddIndex]->isCoding()) {

	if ($self->{_codingExons}->[$exonToAddIndex]->getIndexExon() == 0) { #first coding exon
	  if (!$self->{_codingExons}->[$exonToAddIndex]->startsWithSTARTCodon()) {
	    warn "Error CodingTranscript.pm The first coding exon #".($exonToAddIndex+1)." in  ".$self->getId()." does not start with a START signal\n";
	    $self->{_isReliable}=0;
	    $self->{_codingExons}->[$exonToAddIndex]->setAsUnreliable();
	  }
	}
	if ($self->{_codingExons}->[$exonToAddIndex]->getIndexExon() == $codingExonIndex ) { #last coding exon
	  if($self->{_codingExons}->[$exonToAddIndex]->getGenomicSequenceLength() >= 3){
	    if (!$self->{_codingExons}->[$exonToAddIndex]->endsWithSTOPCodon()) {
	      warn "Error CodingTranscript.pm The last coding exon #".($exonToAddIndex+1)." in  ".$self->getId()." does not end with a STOP signal\n";
	      $self->{_isReliable}=0;
	      $self->{_codingExons}->[$exonToAddIndex]->setAsUnreliable();
	    }
	  }else{
	    warn "Error CodingTranscript.pm The last coding exon #".($exonToAddIndex+1)." in  ".$self->getId()." is not at least 3 bp long\n";
	    $self->{_isReliable}=0;
	    $self->{_codingExons}->[$exonToAddIndex]->setAsUnreliable();
	  }
	}

	if ( ($self->{_codingExons}->[$exonToAddIndex]->getGenomicSequenceLength() % 3) != 0) {
	  $self->{_isReliable}=0;
	  $self->{_codingExons}->[$exonToAddIndex]->setAsUnreliable();
	  push(@exonsNoMult3,($exonToAddIndex+1));
	  $exonsWithNoMult3=1;


	}
      }

      if ($self->{_framesFromFile}->[$exonToAddIndex] != $self->{_framesCalculated}->[$exonToAddIndex]) {
	$self->{_codingExons}->[$exonToAddIndex]->setAsUnreliable();
      }

    }				#end for

    if ($exonsWithNoMult3) {
      if ($#exonsNoMult3 == 0 ) {
	warn "Error CodingTranscript.pm The length of the coding sequence of exon #".$exonsNoMult3[0]." in  ".$self->getId()." is not a multiple of 3\n";
      } else {
	warn "Error CodingTranscript.pm The length of the coding sequence of exons #".join(",",@exonsNoMult3)." in  ".$self->getId()." are not a multiple of 3\n";
      }
    }
    #END FIXING THE GENOMIC SEQUENCES



  } #end if transcript is coding

}#end sub buildExonCoding



=item annotateSNP

  Subroutine to annotate a SNP that falls within the range
  of the CodingTranscript. It starts by calling annotateCoordinate()
  from Transcript.pm. It the coordinate falls within the coding part of
  the gene.
  Arguments:
   coordinate          : The coordinate of the SNP on the chromosome
   refDNA              : The DNA base pair from the reference sequence
   readDNA             : The DNA base pair from the read
   produceFullProtSeq  : Flag to determine whether or not we have to produce the full protein sequence
   noRefBP             : If set to 1, we will not check to see if the refDNA corresponds to our own reference sequence.
  Returns:
   Returns a hash with the following fields:
     key        value meaning
      'code'     1     exon
                 2     intron
                 3     upstream
                 4     downstream
      'partial'  #     flag to determine if the 5' or 3'
                 #     end is partial for up/downstream cases
      'index'    #     index of the exon or intron
                       (only used when code = 1 or 2)
      'distance5P'    #     distance to the 5 prime end
      'distance3P'    #     distance to the 3 prime end
      'reliable'    1 or 0  Whether the CodingTranscript is reliable or not
                        (see documentation of subroutine isReliable())
      'coding'      1 or 0  Whether the coordinate is within the
                            protein coding range
      'utr'         # if 'code'=1 and coding=0, utr=5 indicated 5p utr, utr=3 indicated 3p utr,
      'splice'      # if 'code'=2, if the SNP is located within 2 bp of an splice junction
      'codonRef'    # The codon pertaining to the reference sequence
      'codonRead'   # The codon pertaining to the read
      'aaRef'       # The amino acid pertaining to the reference sequence
      'aaRead'      # The amino acid pertaining to the read
      'synonymous'  # 1 if aaRef aaRead are identical, 0 otherwise
      'seqRef'      # The full protein sequence for the reference sequence (if produceFullProtSeq is set to 1)
      'seqRead'     # The full protein sequence for the read  (if produceFullProtSeq is set to 1)
      'coordInAASeq' # Coordinate of SNP within the amino acid sequence


=cut
sub annotateSNP{
  my($self,$coordinate,$refDNA,$readDNA,$produceFullProtSeq,$noRefBP)=@_;

  if(!$self->isCoding()){
    my $resultsToReturn=$self->SUPER::annotateCoordinate($coordinate);
    if ( $self->{_isReliable} ) { #if the current transcript is reliable
      $resultsToReturn->{'reliable'}=1;
    }else{
      $resultsToReturn->{'reliable'}=0;
    }
    return $resultsToReturn;
  }

  if( !($self->{_codingExons}) ){
    die "Error CodingTranscript.pm Please call the buildExonCoding() subroutine prior to annotateSNP()\n";
  }

  my $resultsToReturn=$self->SUPER::annotateCoordinate($coordinate);

  if ( $self->{_isReliable} ) { #if the current transcript is reliable
    $resultsToReturn->{'reliable'}=1;
    if ( ($resultsToReturn->{'code'} == 1) ) { #is in an exon

      if ( $self->isCoordinateInCodingRange($coordinate) ) {
	$resultsToReturn->{'coding'}=1;


	if($produceFullProtSeq){ #if we need to produce the protein sequence
	  my $proteinSequenceRead="";
	  my $proteinSequenceRef="";


	  for(my $exonToAddIndex=0 ; $exonToAddIndex<$self->getExonCount() ; $exonToAddIndex++){
	    if($exonToAddIndex == ($resultsToReturn->{'index'}-1)){
	      my ($codonRef,$codonRead,$aaRef,$aaRead,$coordInAAExon,$seqRef,$seqRead) = $self->{_codingExons}->[$resultsToReturn->{'index'}-1]->annotateBasePair($coordinate,$refDNA,$readDNA,1,$noRefBP);


	      $resultsToReturn->{'coordInAASeq'}    = length($proteinSequenceRead) + $coordInAAExon ;


	      if($self->{_codingExons}->[$exonToAddIndex]->getBasePairsAtTheStart() != 0){ #if an aa was added at the beginning, we need to remove it
		$resultsToReturn->{'coordInAASeq'}--;
		$proteinSequenceRead .= substr($seqRead,1);
		$proteinSequenceRef  .= substr($seqRef,1);
	      }else{
		$proteinSequenceRead .= $seqRead;
		$proteinSequenceRef  .= $seqRef;
	      }


	      $resultsToReturn->{'codonRef'}        = $codonRef;
	      $resultsToReturn->{'codonRead'}       = $codonRead;
	      $resultsToReturn->{'aaRef'}           = $aaRef;
	      $resultsToReturn->{'aaRead'}          = $aaRead;

	      if($aaRef eq $aaRead){
		$resultsToReturn->{'synonymous'}=1;
	      }else{
		$resultsToReturn->{'synonymous'}=0;
	      }


	    }else{

	      $self->{_codingExons}->[$exonToAddIndex]->translate2Protein();
	      if($self->{_codingExons}->[$exonToAddIndex]->getBasePairsAtTheStart() != 0){ #if an aa was added at the beginning, we need to remove it
		$proteinSequenceRead .= substr($self->{_codingExons}->[$exonToAddIndex]->getProteomicSequence(),1);
		$proteinSequenceRef  .= substr($self->{_codingExons}->[$exonToAddIndex]->getProteomicSequence(),1);

	      }else{
		$proteinSequenceRead .= $self->{_codingExons}->[$exonToAddIndex]->getProteomicSequence();
		$proteinSequenceRef  .= $self->{_codingExons}->[$exonToAddIndex]->getProteomicSequence();
	      }
	    }
	  } #end for each exon

	  $resultsToReturn->{'seqRef'}     = $proteinSequenceRef;
	  $resultsToReturn->{'seqRead'}    = $proteinSequenceRead;



	}else{ #if we do not need to produce the protein sequence
	  my ($codonRef,$codonRead,$aaRef,$aaRead) = $self->{_codingExons}->[$resultsToReturn->{'index'}-1]->annotateBasePair($coordinate,$refDNA,$readDNA,0,$noRefBP);
	  $resultsToReturn->{'codonRef'}   = $codonRef;
	  $resultsToReturn->{'codonRead'}  = $codonRead;
	  $resultsToReturn->{'aaRef'}      = $aaRef;
	  $resultsToReturn->{'aaRead'}     = $aaRead;
	  $resultsToReturn->{'seqRef'}     = 0;
	  $resultsToReturn->{'seqRead'}    = 0;

	  if($aaRef eq $aaRead){
	    $resultsToReturn->{'synonymous'}=1;
	  }else{
	    $resultsToReturn->{'synonymous'}=0;
	  }

	}
      } else { #not coding
	$resultsToReturn->{'coding'}=0;
	if($self->getStrand() eq "+"){

	  if( ($self->get5Prime() <= $coordinate) &&
	      ($coordinate        <= $self->get5PrimeCoding()) ){
	    $resultsToReturn->{'utr'}=5;
	  }

	  if( ($self->get3PrimeCoding() <= $coordinate) &&
	      ($coordinate              <= $self->get3Prime()) ){
	    $resultsToReturn->{'utr'}=3;
	  }

	}elsif($self->getStrand() eq "-"){

	  if( ($self->get5PrimeCoding() <= $coordinate) &&
	      ($coordinate              <= $self->get5Prime() ) ){
	    $resultsToReturn->{'utr'}=5;
	  }

	  if( ($self->get3Prime() <= $coordinate) &&
	      ($coordinate        <= $self->get3PrimeCoding() ) ){
	    $resultsToReturn->{'utr'}=3;
	  }

	}else{
	  die "Strand ".$self->getStrand()." is wrong in ".Dumper($self)."\n";
	}
      }

    }elsif ( ($resultsToReturn->{'code'} == 2) ) { #is in an intron
      #'splice'
      my $intronToCheck=$self->getIntrons()->[$resultsToReturn->{'index'}-1];
      if($intronToCheck->contains($coordinate)){
	if( ($coordinate-$intronToCheck->getStart()) < 2 ||
	    ($intronToCheck->getEnd()-$coordinate) < 2 ){
	  $resultsToReturn->{'splice'}=1;
	}else{
	  $resultsToReturn->{'splice'}=0;
	}
      }else{
	die "The intron index returned is wrong in ".Dumper($self)."\n";
      }
    }else{
      $resultsToReturn->{'splice'}=0;
    }
  } else {
    $resultsToReturn->{'reliable'}=0;
    if ( $self->isCoordinateInCodingRange($coordinate) ) {
      $resultsToReturn->{'coding'}=1;
    } else {
      $resultsToReturn->{'coding'}=0;
    }
  }



  return $resultsToReturn;
}









=item annotateINSERT

  Subroutine to annotate an insert that falls within the range
  of the CodingTranscript. It starts by calling annotateCoordinate()
  from Transcript.pm. It the coordinate falls within the coding part of
  the gene.
  Arguments:
   coordinate          : The coordinate of the insert on the chromosome
   insertSeq           : The DNA sequence of the insert
  Returns:
   Returns a hash with the following fields:
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
      'reliable'    1 or 0  Whether the CodingTranscript is reliable or not
                        (see documentation of subroutine isReliable())
      'coding'      1 or 0  Whether the coordinate is within the
                            protein coding range
      'seqRef'      # The full protein sequence for the reference sequence (if produceFullProtSeq is set to 1)
      'seqRead'     # The full protein sequence for the read with the insert (if produceFullProtSeq is set to 1)
      'coordInAASeq' # Coordinate of the insert within the amino acid sequence


=cut
sub annotateINSERT{
  my($self,$coordinate,$insertSeq)=@_;

  $insertSeq=uc($insertSeq);

  if($self->{_strand} eq "+"){
    foreach my $char (split("",$insertSeq)){
      if( !(exists $complement{$char}) ){
	die "CodingTranscript.pm annotateINSERT() DNA sequence must be either A,C,G,T\n";
      }
    }
  }else{ #take the reverse comp of the string if 
    my $tempS="";
    foreach my $char (reverse(split("",$insertSeq))){
      if( !(exists $complement{$char}) ){
	die "CodingTranscript.pm annotateINSERT() DNA sequence must be either A,C,G,T\n";
      }else{
	$tempS.=$complement{$char};
      }
    }
    $insertSeq=$tempS;
  }

  if(!$self->isCoding()){
    my $resultsToReturn=$self->SUPER::annotateCoordinate($coordinate);
    if ( $self->{_isReliable} ) { #if the current transcript is reliable
      $resultsToReturn->{'reliable'}=1;
    }else{
      $resultsToReturn->{'reliable'}=0;
    }
    return $resultsToReturn;
  }

  if( !($self->{_codingExons}) ){
    die "Error CodingTranscript.pm Please call the buildExonCoding() subroutine prior to annotateSNP()\n";
  }

  my $resultsToReturn=$self->SUPER::annotateCoordinate($coordinate);


  if ( $self->{_isReliable} ) { #if the current transcript is reliable
    $resultsToReturn->{'reliable'}=1;
    if ( ($resultsToReturn->{'code'} == 1) ) { #is in an exon

      if ( $self->isCoordinateInCodingRange($coordinate) ) {

	$resultsToReturn->{'coding'}=1;

	my $proteinSequenceRead="";
	my $proteinSequenceRef="";
	my $bpToAppend="";

	for (my $exonToAddIndex=0 ; $exonToAddIndex<$self->getExonCount() ; $exonToAddIndex++) {
	  if ($exonToAddIndex == ($resultsToReturn->{'index'}-1)) {

	    my ($coordInAAExon,$seqRef,$seqRead,$remainingBP) = $self->{_codingExons}->[$resultsToReturn->{'index'}-1]->annotateInsert($coordinate,$insertSeq);
	    $bpToAppend=$remainingBP;

	    $resultsToReturn->{'coordInAASeq'}    = length($proteinSequenceRead) + $coordInAAExon ;


	    if ($self->{_codingExons}->[$exonToAddIndex]->getBasePairsAtTheStart() != 0) { #if an aa was added at the beginning, we need to remove it
	      $resultsToReturn->{'coordInAASeq'}--;
	      $proteinSequenceRead .= substr($seqRead,1);
	      $proteinSequenceRef  .= substr($seqRef,1);
	    } else {
	      $proteinSequenceRead .= $seqRead;
	      $proteinSequenceRef  .= $seqRef;
	    }



	  } elsif ($exonToAddIndex < ($resultsToReturn->{'index'}-1)) {
	    $self->{_codingExons}->[$exonToAddIndex]->translate2Protein();
	    if ($self->{_codingExons}->[$exonToAddIndex]->getBasePairsAtTheStart() != 0) { #if an aa was added at the beginning, we need to remove it
	      $proteinSequenceRead .= substr($self->{_codingExons}->[$exonToAddIndex]->getProteomicSequence(),1);
	      $proteinSequenceRef  .= substr($self->{_codingExons}->[$exonToAddIndex]->getProteomicSequence(),1);
	    } else {
	      $proteinSequenceRead .= $self->{_codingExons}->[$exonToAddIndex]->getProteomicSequence();
	      $proteinSequenceRef  .= $self->{_codingExons}->[$exonToAddIndex]->getProteomicSequence();
	    }

	  }else{

	    my ($protseq,$remainingBP)=$self->{_codingExons}->[$exonToAddIndex]->translate2ProteinCustom($bpToAppend);
	    $bpToAppend=$remainingBP;
	    $self->{_codingExons}->[$exonToAddIndex]->translate2Protein();

	    if ($self->{_codingExons}->[$exonToAddIndex]->getBasePairsAtTheStart() != 0) { #if an aa was added at the beginning, we need to remove it
	      $proteinSequenceRead .= $protseq;
	      $proteinSequenceRef  .= substr($self->{_codingExons}->[$exonToAddIndex]->getProteomicSequence(),1);
	    } else {
	      $proteinSequenceRead .= $protseq;
	      $proteinSequenceRef  .= $self->{_codingExons}->[$exonToAddIndex]->getProteomicSequence();
	    }

	  }
	}	#end for each exon

	$resultsToReturn->{'seqRef'}     = $proteinSequenceRef;
	$resultsToReturn->{'seqRead'}    = $proteinSequenceRead;

      } else {
	$resultsToReturn->{'coding'}=0;
      }

    }
  } else {
    $resultsToReturn->{'reliable'}=0;
    if ( $self->isCoordinateInCodingRange($coordinate) ) {
      $resultsToReturn->{'coding'}=1;
    } else {
      $resultsToReturn->{'coding'}=0;
    }
  }



  return $resultsToReturn;
}





=item annotateDELETION

  Subroutine to annotate a deletion that falls within the range
  of the CodingTranscript. It starts by calling annotateCoordinate()
  from Transcript.pm on coordinate1 and coordinate2. It then goes through
  the exons to determine if they are:
    1) Downstream of both coordinates
    2) Upstream   of both coordinates
    3) Completely encompassed within the deletion
    4) One coordinate within the exon and the other is upstream
    5) One coordinate within the exon and the other is downstream
    6) Both coordinate are contained within the exon
  It will return 2 sequences, the original protein sequence and the
  one modified by the deletion
  Arguments:
   coordinate1         : The lowest coordinate on the chromosome of the deletion
   coordinate2         : The highest coordinate on the chromosome of the deletion
  Returns:
   Returns a hash with the following fields:
     key        value meaning
      'coding'      1 or 0  Whether the coordinate is within the
                            protein coding range
      'seqRef'      # The full protein sequence for the reference sequence (if produceFullProtSeq is set to 1)
      'seqRead'     # The full protein sequence for the read with the insert (if produceFullProtSeq is set to 1)


=cut
sub annotateDELETION{
  my($self,$coordinate1,$coordinate2)=@_;

  my $deletionRange=Range->new(start => $coordinate1,
			       end   => $coordinate2);

  my $resultsToReturn;
  my $coordinate5p;
  my $coordinate3p;
  if ($self->getStrand() eq "+") {
    $coordinate5p=$coordinate1;
    $coordinate3p=$coordinate2;
  }elsif ($self->getStrand() eq "-") {
    $coordinate5p=$coordinate2;
    $coordinate3p=$coordinate1;
  }

  my $results1=$self->SUPER::annotateCoordinate($coordinate5p);
  my $results2=$self->SUPER::annotateCoordinate($coordinate3p);

  if ( $self->{_isReliable} ) { #if the current transcript is reliable
    $resultsToReturn->{'reliable'}=1;

    if ( $self->isRangeInCodingRange($deletionRange)  ) {
      $resultsToReturn->{'coding'}=1;

      my $proteinSequenceRead ="";
      my $proteinSequenceRef  ="";
      my $bpToAppendRead      ="";
      my $bpToAppendRef       ="";

      for (my $exonToAddIndex=0 ; $exonToAddIndex<$self->getExonCount() ; $exonToAddIndex++) {

	if(!$self->{_codingExons}->[$exonToAddIndex]->isCoding()){
	  next;
	}

	if (       $self->{_codingExons}->[$exonToAddIndex]->isCoordDownstream($coordinate5p) &&
		   $self->{_codingExons}->[$exonToAddIndex]->isCoordDownstream($coordinate3p) ) {
	  my ($protseqread,$remainingBPRead) = $self->{_codingExons}->[$exonToAddIndex]->translate2ProteinCustom($bpToAppendRead);
	  my ($protseqref ,$remainingBPRef)  = $self->{_codingExons}->[$exonToAddIndex]->translate2ProteinCustom($bpToAppendRef);
	  $bpToAppendRead = $remainingBPRead;
	  $bpToAppendRef  = $remainingBPRef;
	  $proteinSequenceRead .= $protseqread;
	  $proteinSequenceRef  .= $protseqref;

	} elsif ( $self->{_codingExons}->[$exonToAddIndex]->isCoordUpstream($coordinate5p) &&
		  $self->{_codingExons}->[$exonToAddIndex]->isCoordUpstream($coordinate3p) ) {
	  my ($protseqread,$remainingBPRead) = $self->{_codingExons}->[$exonToAddIndex]->translate2ProteinCustom($bpToAppendRead);
	  my ($protseqref,$remainingBPRef)   = $self->{_codingExons}->[$exonToAddIndex]->translate2ProteinCustom($bpToAppendRef);
	  $bpToAppendRead = $remainingBPRead;
	  $bpToAppendRef  = $remainingBPRef;
	  $proteinSequenceRead .= $protseqread;
	  $proteinSequenceRef  .= $protseqref;


	} elsif ( $self->{_codingExons}->[$exonToAddIndex]->isCoordUpstream($coordinate5p)   &&
		  $self->{_codingExons}->[$exonToAddIndex]->isCoordDownstream($coordinate3p) ) { #completely encompassed
	  my ($protseqref ,$remainingBPRef)  = $self->{_codingExons}->[$exonToAddIndex]->translate2ProteinCustom($bpToAppendRef);
	  $bpToAppendRef  = $remainingBPRef;
	  $proteinSequenceRef  .= $protseqref;

	} elsif (   $self->{_codingExons}->[$exonToAddIndex]->isCoordUpstream($coordinate5p)   &&
		    $self->{_codingExons}->[$exonToAddIndex]->isCoordContained($coordinate3p) ) {
	  my ($protseqread,$remainingBPRead) = $self->{_codingExons}->[$exonToAddIndex]->annotateDeletion($coordinate3p,$self->{_codingExons}->[$exonToAddIndex]->get5PrimeEnd(),$bpToAppendRead);
	  my ($protseqref ,$remainingBPRef)  = $self->{_codingExons}->[$exonToAddIndex]->translate2ProteinCustom($bpToAppendRef);

	  $bpToAppendRead = $remainingBPRead;
	  $bpToAppendRef  = $remainingBPRef;

	  $proteinSequenceRead .= $protseqread;
	  $proteinSequenceRef  .= $protseqref;

	} elsif (   $self->{_codingExons}->[$exonToAddIndex]->isCoordContained($coordinate5p)   &&
		    $self->{_codingExons}->[$exonToAddIndex]->isCoordDownstream($coordinate3p) ) {
	  my ($protseqread,$remainingBPRead) = $self->{_codingExons}->[$exonToAddIndex]->annotateDeletion($coordinate5p,
													  $self->{_codingExons}->[$exonToAddIndex]->get3PrimeEnd(),
													  $bpToAppendRead);
	  my ($protseqref ,$remainingBPRef)  = $self->{_codingExons}->[$exonToAddIndex]->translate2ProteinCustom($bpToAppendRef);

	  $bpToAppendRead = $remainingBPRead;
	  $bpToAppendRef  = $remainingBPRef;

	  $proteinSequenceRead .= $protseqread;
	  $proteinSequenceRef  .= $protseqref;

	} elsif (   $self->{_codingExons}->[$exonToAddIndex]->isCoordContained($coordinate5p)   &&
		    $self->{_codingExons}->[$exonToAddIndex]->isCoordContained($coordinate3p) ) {

	  my ($protseqread,$remainingBPRead) = $self->{_codingExons}->[$exonToAddIndex]->annotateDeletion($coordinate1,$coordinate2,$bpToAppendRead);
	  my ($protseqref ,$remainingBPRef)  = $self->{_codingExons}->[$exonToAddIndex]->translate2ProteinCustom($bpToAppendRef);

	  $bpToAppendRead = $remainingBPRead;
	  $bpToAppendRef  = $remainingBPRef;

	  $proteinSequenceRead .= $protseqread;
	  $proteinSequenceRef  .= $protseqref;

	} else {
	  die "CodingTranscript.pm annotateDELETION() wrong state for coordinates $coordinate5p - $coordinate3p ".Dumper($self)."\n";
	}
      }				#end for each exon

      $resultsToReturn->{'seqRef'}     = $proteinSequenceRef;
      $resultsToReturn->{'seqRead'}    = $proteinSequenceRead;
    } else {
      $resultsToReturn->{'coding'}=0;
    }

  } else {
    $resultsToReturn->{'reliable'}=0;
    if ( $self->isRangeInCodingRange($deletionRange)  ) {
      $resultsToReturn->{'coding'}=1;
    } else {
      $resultsToReturn->{'coding'}=0;
    }
  }



  return $resultsToReturn;


}



=item codingExonsAvailable

  Subroutine to verify if buildExonCoding()
  has been called.

=cut
sub codingExonsAvailable{
  my($self)=@_;
  if ( !($self->{_codingExons}) ) {
    return 0;
  }else{
    return 1;
  }

}


=item returnProteinSequence

  Subroutine to return the protein sequence associated with the
  CodingTranscript.

=cut
sub returnProteinSequence{
  my($self)=@_;

  my $proteinSequenceRead="";

  if(!$self->isCoding()){
      #not coding
  } else {
    if ( !($self->{_codingExons}) ) {
      die "Error CodingTranscript.pm Please call the buildExonCoding() subroutine prior to annotateSNP()\n";
    }

    if ( $self->{_isReliable} ) { #if the current transcript is reliable

      for (my $exonToAddIndex=0 ; $exonToAddIndex<$self->getExonCount() ; $exonToAddIndex++) {
	$self->{_codingExons}->[$exonToAddIndex]->translate2Protein();
	if ($self->{_codingExons}->[$exonToAddIndex]->getBasePairsAtTheStart() != 0) { #if an aa was added at the beginning, we need to remove it
	  $proteinSequenceRead .= substr($self->{_codingExons}->[$exonToAddIndex]->getProteomicSequence(),1);
	} else {
	  $proteinSequenceRead .= $self->{_codingExons}->[$exonToAddIndex]->getProteomicSequence();
	}
      }
    }
  }


  return $proteinSequenceRead;
}



=item annotateRange

  Subroutine to add the following field:
  'coding'      #     to know whether the coordinate
                      is within the coding range for those
                      landing in the gene
                      code:
                        1 not in coding range
                        2 in coding range
                        3 overlaps a coding start/end


=cut
sub annotateRange{
  my($self,$coordinate1,$coordinate2)=@_;
  my $resultsToReturn = $self->SUPER::annotateRange($coordinate1,$coordinate2);

  if($resultsToReturn->{'code'} == 9 || $resultsToReturn->{'code'} == 10){
    $resultsToReturn->{'coding'}=0;
  }else{

    if(      $self->getCodingRange()->contains($coordinate1) &&  $self->getCodingRange()->contains($coordinate2) ){
      $resultsToReturn->{'coding'}=2;
    }elsif( !$self->getCodingRange()->contains($coordinate1) &&  $self->getCodingRange()->contains($coordinate2) ){
      $resultsToReturn->{'coding'}=3;
    }elsif(  $self->getCodingRange()->contains($coordinate1) && !$self->getCodingRange()->contains($coordinate2) ){
      $resultsToReturn->{'coding'}=3;
    }elsif( !$self->getCodingRange()->contains($coordinate1) && !$self->getCodingRange()->contains($coordinate2) ){
      $resultsToReturn->{'coding'}=1;
    }

  }

  return $resultsToReturn;
}


1;
