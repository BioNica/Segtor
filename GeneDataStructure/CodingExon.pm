=head1 NAME

   GeneDataStructure::CodingExon.pm


=head1 Disclaimer

This software and all the data therein were developed by the Bioinformatics and Computational Biology Lab (Laboratorio de Bioinformatica e Biologia Computational (LBBC)) of the Brazilian National Cancer Institute (Instituto Nacional de Cancer (INCA)) in Rio de Janeiro, RJ, Brazil. The INCA is part of the Health Ministry (Ministerio da Saude) of the Federal Government of Brazil. This software was developed by the staff of the LBBC as part of their official duties. According to federal regulations, this software is therefore public domain.

This software is provided as is and the INCA does make any guarantee regarding the accuracy/quality/reliability of the results. The INCA cannot be held accountable for any hardware or software damage that might be caused through the use or misuse of this software.

The LBBC encourages researchers using our software to acknowledge us. The software can be distributed freely and modified given that any derivative bears an acknowledgment that it was derived from our original software.

=head1 Author

Gabriel Renaud grenaud@inca.gov.br



=cut


package CodingExon;


use strict;
use warnings;

use Data::Dumper;

use POSIX qw(ceil floor);

my $stopSymbol='*';

my %codon2aa=(
	      "AAA"   =>  'K',
	      "AAC"   =>  'N',
	      "AAG"   =>  'K',
	      "AAT"   =>  'N',
	      "ACA"   =>  'T',
	      "ACC"   =>  'T',
	      "ACG"   =>  'T',
	      "ACT"   =>  'T',
	      "AGA"   =>  'R',
	      "AGC"   =>  'S',
	      "AGG"   =>  'R',
	      "AGT"   =>  'S',
	      "ATA"   =>  'I',
	      "ATC"   =>  'I',
	      "ATG"   =>  'M',
	      "ATT"   =>  'I',
	      "CAA"   =>  'Q',
	      "CAC"   =>  'H',
	      "CAG"   =>  'Q',
	      "CAT"   =>  'H',
	      "CCA"   =>  'P',
	      "CCC"   =>  'P',
	      "CCG"   =>  'P',
	      "CCT"   =>  'P',
	      "CGA"   =>  'R',
	      "CGC"   =>  'R',
	      "CGG"   =>  'R',
	      "CGT"   =>  'R',
	      "CTA"   =>  'L',
	      "CTC"   =>  'L',
	      "CTG"   =>  'L',
	      "CTT"   =>  'L',
	      "GAA"   =>  'E',
	      "GAC"   =>  'D',
	      "GAG"   =>  'E',
	      "GAT"   =>  'D',
	      "GCA"   =>  'A',
	      "GCC"   =>  'A',
	      "GCG"   =>  'A',
	      "GCT"   =>  'A',
	      "GGA"   =>  'G',
	      "GGC"   =>  'G',
	      "GGG"   =>  'G',
	      "GGT"   =>  'G',
	      "GTA"   =>  'V',
	      "GTC"   =>  'V',
	      "GTG"   =>  'V',
	      "GTT"   =>  'V',
	      "TAA"   =>  $stopSymbol,
	      "TAC"   =>  'Y',
	      "TAG"   =>  $stopSymbol,
	      "TAT"   =>  'Y',
	      "TCA"   =>  'S',
	      "TCC"   =>  'S',
	      "TCG"   =>  'S',
	      "TCT"   =>  'S',
	      "TGA"   =>  $stopSymbol,
	      "TGC"   =>  'C',
	      "TGG"   =>  'W',
	      "TGT"   =>  'C',
	      "TTA"   =>  'L',
	      "TTC"   =>  'F',
	      "TTG"   =>  'L',
	      "TTT"   =>  'F');

my %dna2RC=('A'=> 'T',
	    'C'=> 'G',
	    'G'=> 'C',
	    'T'=> 'A');

=item new

   This object represents an exon along with its sequence.
   Arguments:
    start           start coordinate on the chromosome
    end             end coordinate on the chromosome
    frame           frame of the exon


=cut
sub new{
  my ($class,%arg)=(@_);

  my ($self) =bless{_start                       => $arg{start},
                    _end                         => $arg{end},
                    _frame                       => $arg{frame},
		    _strand                      => $arg{strand},
		    _indexExon                   => $arg{indexExon}}, $class;

  if($arg{frame} == -1){
    $self->{_isCoding}=0;
  }elsif($arg{frame} >= 0 && $arg{frame} <= 2){
    $self->{_isCoding}=1;
  }else{
    die "Error CodingExon.pm The frame should either be -1,0,1,2\n";
  }

  if($self->{_start} > $self->{_end} ){
    die "Error CodingExon.pm The start ".$self->{_start}." is greater than the end ".$self->{_end}."\n";
  }

  if($self->{_strand} eq "+"){
    $self->{_coord5p}=$self->{_start};
    $self->{_coord3p}=$self->{_end};
  }elsif($self->{_strand} eq "-"){
    $self->{_coord5p}=$self->{_end};
    $self->{_coord3p}=$self->{_start};
  }else{
    die "Error CodingExon.pm In subroutine new(), the strand for the CodingExon is wrong\n";
  }

  $self->{_sequenceAvailable} = 0;
  $self->{_proteinSequenceAvailable} = 0;
  $self->{_genomicSequence}   = "";

  $self->{_basePairsAtTheStart}   = 0;
  $self->{_basePairsAtTheEnd}     = 0;

  return $self;
}


=item get5PrimeEnd

    Get the 5' end coordinate of the coding exon
    on the chromosome

=cut
sub get5PrimeEnd{
  my ($self)=@_;
  return $self->{_coord5p};
}

=item get3PrimeEnd

    Get the 3' end coordinate of the coding exon
    on the chromosome

=cut
sub get3PrimeEnd{
  my ($self)=@_;
  return $self->{_coord3p};
}

=item getStart

    Get the lowest coordinate on the chromosome

=cut
sub getStart{
  my ($self)=@_;
  return $self->{_start};
}

=item getEnd

    Get the highest coordinate on the chromosome

=cut
sub getEnd{
  my ($self)=@_;
  return $self->{_end};
}

=item startsWithSTARTCodon

    Returns 1 if the coding exon starts with a
    start codon, 0 otherwise

=cut
sub startsWithSTARTCodon{
  my ($self)=@_;
  my $codon=substr($self->{_genomicSequence},0,3);
  if (exists $codon2aa{$codon}) {
    if($codon2aa{$codon} eq "M"){
      return 1;
    }else{
      return 0;
    }
  } else {
    warn "Error CodingExon.pm startsWithSTARTCodon() Codon ".$codon." not found\n";
    return 0;
  }
}


=item startsWithSTOPCodon

    Returns 1 if the coding exon ends with a
    stop codon, 0 otherwise

=cut
sub endsWithSTOPCodon{
  my ($self)=@_;
  my $codon=substr($self->{_genomicSequence},-3);
  if (exists $codon2aa{$codon}) {
    if($codon2aa{$codon} eq $stopSymbol){
      return 1;
    }else{
      return 0;
    }
  } else {
    warn "Error CodingExon.pm ".Dumper($self). "endsWithSTOPCodon() Codon ".$codon." not found\n";
    return 0;
  }

}



=item getIndexExon

    Returns the index of the coding exon
    the transcript

=cut
sub getIndexExon{
  my ($self)=@_;
  return $self->{_indexExon};
}


=item isCoding

    Returns 1 if the coding exon codes for
    a protein, 0 otherwise

=cut
sub isCoding{
  my ($self)=@_;
  return $self->{_isCoding};
}

=item getFrame

    Returns the numerical of
    the frame

=cut
sub getFrame{
  my ($self)=@_;
  return $self->{_frame};
}

=item setGenomicSequence

    Called by CodingTranscript.pm to set
    the original sequence. The internal sequence
    can be modified using appendStartGenomicSequence()
    and appendEndGenomicSequence() to add the base pairs
    from the previous and next coding exons.
    Arguments:
      genomicSequence : The original genomic argument

=cut
sub setGenomicSequence{
  my ($self,$genomicSequence)=@_;
  $self->{_genomicSequence}=$genomicSequence;

  $self->{_firstCharOfOrigSeq}  = substr($genomicSequence,0,1);
  $self->{_first2CharOfOrigSeq} = substr($genomicSequence,0,2);
  $self->{_lastCharOfOrigSeq}   = substr($genomicSequence,-1);
  $self->{_last2CharOfOrigSeq}  = substr($genomicSequence,-2);

  $self->{_isReliable}=1;

  $self->{_sequenceAvailable}=1;
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

=item setAsUnreliable

  Subroutine called by CodingTranscript.pm
  if it deems that the current coding exon
  is unrealiable.

=cut
sub setAsUnreliable{
  my($self)=@_;
  $self->{_isReliable}=0;
}

=item getGenomicSequence

  Subroutine to return the genomic
  sequence for this coding exon.

=cut
sub getGenomicSequence{
  my ($self)=@_;
  if($self->{_sequenceAvailable}){
    return $self->{_genomicSequence};
  }else{
    die "Error: CodingExon.pm in subroutine getGenomicSequence The genomic sequence was not set\n";
  }
}

=item getGenomicSequenceLength

  Subroutine to return the length of
  the genomic sequence for this coding exon.

=cut
sub getGenomicSequenceLength{
  my ($self)=@_;
  if($self->{_sequenceAvailable}){
    return length($self->{_genomicSequence});
  }else{
    die "Error: CodingExon.pm in subroutine getGenomicSequenceLength The genomic sequence was not set\n";
  }
}

=item getProteomicSequence

  Subroutine to return the translated
  sequence for the genomic sequence
  for this coding exon.

=cut
sub getProteomicSequence{
  my ($self)=@_;
  if($self->{_proteinSequenceAvailable}){
    return $self->{_proteinSequence};
  }else{
    die "Error: CodingExon.pm in subroutine getProteomicSequence The proteomic sequence was not set\n";
  }
}


=item appendStartGenomicSequence

  Subroutine used by CodingExon.pm to
  add base pairs from the previous exon
  at the beginning of the genomic sequence.

=cut
sub appendStartGenomicSequence{
  my ($self,$sequenceToAppend)=@_;
  if($self->{_sequenceAvailable}){
    $self->{_genomicSequence}=$sequenceToAppend.$self->{_genomicSequence};
    $self->{_basePairsAtTheStart} += length($sequenceToAppend);
  }else{
    die "Error: CodingExon.pm in subroutine appendStartGenomicSequence() The genomic sequence was not set\n";
  }
}

=item appendEndGenomicSequence

  Subroutine used by CodingExon.pm to
  add base pairs from the next exon
  at the end of the genomic sequence.

=cut
sub appendEndGenomicSequence{
  my ($self,$sequenceToAppend)=@_;
  if($self->{_sequenceAvailable}){
    $self->{_genomicSequence}=$self->{_genomicSequence}.$sequenceToAppend;
    $self->{_basePairsAtTheEnd} += length($sequenceToAppend);
  }else{
    die "Error: CodingExon.pm in subroutine appendEndGenomicSequence The genomic sequence was not set\n";
  }
}


=item getBasePairsAtTheStart

  Subroutine to return the length of
  the base pairs that were appended using
  appendStartGenomicSequence().

=cut
sub getBasePairsAtTheStart{
  my ($self)=@_;
  return $self->{_basePairsAtTheStart};
}

=item getBasePairsAtTheEnd

  Subroutine to return the length of
  the base pairs that were appended using
  appendEndGenomicSequence().

=cut
sub getBasePairsAtTheEnd{
  my ($self)=@_;
  return $self->{_basePairsAtTheEnd};
}


=item getFirstCharOfOrigSeq

  Subroutine to return the first char
  of the genomic sequence before
  appendStartGenomicSequence() was called.

=cut
sub getFirstCharOfOrigSeq{
  my ($self)=@_;

  return $self->{_firstCharOfOrigSeq};
}

=item getFirst2CharOfOrigSeq

  Subroutine to return the first 2 chars
  of the genomic sequence before
  appendStartGenomicSequence() was called.

=cut
sub getFirst2CharOfOrigSeq{
  my ($self)=@_;
  return $self->{_first2CharOfOrigSeq};
}


=item getLastCharOfOrigSeq

  Subroutine to return the last char
  of the genomic sequence before
  appendEndGenomicSequence() was called.

=cut
sub getLastCharOfOrigSeq{
  my ($self)=@_;
  return $self->{_lastCharOfOrigSeq};
}

=item getLast2CharOfOrigSeq

  Subroutine to return the last 2 chars
  of the genomic sequence before
  appendEndGenomicSequence() was called.

=cut
sub getLast2CharOfOrigSeq{
  my ($self)=@_;
  return $self->{_last2CharOfOrigSeq};
}



=item isCoordContained

  Subroutine to determine if the genomic coordinate is
  within the coding exon or not.
  Arguments:
   coordinate : coordinate on the chromosome

=cut
sub isCoordContained{
  my ($self,$coordinate)=@_;
  if( ($coordinate   < $self->{_start}) ||
      ($self->{_end} < $coordinate) ){
    return 0;
  }
  return 1;
}


=item isCoordUpstream

  Subroutine to determine if the genomic coordinate is
  upstream of the the coding exon or not.
  Arguments:
   coordinate : coordinate on the chromosome

=cut
sub isCoordUpstream{
  my ($self,$coordinate)=@_;
  if($self->{_strand} eq "+"){
    if($coordinate < $self->{_start}){
      return 1;
    }
  }elsif($self->{_strand} eq "-"){
    if($self->{_end} < $coordinate ){
      return 1;
    }
  }else{
    die "Error CodingExon.pm In subroutine isCoordUpstream(), the strand for the CodingExon is wrong\n";
  }
  return 0;
}

=item isCoordDownstream

  Subroutine to determine if the genomic coordinate is
  downstream of the the coding exon or not.
  Arguments:
   coordinate : coordinate on the chromosome

=cut
sub isCoordDownstream{
  my ($self,$coordinate)=@_;
  if($self->{_strand} eq "+"){
    if($self->{_end} < $coordinate ){
      return 1;
    }
  }elsif($self->{_strand} eq "-"){
    if($coordinate < $self->{_start}){
      return 1;
    }
  }else{
    die "Error CodingExon.pm In subroutine isCoordDownstream(), the strand for the CodingExon is wrong\n";
  }
  return 0;
}



=item translate2Protein

  Subroutine to compute the translated sequence
  using the genomic sequence for this coding exon.
  appendStartGenomicSequence and appendEndGenomicSequence
  must have been called prior to calling this subroutine
  to have a genomic sequence that has a length of a
  multiple of 3.

=cut
sub translate2Protein{
  my ($self)=@_;
  $self->{_proteinSequence}="";

  for (my $indexSeq=0;$indexSeq<length($self->{_genomicSequence});$indexSeq+=3) {
    if (exists $codon2aa{ substr($self->{_genomicSequence},$indexSeq,3) } ) {
      $self->{_proteinSequence}.=$codon2aa{ substr($self->{_genomicSequence},$indexSeq,3) };
    } else {
      warn "Error CodingExon.pm in translate2Protein() Codon ".substr($self->{_genomicSequence},$indexSeq,3)." not found, skipping\n";
    }
  }

  $self->{_proteinSequenceAvailable}=1;
}

=item translate2ProteinCustom

  Subroutine to compute the translated sequence
  using the genomic sequence for this coding exon.
  It uses the genomic sequence with the bases pairs
  appended at the beginning/end by
  appendStartGenomicSequence and appendEndGenomicSequence.
  It requires the base pairs of the codon in the 
  previous exon.
  Argument:
    sequenceToAddBeginning : The 1 or 2 base pairs to add
                             for the first codon if they are
                             on the previous coding exon.
  Returns:
   protSeq : The protein sequence for the current coding exon
   remainingBP : The remaining base pairs that were not matched to a codon


=cut
sub translate2ProteinCustom{
  my ($self,$sequenceToAddBeginning)=@_;

  my $genSeq;
  if($self->{_basePairsAtTheEnd} != 0){
    $genSeq  =$sequenceToAddBeginning.substr($self->{_genomicSequence}, $self->{_basePairsAtTheStart},-1*$self->{_basePairsAtTheEnd});
  }else{
    $genSeq  =$sequenceToAddBeginning.substr($self->{_genomicSequence}, $self->{_basePairsAtTheStart});
  }

  my $protSeq="";

  for (my $indexSeq=0;($indexSeq+2)<length($genSeq);$indexSeq+=3) {
    if (exists $codon2aa{ substr($genSeq,$indexSeq,3) } ) {
      $protSeq.=$codon2aa{ substr($genSeq,$indexSeq,3) };
    } else {
      die "Error CodingExon.pm in translate2ProteinCustom() Codon ".substr($genSeq,$indexSeq,3)." not found\n";
    }
  }

  my $remainingBP="";
  if(length($genSeq)%3 != 0 ){
    $remainingBP=substr($genSeq,-1*(length($genSeq)%3));
  }

  return ($protSeq,$remainingBP);
}



=item annotateBasePair

  Subroutine used by CodingExon.pm to
  annotate a snp once it has been determined
  that the snp lies within the current coding
  exon.
  Arguments:
   coordinate : coordinate on the chromosome
   refDNA     : DNA base pair of the reference sequence
   readDNA    : DNA base pair of the read
   produceFullProtSeq : flag, 0 = do not return full protein sequence, 1 do return
   noRefBP    : if set to 1, do not check the reference base pair
  Returns
   codonRef       : codon of the reference sequence
   codonRead      : codon of the read
   aaRef          : amino acid produced by codon of the reference sequence
   aaRead         : amino acid produced by codon of the read with the snp
   coordInAAExon  : coordinate of the snp in the amino acid sequence (if noRefBP=1)
   seqRef         : full protein sequence produced by the refence sequence (if noRefBP=1)
   seqRead        : full protein sequence produced by the read with the snp  (if noRefBP=1)

=cut
sub annotateBasePair{
  my ($self,$coordinate,$refDNA,$readDNA,$produceFullProtSeq,$noRefBP)=@_;

  if( ! $self->{_sequenceAvailable} ){
    die "Error CodingExon.pm Cannot call annotateBasePair() without setting the genomic sequence\n";
  }

  if(! $self->isCoordContained($coordinate) ){
    die "Error CodingExon.pm In subroutine annotateBasePair(), the coordinate $coordinate is not contained in the exon\n";
  }

  #verify if the refDNA is the same as the genomicSequence
  my $coordinateInString;
  if($self->{_strand} eq "+"){
    $coordinateInString  = ($coordinate - $self->{_start}) + $self->{_basePairsAtTheStart};
  }elsif($self->{_strand} eq "-"){
    $coordinateInString  = ($self->{_end} - $coordinate)   + $self->{_basePairsAtTheStart};

    if($noRefBP){
      $refDNA=substr($self->{_genomicSequence},$coordinateInString,1);
    }

    if(exists $dna2RC{$refDNA}){
      $refDNA=$dna2RC{$refDNA};
    }else{
      die "Error CodingExon.pm In subroutine annotateBasePair(), the reference DNA is an invalid base-pair ".$refDNA."\n";
    }

    if(exists $dna2RC{$readDNA}){
      $readDNA=$dna2RC{$readDNA};
    }else{
      die "Error CodingExon.pm In subroutine annotateBasePair(), the read DNA is an invalid base-pair ".$readDNA."\n";
    }


  }else{
    die "Error CodingExon.pm In subroutine annotateBasePair(), the strand for the CodingExon is wrong\n";
  }

  if($noRefBP){
    #do not check ref bp
    $refDNA = substr($self->{_genomicSequence},$coordinateInString,1);
  }else{
    if(substr($self->{_genomicSequence},$coordinateInString,1) ne $refDNA){
      die "Error CodingExon.pm In subroutine annotateBasePair(), the reference DNA ".$refDNA." is not the same as the ref ".substr($self->{_genomicSequence},$coordinateInString,1)." \n";
    }
  }
  my $coordMod3=  $coordinateInString%3;
  my $codonRef  ;
  my $codonRead ;
  my $aaRef  ;
  my $aaRead ;


  if($coordMod3 == 0){
    $codonRef  = $refDNA .substr($self->{_genomicSequence},$coordinateInString+1,2);
    $codonRead = $readDNA.substr($self->{_genomicSequence},$coordinateInString+1,2);
  }elsif($coordMod3 == 1){
    $codonRef  = substr($self->{_genomicSequence},$coordinateInString-1,1).$refDNA .substr($self->{_genomicSequence},$coordinateInString+1,1);
    $codonRead = substr($self->{_genomicSequence},$coordinateInString-1,1).$readDNA.substr($self->{_genomicSequence},$coordinateInString+1,1);
  }elsif($coordMod3 == 2){
    $codonRef  = substr($self->{_genomicSequence},$coordinateInString-2,2).$refDNA;
    $codonRead = substr($self->{_genomicSequence},$coordinateInString-2,2).$readDNA;
  }


  if (exists $codon2aa{$codonRef}) {
    $aaRef=$codon2aa{$codonRef};
  } else {
    die "Error CodingExon.pm Codon ".$codonRef." not found\n";
  }

  if (exists $codon2aa{$codonRead}) {
    $aaRead=$codon2aa{$codonRead};
  } else {
    die "Error CodingExon.pm Codon ".$codonRead." not found\n";
  }

  if($produceFullProtSeq){
    my $seqRef  ="";
    my $seqRead ="";
    my $coordInAAExon=0;

    for(my $indexString=0;$indexString<length($self->{_genomicSequence});$indexString+=3){

      if(  ($indexString<=$coordinateInString) &&
	   ($coordinateInString-$indexString)<= 2){
	$seqRef  .= $aaRef;
	$seqRead .= $aaRead;
	$coordInAAExon=length($seqRead);
      }else{
	my $aaToADD;
	if (exists $codon2aa{ substr($self->{_genomicSequence},$indexString,3) } ) {
	  $aaToADD=$codon2aa{ substr($self->{_genomicSequence},$indexString,3) };
	} else {
	  die "Error CodingExon.pm Codon ".$codonRef." not found\n";
	}
	$seqRef  .= $aaToADD;
	$seqRead .= $aaToADD;
      }
    }

    return ($codonRef,$codonRead,$aaRef,$aaRead,$coordInAAExon,$seqRef,$seqRead);
  }else{
    return ($codonRef,$codonRead,$aaRef,$aaRead);
  }

}



=item annotateInsert

  Subroutine used by CodingExon.pm to
  annotate an insertion within the exon.
  Arguments:
   coordinate : coordinate on the chromosome
   dnaSeq     : DNA sequence of the insertion
  Returns
   coordInAAExon  : Coordinate the insert within the amino acid sequence
   seqRef         : full protein sequence produced by the refence sequence (if noRefBP=1)
   seqRead        : full protein sequence produced by the read with the insert  (if noRefBP=1)
   untranslatedBp : the remaining base pairs that were not used to form a full codon
                    and need to be appended to the codon from the next exon

=cut
sub annotateInsert{
  my ($self,$coordinate,$dnaSeq)=@_;

  if( ! $self->{_sequenceAvailable} ){
    die "Error CodingExon.pm Cannot call annotateInsert() without setting the genomic sequence\n";
  }

  if(! $self->isCoordContained($coordinate) ){
    die "Error CodingExon.pm In subroutine annotateInsert(), the coordinate $coordinate is not contained in the exon\n";
  }

  #verify if the refDNA is the same as the genomicSequence
  my $coordinateInString;
  if($self->{_strand} eq "+"){
    $coordinateInString  = ($coordinate - $self->{_start}) + $self->{_basePairsAtTheStart};
  }elsif($self->{_strand} eq "-"){
    $coordinateInString  = ($self->{_end} - $coordinate)   + $self->{_basePairsAtTheStart};
  }else{
    die "Error CodingExon.pm In subroutine annotate(), the strand for the CodingExon is wrong\n";
  }

  my $modifiedDNAseq; #sequence without the base pairs are the end
  if($self->{_basePairsAtTheEnd} != 0){
    $modifiedDNAseq  = substr($self->{_genomicSequence},0,-1*$self->{_basePairsAtTheEnd});
  }else{
     $modifiedDNAseq = substr($self->{_genomicSequence},0); #does nothing
  }

  my $coordInAAExon=floor($coordinateInString/3);

  $modifiedDNAseq=substr($modifiedDNAseq,0,$coordinateInString).$dnaSeq.substr($modifiedDNAseq,$coordinateInString);

  my $seqRef  ="";
  my $seqRead ="";

  for(my $indexString=0;$indexString<length($self->{_genomicSequence});$indexString+=3){
      my $aaToADD;
      if (exists $codon2aa{ substr($self->{_genomicSequence},$indexString,3) } ) {
	$aaToADD=$codon2aa{ substr($self->{_genomicSequence},$indexString,3) };
      } else {
	die "Error CodingExon.pm Codon ".substr($self->{_genomicSequence},$indexString,3)." not found\n";
      }
      $seqRef  .= $aaToADD;
  }

  for(my $indexString=0;($indexString+2)<length($modifiedDNAseq);$indexString+=3){
      my $aaToADD;
      if (exists $codon2aa{ substr($modifiedDNAseq,$indexString,3) } ) {
	$aaToADD=$codon2aa{ substr($modifiedDNAseq,$indexString,3) };
      } else {
	die "Error CodingExon.pm Codon ".substr($modifiedDNAseq,$indexString,3)." not found\n";
      }
      $seqRead  .= $aaToADD;
  }



  my $untranslatedBp="";
  my $untranslatedBpLength=length($modifiedDNAseq)%3;

  if($untranslatedBpLength != 0 ){
    $untranslatedBp=substr($modifiedDNAseq, -1*( $untranslatedBpLength ));
  }

  return ($coordInAAExon,$seqRef,$seqRead,$untranslatedBp);

}


=item annotateInsert

  Subroutine used by CodingExon.pm to
  annotate a deletion within the exon.
  Arguments:
   coordinate1 : lowest coordinate of the deletion
   coordinate2 : highest coordinate of the deletion
   bpToAppend  : base pairs to append to the original genomic sequence to
                 complete the initial codon.
  Returns
   seqRead        : The protein sequence with the deletion
   untranslatedBp : the remaining base pairs that were not used to form a full codon
                    and need to be appended to the codon from the next exon

=cut

sub annotateDeletion{
  my ($self,$coordinate1,$coordinate2,$bpToAppend)=@_;

  if( ! $self->{_sequenceAvailable} ){
    die "Error CodingExon.pm ".Dumper($self)." Cannot call annotateDeletion() without setting the genomic sequence\n";
  }

  if( ($coordinate1   < $self->{_start}) ||
      ($self->{_end} < $coordinate2) ){
    die "Error CodingExon.pm In subroutine annotateDeletion(), the coordinate $coordinate1 is not contained in the exon\n";
  }
  if( ($coordinate2   < $self->{_start}) ||
      ($self->{_end} < $coordinate2) ){
    die "Error CodingExon.pm In subroutine annotateDeletion(), the coordinate $coordinate2 is not contained in the exon\n";
  }

  my $coordinateInString1;
  if($self->{_strand} eq "+"){
    $coordinateInString1  = ($coordinate1 - $self->{_start}) ;#+ $self->{_basePairsAtTheStart};
  }elsif($self->{_strand} eq "-"){
    $coordinateInString1  = ($self->{_end} - $coordinate1)   ;#+ $self->{_basePairsAtTheStart};
  }else{
    die "Error CodingExon.pm In subroutine annotateDeletion(), the strand for the CodingExon is wrong\n";
  }

  my $coordinateInString2;
  if($self->{_strand} eq "+"){
    $coordinateInString2  = ($coordinate2 - $self->{_start}) ;#+ $self->{_basePairsAtTheStart};
  }elsif($self->{_strand} eq "-"){
    $coordinateInString2  = ($self->{_end} - $coordinate2)   ;#+ $self->{_basePairsAtTheStart};
  }else{
    die "Error CodingExon.pm In subroutine annotateDeletion(), the strand for the CodingExon is wrong\n";
  }

  if($coordinateInString1 > $coordinateInString2){
    my $temp=$coordinateInString2;
    $coordinateInString2=$coordinateInString1;
    $coordinateInString1=$temp;
  }

  my $seqRead ="";

  my $modifiedDNAseq; #sequence without the base pairs are the end
  if($self->{_basePairsAtTheStart} != 0){
    $modifiedDNAseq  = substr($self->{_genomicSequence},$self->{_basePairsAtTheStart});
  }else{
    $modifiedDNAseq = $self->{_genomicSequence};
  }

  if($self->{_basePairsAtTheEnd} != 0){
    $modifiedDNAseq  = substr($modifiedDNAseq,0,-1*$self->{_basePairsAtTheEnd});
  }else{
    $modifiedDNAseq = substr($modifiedDNAseq,0); #does nothing
  }

  $modifiedDNAseq = $bpToAppend. substr($modifiedDNAseq,0,$coordinateInString1) .substr($modifiedDNAseq,$coordinateInString2+1);

  for(my $indexString=0;($indexString+2)<length($modifiedDNAseq);$indexString+=3){
      my $aaToADD;
      if (exists $codon2aa{ substr($modifiedDNAseq,$indexString,3) } ) {
	$aaToADD=$codon2aa{ substr($modifiedDNAseq,$indexString,3) };
      } else {
	die "Error CodingExon.pm Codon ".substr($modifiedDNAseq,$indexString,3)." not found\n";
      }
      $seqRead  .= $aaToADD;
  }

  my $untranslatedBp="";
  my $untranslatedBpLength=length($modifiedDNAseq)%3;

  if($untranslatedBpLength != 0 ){
    $untranslatedBp=substr($modifiedDNAseq, -1*( $untranslatedBpLength ));
  }

  return ($seqRead,$untranslatedBp);

}





1;

