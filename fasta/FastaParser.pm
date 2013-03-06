=head1 NAME

   FastaParser.pm


=head1 Disclaimer

This software and all the data therein were developed by the Bioinformatics and Computational Biology Lab (Laboratorio de Bioinformatica e Biologia Computational (LBBC)) of the Brazilian National Cancer Institute (Instituto Nacional de Cancer (INCA)) in Rio de Janeiro, RJ, Brazil. The INCA is part of the Health Ministry (Ministerio da Saude) of the Federal Government of Brazil. This software was developed by the staff of the LBBC as part of their official duties. According to federal regulations, this software is therefore public domain.

This software is provided as is and the INCA does make any guarantee regarding the accuracy/quality/reliability of the results. The INCA cannot be held accountable for any hardware or software damage that might be caused through the use or misuse of this software.

The LBBC encourages researchers using our software to acknowledge us. The software can be distributed freely and modified given that any derivative bears an acknowledgment that it was derived from our original software.

=head1 Author

Gabriel Renaud grenaud@inca.gov.br



=cut


package FastaParser;

use strict;
use warnings;

use Data::Dumper;

use FastaRecord;

my $FASTAFILEPTR;



=item new

  Constructor, initializes the file pointer and
  calls parseFastaSequence() to get the next record

=cut

sub new{
  my ($class,%arg)=(@_);

  if(!$arg{filename}){
    die "Please enter the filename\n";
  }

  if( !(-e $arg{filename}) ){
    die "FastaParser.pm File ".$arg{filename}." does not exist\n";
  }

  my ($self) =bless{_filename  => $arg{filename}
		   }, $class;

  if( !(-e $self->{_filename}) ){
    die "FastaParser.pm File ".$self->{_filename}." does not exist\n";
  }

  $self->{_nextNextToReturn}=undef;
  $self->{_nextToReturn}=undef;
  $self->{_recordToReturn}=0;
  $self->{_recordLeft}=1;


  open($FASTAFILEPTR, $self->{_filename}) or die ("Can't open ".$self->{_filename}."\n");
  detectFileFormat($self );#, $self->{_filename} );
  open($FASTAFILEPTR, $self->{_filename}) or die ("Can't open ".$self->{_filename}."\n");
  parseFastaSequence($self,$self->{_filename});

  if(!$self->{_recordToReturn}){
    warn "FastaParser.pm  File ".$self->{_filename}." does not contain any record\n";
  }

  return $self;
}






=item hasRecord

   Subroutine to verify if there is a record to
   return using getNextRecord()

=cut
sub hasRecord {
  my($self)=@_;
  return $self->{_recordToReturn};
}





=item getNextRecord

  Subroutine to return the next record and parse the file to find
  the next record.

=cut
sub getNextRecord {
  my($self)=@_;
  $self->{_nextToReturn}=  $self->{_nextNextToReturn};
  if($self->{_recordLeft}){
    parseFastaSequence($self);
  }else{
    $self->{_recordToReturn}=0;
  }

  return $self->{_nextToReturn};
}




=item parseFastaSequence

  Method to go through the file (given that myfastafile has been initialized)
  Creates a FastaRecord and stores it in nextNextToReturn

=cut
sub parseFastaSequence {
  my($self)=@_;
  my $line;

  my $defline;
  my $sequence;
  my $filename=$self->{_filename};

  my $waitingForDefline   =1;
  my $firstLineOfSequence =0;
  my $inSequence          =0;
  my $nextLineShouldBeLast=0;
  my $foundOne=0;


  my $lengthOfLine;

  my $keepLooping=1;




  while ($keepLooping) {
    my $testEOF=($line = <$FASTAFILEPTR>);

    #If we have reached the bottom of the file
    if (!$testEOF) {
      $keepLooping=0;
      $self->{_recordLeft}=0;
      close($FASTAFILEPTR);

      #If we found a record, store it
      if($foundOne){
	if(!$sequence){
	  die "FastaParser.pm  No sequence was found in ".$filename." for sequence ".$defline."\n";
	}
	$self->{_nextNextToReturn}=FastaRecord->new(defline  => $defline,
						    sequence => $sequence);
      }

      #If we were expecting the first line and got to the end, this means the sequence was empty
      if($firstLineOfSequence){
	die "FastaParser.pm  No sequence was found in ".$filename." for sequence ".$defline."\n";
      }

    } else {
      #Trim \n at the end
      chomp($line);
      #print "line = #".$line."# $waitingForDefline $firstLineOfSequence $inSequence $foundOne\n";
      if(substr($line,0,1) eq "#"){
	next;
      }

      #If we found the defline and the first sequence line
      if ($inSequence) {
	my $found=index($line,">");

	#if we found a >, the start of a new record
	if ($found != -1) {
	  #The > was not at index 0, generate error
	  if ($found !=0) {
	    die "FastaParser.pm Error while parsing ".$filename." '>' found in unexpected index\n";
	  } else {
	    #Store the current record and re-position the pointer for next time around
	    $self->{_nextNextToReturn}=FastaRecord->new(defline  => $defline,
							sequence => $sequence);
	    $keepLooping=0;

	    seek($FASTAFILEPTR,tell($FASTAFILEPTR)-(length($line)+1),0);
	  }
	} else {
	  my $lineToAdd=trimAndCheckLine($line);

	  #The line is empty, we reached the end of the current record
	  if (length($lineToAdd) == 0) {
	    $self->{_nextNextToReturn}=FastaRecord->new(defline  => $defline,
							sequence => $sequence);
	    $keepLooping=0;

	  } else {#normal case, we found another sequence line
	    if($nextLineShouldBeLast){
	      die "FastaParser: Potential error in ".$filename." wrong length with line ".$line."\n";
	    }

	    $sequence.=$lineToAdd;
	    #If the length of the current line is greater than what we found so far
	    if ($lengthOfLine < length($lineToAdd)) {
	      die "FastaParser: Potential error in ".$filename." wrong length with line ".$line."\n";
	    } else {
	      #If the length of the line is lesser that we we found so far, it becomes the new maximum
	      if ( $lengthOfLine > length($lineToAdd) ) {
		#$lengthOfLine = length($lineToAdd);
		$nextLineShouldBeLast=1;
	      }
	    }
	  }
	}

      }


      #First line of sequences
      if ($firstLineOfSequence) {
	my $found=index($line,">");
	#if we found another >, the previous sequence was empty
	if ($found != -1) {
	  die "FastaParser.pm  No sequence was found in ".$filename." for sequence ".$defline."\n";
	} else {
	  my $lineToAdd=trimAndCheckLine($line);
	  #If the line is empty, throw error
	  if (length($lineToAdd) == 0) {
	    die "FastaParser.pm No sequence was found in ".$filename." for sequence ".$defline."\n";
	  } else {
	    #If the line is not empty, measure
	    $sequence=$lineToAdd;
	    $lengthOfLine=length($sequence);
	    $firstLineOfSequence=0;
	    $inSequence=1;
	  }
	}
      }


      #Waiting to find a >
      if ($waitingForDefline) {
	my $found=index($line,">");
	#Found one
	if ($found != -1) {
	  #The > was not found at the index 0
	  if ($found !=0) {
	    die "FastaParser.pm error while parsing file ".$filename." '>' found at unexpected index\n";
	  } else {
	    $defline=trimLine($line);
	    $waitingForDefline   =0;
	    $firstLineOfSequence =1;
	    $foundOne            =1;
	    $nextLineShouldBeLast=0;
	  }
	} else {
	  #If the line does not have a > and is not blank
	  if (length(trimLine($line)) != 0) {
	    die "FastaParser.pm Error while parsing ".$filename." non-white char while expecting '>'\n";
	  }
	}

      }
    }

  }

  $self->{_recordToReturn}=$foundOne;
}













=item detectFileFormat

  This subroutine makes sure that the file is in UNIX format
  (NO CR, just NL)

=cut
sub detectFileFormat  {
  my($self,$filename)=@_;
  my $line = <$FASTAFILEPTR>;
  #die Dumper($self);
  #print "filename $filename #\n";
  if( ord(substr($line,length($line)-2,1)) == 13){
    die "Please make sure that the file ".$self->{_filename}." is in UNIX format\n";
  }

}









=item trimAndCheckLine

  This subroutine calls trimLine on lineToCheck
  and makes sure that no white spaces are present
  in lineToCheck

=cut
sub trimAndCheckLine  {
  my($lineToCheck)=@_;

  $lineToCheck =trimLine($lineToCheck);
  if( $lineToCheck =~ /\s/ ){
    die "The line $lineToCheck has white spaces\n";
  }

  return $lineToCheck;
}



=item trimLine

  This subroutine trims the the white spaces/tabs
  that may exist at the beginning/end of a string
  and returns it.

=cut
sub trimLine  {
  my($lineToCheck)=@_;

  $lineToCheck =~ s/^\s+//;
  $lineToCheck =~ s/\s+$//;

  return $lineToCheck;
}






1;





