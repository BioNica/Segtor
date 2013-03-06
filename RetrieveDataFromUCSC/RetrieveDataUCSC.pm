=head1 NAME

   RetrieveDataUCSC::RetrieveDataUCSC.pm


=head1 Disclaimer

This software and all the data therein were developed by the Bioinformatics and Computational Biology Lab (Laboratorio de Bioinformatica e Biologia Computational (LBBC)) of the Brazilian National Cancer Institute (Instituto Nacional de Cancer (INCA)) in Rio de Janeiro, RJ, Brazil. The INCA is part of the Health Ministry (Ministerio da Saude) of the Federal Government of Brazil. This software was developed by the staff of the LBBC as part of their official duties. According to federal regulations, this software is therefore public domain.

This software is provided as is and the INCA does make any guarantee regarding the accuracy/quality/reliability of the results. The INCA cannot be held accountable for any hardware or software damage that might be caused through the use or misuse of this software.

The LBBC encourages researchers using our software to acknowledge us. The software can be distributed freely and modified given that any derivative bears an acknowledgment that it was derived from our original software.

=head1 Author

Gabriel Renaud grenaud@inca.gov.br



=cut



package RetrieveDataUCSC;


use strict;
use warnings;

use Data::Dumper;
use Time::Local;
use LWP;




local $| = 1;


my $finalSize;
my $FILEOUT;
my $contentSize;
my $content;

my $sleepTimeSeconds=30;



=item printBar

   Returns a string with the progress with
   $current done out of $target

=cut
sub printBar{
  my ($current,$target)=@_;
  my $numberOfEquals=40;

  my $stringToReturn="";

  $stringToReturn.="|";
  my $numberOfEqualsToPrint;
  if($target == 0 ){
    $numberOfEqualsToPrint=$numberOfEquals;
  }else{
    $numberOfEqualsToPrint=int(($current/$target)*$numberOfEquals);
  }

  for(my $equals=0;$equals<$numberOfEqualsToPrint;$equals++){
    $stringToReturn.="=";
  }
  for(my $equals=$numberOfEqualsToPrint;$equals<$numberOfEquals;$equals++){
    $stringToReturn.=" ";
  }
  $stringToReturn.="|";
  if($target == 0 ){
    $stringToReturn.="  ".int(100)."%";
  }else{
    $stringToReturn.="  ".int(($current/$target)*100)."%";
  }
  return $stringToReturn;
}


=item contentCallBACK

   As custom callback function, puts the size of the content in
   $contentSize and the content in $content. Called by retrieveURLNoFinalSize()
   Does not display the progress bar since it does not have information about the final size
   Needs to have $contentSize to zero and $content initialized prior to dowloading

=cut
sub contentCallBACK {
  my ($dataChunk, $resp, $prot) = @_;
  $contentSize+= length $dataChunk;
  $content.=$dataChunk;
}



=item retrieveURLNoFinalSize

    Retrieve a URL without information as to the final size.
    Must initialize $content and $contentSize

=cut

sub retrieveURLNoFinalSize {
  my ($self,$url)=@_;
  my $req;
  my $response;

  $content="";
  $contentSize=0;

  $req      = new HTTP::Request 'GET' => $url;
  $response = $self->{_ua}->request($req,\&contentCallBACK);
  $req->authorization_basic('anonymous', 'anonymous');

  return $response;
}






=item testInternet

    Returns 0 if the module fails to retrieve $self->{_ftpTarget}, 1 otherwise


=cut

sub testInternet {
  my ($self)=@_;

  my $response = retrieveURLNoFinalSize($self,$self->{_ftpTarget});

  if(! ($response->is_success) ){
    return 0;
  }

  if (!$content) {
    return 0;
  }

  return 1;
}





=item contentCallBACKWrite

   As custom callback function, puts the size of the content in
   $contentSize writes the content in FILEOUT. Called by retrieveURLFinalSize()
   Needs to have FILEOUT and $finalSize to be set and $contentSize to zero

=cut

sub contentCallBACKWrite {
  my ($dataChunk, $resp, $prot) = @_;
  $contentSize+= length $dataChunk;
  print FILEOUT $dataChunk;

  #display progress bar
  my $stringToPrint=printBar($contentSize,$finalSize);
  print "\r";
  print "\t";
  print $stringToPrint;
}

=item retrieveURLFinalSize

    Retrieve a URL without information as to the final size.
    Must initialize $content and $contentSize

=cut

sub retrieveURLFinalSize {
  my ($self,$url,$finalSize2,$fileChrToWriteTo)=@_;
  my $req;
  my $response;

  #sets the final size
  $finalSize=$finalSize2;
  $contentSize=0;

  #Initialized the file pointer
  open(FILEOUT,">".$fileChrToWriteTo) or die "Cannot open $fileChrToWriteTo\n";

  $req      = new HTTP::Request 'GET' => $url;
  $response = $self->{_ua}->request($req,\&contentCallBACKWrite);
  $req->authorization_basic('anonymous', 'anonymous');
  #close the file pointer
  close(FILEOUT);
  print "\n";
  return $response;
}




=item new

    Constructor for the object that requires the following attributes:
    'ftpTarget'        ftp address of UCSC site (usually ftp://hgdownload.cse.ucsc.edu)
    'httpProxy'        ip address of the http proxy if any
    'ftpProxy'         ip address of the ftp proxy if any
    'waisProxy'        ip address of the wais proxy if any

    It creates an instante of a LWP::UserAgent for ftp retrieval

=cut

sub new{
  my ($class,%arg)=(@_);
  my ($self) =bless{_ftpTarget  => $arg{ftpTarget},
                    _httpProxy  => $arg{httpProxy},
		    _ftpProxy   => $arg{ftpProxy},
		    _waisProxy  => $arg{waisProxy}}, $class;



  $self->{_ua}=new LWP::UserAgent;
  $self->{_ua}->agent("Mozilla/3.0");

  if($self->{_httpProxy}){
    $self->{_ua}->proxy(http => $self->{_httpProxy});
  }

  if($self->{_ftpProxy}){
    $self->{_ua}->proxy(ftp  => $self->{_ftpProxy});
  }

  if($self->{_waisProxy}){
    $self->{_ua}->proxy(wais => $self->{_waisProxy});
  }

  #detect path of gunzip/unzip
  if(-e "/usr/bin/gunzip"){
    $self->{_gunzip}="/usr/bin/gunzip";
  }elsif(-e "/bin/gunzip"){
    $self->{_gunzip}="/bin/gunzip";
  }else{
    die "Cannot detect gunzip on this system\n";
  }

  if(-e "/usr/bin/unzip"){
    $self->{_unzip}="/usr/bin/unzip";
  }elsif(-e "/bin/unzip"){
    $self->{_unzip}="/bin/unzip";
  }else{
    die "Cannot detect unzip on this system\n";
  }


  return $self;
}


=item getSizeDatabase

    Subroutine to retrieve the size of the database file
    named $databaseName based on the database listing in
    $url.

    Returns ($log,$warnLog,$dieLog,$finalSizeFound);

=cut

sub getSizeDatabase{
  my ($self,$speciesCode,$databaseName)=@_;
  my $log="";
  my $warnLog="";
  my $dieLog;
  my $url="ftp://hgdownload.cse.ucsc.edu/goldenPath/".$speciesCode."/database/";

  $log.="Retrieving directory listing from $url\n";
  print "Retrieving directory listing from $url\n";
  my $response;

  #Calls retrieveURLNoFinalSize and stores the content in $content
  while (1) {
    $response=retrieveURLNoFinalSize($self,$url);

    if ($response->is_success) {
    } else {
      $warnLog.="Unable to retrieve $url\n"."Status line ".$response->status_line."\nwaiting $sleepTimeSeconds seconds to retry again\n";
      warn "Unable to retrieve $url\n"."Status line ".$response->status_line."\nwaiting $sleepTimeSeconds seconds to retry again\n";
      sleep($sleepTimeSeconds);
    }

    if (!$content) {
      $log.="Content is empty\n"."Status line ".$response->status_line."\nwaiting $sleepTimeSeconds seconds to retry again\n";
      sleep($sleepTimeSeconds);
    } else {
      last;
    }
  }



  my $waitingForParent=1;
  my $inChr           =0;
  my $doneSeeingChr   =0;

  # reads the directory listing in $content until the
  # $databaseName is found
  foreach my $line (split("\n",$content)) {
    if ($doneSeeingChr) {
      if ($line =~ /<A HREF/i) {
	$dieLog="Problem#2 in parsing $url\n";
	return ($log,$warnLog,$dieLog,-1);
      }
    }

    if ($inChr) {
      #                                    1                                        2
      #-r--r--r--   1 ftp      ftp         13635 Mar 20  2009 <a href="/goldenPath/hg19/chromosomes/chrUn_gl000240.fa.gz">chrUn_gl000240.fa.gz</a>
      #-r--r--r--   1 ftp      ftp         13635 Mar 20  2009 chrUn_gl000240.fa.gz
      #                   1                                      2
      if( ($line =~ /\s+(\d+)\s+\S+\s+\d+\s+[\d\:]+\s+<A HREF=\"(\S+)\">/i) ||
	  ($line =~ /\s+(\d+)\s+\S+\s+\d+\s+[\d\:]+\s+(\S+)/) ){

	my $fileChrName=$2;
	my $fileChrToWriteTo;
	my $fileChrNameFasta;
	my $finalSizeFound=$1;

	if( ($fileChrName =~ /^[\w\/]*?\/?(\w+)\.txt\.gz$/  ) ||
	    ($fileChrName =~ /^[\w\/]*?\/?(\w+)\.txt\.zip$/  ) ) {
	  $fileChrToWriteTo = $1;

	  #if the $databaseName is found, return it
	  if($fileChrToWriteTo eq $databaseName){
	    return ($log,$warnLog,$dieLog,$finalSizeFound);
	  }
	} else {
	  next;
	}

      } else {
	if ($line =~ /<\/PRE>/i) {
	  $inChr=0;
	  $doneSeeingChr=1;
	} else {
	  $dieLog="Problem#1 in parsing $url with line $line\n";
	  return ($log,$warnLog,$dieLog,-1);
	}
      }
    }

    if ($waitingForParent) {
      #if ($line =~ /<A HREF=\"\S+\">Up to higher level directory<\/A>/i) {
      if ($line =~ /\.\./) {
	$waitingForParent=0;
	$inChr=1;
      }
    }
  }


  return ($log,$warnLog,$dieLog,-1);
}












=item listChromosomes

    Subroutine to retrieve the list of chromosomes
    given the $speciesCode.

    Returns an array of hashes where each hash
    has the following the properties:
    'size'     The size of the file
    'chr'      The name of the chromosome file
    'chrFull'  The name of the chromosome file

=cut

sub listChromosomes{
  my ($self,$speciesCode)=@_;
  my @arrayOfChromosomes;
  my $log="";
  my $warnLog="";
  my $dieLog;



  my $url=$self->{_ftpTarget}."/goldenPath/".$speciesCode."/chromosomes/";
  my $chrPath="/goldenPath/".$speciesCode."/chromosomes/";
  $log.="Retrieving directory listing from $url\n";

  my $response;

  while (1) {

    $response=retrieveURLNoFinalSize($self,$url);

    if ($response->is_success) {
    } else {
      $warnLog.="Unable to retrieve $url\n"."Status line ".$response->status_line."\nwaiting $sleepTimeSeconds seconds to retry again\n";
      warn "Unable to retrieve $url\n"."Status line ".$response->status_line."\nwaiting $sleepTimeSeconds seconds to retry again\n";
      sleep($sleepTimeSeconds);
    }

    if (!$content) {
      $log.="Content is empty\n"."Status line ".$response->status_line."\nwaiting $sleepTimeSeconds seconds to retry again\n";
      sleep($sleepTimeSeconds);
    } else {
      last;
    }
  }


  my $waitingForParent=1;
  my $inChr           =0;
  my $doneSeeingChr   =0;


  foreach my $line (split("\n",$content)) {
    if ($doneSeeingChr) {
      if ($line =~ /<A HREF/i) {
	$dieLog="Problem#2 in parsing $url\n";
	return ($log,$warnLog,$dieLog,@arrayOfChromosomes);
      }
    }

    if ($inChr) {

      #                   1                                      2
      #if ($line =~ /\s+(\d+)\s+\S+\s+\d+\s+[\d\:]+\s+<A HREF=\"(\S+)\">/i) {
      if(  ($line =~ /\s+(\d+)\s+\S+\s+\d+\s+[\d\:]+\s+<A HREF=\"(\S+)\">/i) ||
	   ($line =~ /\s+(\d+)\s+\S+\s+\d+\s+[\d\:]+\s+(\S+)/) ){
	my $fileChrName=$2;
	my $fileChrToWriteTo;
	my $fileChrNameFasta;
	my $finalSizeFound=$1;
	#print "fileChrName $fileChrName\n";
	if( ($fileChrName =~ /^[\w\/]*?\/?([\w\.]+\.fa\.gz)$/) ||  ($fileChrName =~ /^[\w\/]*?\/?([\w\.]+\.fa\.zip)$/) ) {
	  $fileChrToWriteTo = $1;
	} else {
	  $log.="Skipping $fileChrName\n" ;
	  next;
	}

	# if it is a chromosome, create a hash and store it in the array
	if($fileChrToWriteTo eq $fileChrName){
	  $fileChrName=$chrPath.$fileChrName;
	}

	my $hashChr={'size'     => $finalSizeFound,
		     'chr'      => $fileChrToWriteTo,
		     'chrFull'  => $fileChrName};

	push(@arrayOfChromosomes,$hashChr);
      } else {
	if ($line =~ /<\/PRE>/i) {
	  $inChr=0;
	  $doneSeeingChr=1;
	} else {
	  $dieLog="Problem#1 in parsing $url with line $line\n";
	  return ($log,$warnLog,$dieLog,@arrayOfChromosomes);
	}
      }
    }

    if ($waitingForParent) {
      #if ($line =~ /<A HREF=\"\S+\">Up to higher level directory<\/A>/i) {
      if ($line =~ /\.\./) {
	$waitingForParent=0;
	$inChr=1;
      }
    }
  }


  return ($log,$warnLog,$dieLog,@arrayOfChromosomes);
}



















=item listDatabases

    Subroutine to retrieve the list of databases
    given the $speciesCode.

    Returns ($log,$warnLog,$dieLog,$finalSizeFound);

=cut

sub listDatabases{
  my ($self,$speciesCode)=@_;
  my $log="";
  my $warnLog="";
  my $dieLog;
  my $url="ftp://hgdownload.cse.ucsc.edu/goldenPath/".$speciesCode."/database/";
  my $dataPath="/goldenPath/".$speciesCode."/database/";
  my @arrayOfDatabases;


  $log.="Retrieving directory listing from $url\n";
  print "Retrieving directory listing from $url\n";
  my $response;

  #Calls retrieveURLNoFinalSize and stores the content in $content
  while (1) {
    $response=retrieveURLNoFinalSize($self,$url);

    if ($response->is_success) {
    } else {
      $warnLog.="Unable to retrieve $url\n"."Status line ".$response->status_line."\nwaiting $sleepTimeSeconds seconds to retry again\n";
      warn "Unable to retrieve $url\n"."Status line ".$response->status_line."\nwaiting $sleepTimeSeconds seconds to retry again\n";
      sleep($sleepTimeSeconds);
    }

    if (!$content) {
      $log.="Content is empty\n"."Status line ".$response->status_line."\nwaiting $sleepTimeSeconds seconds to retry again\n";
      sleep($sleepTimeSeconds);
    } else {
      last;
    }
  }



  my $waitingForParent=1;
  my $inChr           =0;
  my $doneSeeingChr   =0;

  # reads the directory listing in $content until the
  # $databaseName is found
  foreach my $line (split("\n",$content)) {
    if ($doneSeeingChr) {
      if ($line =~ /<A HREF/i) {
	$dieLog="Problem#2 in parsing $url\n";
	return ($log,$warnLog,$dieLog,-1);
      }
    }

    if ($inChr) {

      #                   1                                      2
      #if ($line =~ /\s+(\d+)\s+\S+\s+\d+\s+[\d\:]+\s+<A HREF=\"(\S+)\">/i) {
      if(  ($line =~ /\s+(\d+)\s+\S+\s+\d+\s+[\d\:]+\s+<A HREF=\"(\S+)\">/i) ||
	   ($line =~ /\s+(\d+)\s+\S+\s+\d+\s+[\d\:]+\s+(\S+)/) ){
	my $fileChrName=$2;
	my $fileChrToWriteTo;
	my $fileChrNameFasta;
	my $finalSizeFound=$1;

	if( ($fileChrName =~ /^[\w\/]*?\/?(\w+)\.txt\.gz$/) ||  ($fileChrName =~ /^[\w\/]*?\/?(\w+)\.txt\.zip$/) ) {
	  $fileChrToWriteTo = $1;

	  if( index($fileChrName,$fileChrToWriteTo) == 0 ){
	    $fileChrName=$dataPath.$fileChrName;
	  }

	  my $hashChr={'size'      => $finalSizeFound,
		       'name'      => $fileChrToWriteTo,
		       'nameFull'  => $fileChrName};
	  push(@arrayOfDatabases,$hashChr);
	} else {
	  next;
	}

      } else {
	if ($line =~ /<\/PRE>/i) {
	  $inChr=0;
	  $doneSeeingChr=1;
	} else {
	  $dieLog="Problem#1 in parsing $url with line $line\n";
	  return ($log,$warnLog,$dieLog,@arrayOfDatabases);
	}
      }
    }

    if ($waitingForParent) {
      #if ($line =~ /<A HREF=\"\S+\">Up to higher level directory<\/A>/i) {
      if ($line =~ /\.\./) {
	$waitingForParent=0;
	$inChr=1;
      }
    }
  }

  return ($log,$warnLog,$dieLog,@arrayOfDatabases);
}






=item downloadGenome

  Subroutine to retrieve the entire genome
  for a $speciesCode. It retrieves the name of the
  chromosomes and their size using listChromosomes()
  It calls retrieveURLFinalSize() on each chromosome
  It returns ($log,$warnLog,$dieLog,@chromosomeNames)

=cut

sub downloadGenome{
  my ($self,$speciesCode)=@_;


  my $log="";
  my $warnLog="";
  my $dieLog;

  my @chromosomeNames;
  my ($logTemp,$warnLogTemp,$dieLogTemp,@arrayOfChromosomes)=listChromosomes($self,$speciesCode);
  $log.=$logTemp;
  $warnLog.=$warnLogTemp;
  if($dieLogTemp){
    $dieLog = $dieLogTemp;
    return ($log,$warnLog,$dieLog,@chromosomeNames);
  }

  my $indexFile=0;
  foreach my $hashRecord (@arrayOfChromosomes) {
    $indexFile++;
    my $fileChrToWriteTo = $hashRecord->{'chr'};
    my $fileChrName      = $hashRecord->{'chrFull'};
    my $finalSizeFound   = $hashRecord->{'size'};
    my $fileChrNameFasta;

    my $url2=$self->{_ftpTarget}.$fileChrName;
    $log.="Retrieving chromosomes $fileChrName ".$url2."\n" ;

    while (1) {
      print "RetrieveDataUCSC: Retrieving chromosome $url2 file ".$indexFile." of ".($#arrayOfChromosomes+1)."\n";

      #retrieves the content of $urls using retrieveURLFinalSize
      my $response2=retrieveURLFinalSize($self,$url2,$finalSizeFound,$fileChrToWriteTo);

      if ($response2->is_success) {
      } else {
	$warnLog.="Unable to retrieve $url2\n"."Status line ".$response2->status_line."\nwaiting $sleepTimeSeconds seconds to retry again\n";
	warn "Unable to retrieve $url2\n"."Status line ".$response2->status_line."\nwaiting $sleepTimeSeconds seconds to retry again\n";
	sleep($sleepTimeSeconds);
      }


      if ($contentSize == 0) {
	$log.="Content is empty\n"."Status line ".$response2->status_line."\nwaiting $sleepTimeSeconds seconds to retry again\n";
	sleep($sleepTimeSeconds);
      } else {
	my $extention;
	if( ($fileChrName =~ /^.*\/([^\/]+)\.(gz)$/) || ($fileChrName =~ /^.*\/([^\/]+)\.(zip)$/) ) {
	  $fileChrNameFasta=$1;
	  $extention=$2;
	} else {
	  $dieLog="Wrong file patern\n";
	  return ($log,$warnLog,$dieLog,@chromosomeNames);
	}
	my $command;
	if($extention eq "gz"){
	  #$command="/usr/bin/gunzip --force $fileChrToWriteTo";
	  $command=$self->{_gunzip}." --force $fileChrToWriteTo";
	}else{
	  #$command="/usr/bin/unzip -o $fileChrToWriteTo && /bin/rm -f $fileChrToWriteTo";
	  $command=$self->{_unzip}." -o $fileChrToWriteTo && /bin/rm -f $fileChrToWriteTo";
	}

	my $output=system($command);
	if ($output != 0 ) {
	  $log.="Command $command failed restarting\n";
	} else {
	  push(@chromosomeNames,$fileChrNameFasta);
	  last;
	}
      }
    }
  }

  return ($log,$warnLog,$dieLog,@chromosomeNames);
}





=item downloadGenome

  Subroutine to retrieve a specific database $databaseName
  for a $speciesCode. It retrieves the size of the database
  using getSizeDatabase().
  It calls retrieveURLFinalSize() on each chromosome
  It returns ($log,$warnLog,$dieLog,@chromosomeNames)

=cut

sub downloadDatabase{
  my ($self,$speciesCode,$databaseName)=@_;


  my $log="";
  my $warnLog="";
  my $dieLog;


  if($databaseName =~ /(\w+)(\.txt)?(\.gz)?/){
    $databaseName=$1;
  }else{
    $dieLog="Please verify the database name\n" ;
    return ($log,$warnLog,$dieLog);
  }

  my $fileDatabaseName=$databaseName.".txt.gz";


  my $url2="ftp://hgdownload.cse.ucsc.edu/goldenPath/".$speciesCode."/database/";
  my ($logTemp,$warnLogTemp,$dieLogTemp,$databaseNameSize)=getSizeDatabase($self,$speciesCode,$databaseName);

  $log     .= $logTemp;
  $warnLog .= $warnLogTemp;
  if($dieLogTemp){
    $dieLog = $dieLogTemp;
    return ($log,$warnLog,$dieLog);
  }




  if($databaseNameSize == -1 ){
    $dieLog="Unable to get database size\n" ;
    return ($log,$warnLog,$dieLog);
  }

  $url2.=$fileDatabaseName;
  $log.="Retrieving database $databaseName ".$url2."\n" ;


  my $req2;
  my $response2;
  my $attempt=0;

  while (1) {
    print "RetrieveDataUCSC: Retrievining database $url2\n";
    $response2=retrieveURLFinalSize($self,$url2,$databaseNameSize,$fileDatabaseName);

    if ($response2->is_success) {
    } else {
      $warnLog.="Unable to retrieve $url2\n"."Status line ".$response2->status_line."\nwaiting $sleepTimeSeconds seconds to retry again\n";
      warn "Unable to retrieve $url2\n"."Status line ".$response2->status_line."\nwaiting $sleepTimeSeconds seconds to retry again\n";
      sleep($sleepTimeSeconds);
    }

    if( $contentSize == 0){
      $log.="Content is empty\n"."Status line ".$response2->status_line."\nwaiting $sleepTimeSeconds seconds to retry again\n";
      if($attempt >  2){
	$log.="Content is empty, skipping\n"."Status line ".$response2->status_line."\n";
	exit 0;
      }
      sleep($sleepTimeSeconds);
      $attempt++;
    } else {


      #my $command="/usr/bin/gunzip --force $fileDatabaseName";
      my $command=$self->{_gunzip}." --force $fileDatabaseName";
      my $output=system($command);

      if ($output != 0) {
	$log.="Command $command failed restarting\n";
      }else{
	last;
      }
    }

  }


  return ($log,$warnLog,$dieLog);
}



1;
