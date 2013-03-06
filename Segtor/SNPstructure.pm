package SNPstructure;

use strict;
use warnings;
use Data::Dumper;

#
#    513 cDNA deletion
#     77 cDNA in-del
#    279 cDNA insertion
#      4 cDNA microsatellite
#     10 cDNA mixed
#      5 cDNA mnp
#  91469 cDNA single

#1721204 genomic deletion
#      4 genomic het
# 433467 genomic in-del
#2564932 genomic insertion
#   4878 genomic microsatellite
# 103802 genomic mixed
#  68498 genomic mnp
#  37327 genomic named
#13807054 genomic single


#=item new
#
#   Constructor
#
#=

sub new{
  my ($class,%arg)=(@_);
  my ($self) =bless{_genomebinsize       =>  $arg{genomebinsize}}, $class;

  $self->{_binnedHashCoord}=[];
  $self->{_chr}=undef;

  return $self;
}

#0      1       2       3       4               5       6       7       8       9       10      11      12      13      14      15
#585	chr1	582	583	rs58108140	0	+	G	G	A/G	genomic	single	unknown	0	0	unknown	exact	1

sub addLine{
  my ($self,$line)=@_;

  my @arrayOfFields=split("\t",$line);

  if ($self->{_chr}) {
    if ($self->{_chr} ne $arrayOfFields[1]) {
      die "SNPstructure.pm the chromsomes cannot be different\n";
    }
  } else {
    $self->{_chr}=$arrayOfFields[1];
  }

  if($arrayOfFields[11] eq "single"){
    my $hashToInsert;
    my $arrayOfObserved=[];
    my $bpRef;

    if($arrayOfFields[7] ne $arrayOfFields[8]){
      return "SNPstructure.pm fields 7 & 8 should be the same for snpid = ".$arrayOfFields[4]."\n";
      # die;
    }
    if(($arrayOfFields[2]+1) != $arrayOfFields[3]){
      return  "SNPstructure.pm  2 (minus 1) & 3 should be the same for ".$arrayOfFields[4]."\n";
    }

    if($arrayOfFields[6] eq "+"){
      foreach my $obsBP (split("/",$arrayOfFields[9])){
	push(@{$arrayOfObserved},$obsBP);
      }
      $bpRef=$arrayOfFields[7];

      foreach my $obsBP (split("/",$arrayOfFields[9])){
	if($obsBP     eq "A"){
	}elsif($obsBP eq "C"){
	}elsif($obsBP eq "G"){
	}elsif($obsBP eq "T"){
	}else{
	  return "SNPstructure.pm wrong base pair, field #9 =  ".$arrayOfFields[9]."  for snpid = ".$arrayOfFields[4]."\n";
	}
      }

      if ($arrayOfFields[7]      eq "A") {
      } elsif ($arrayOfFields[7] eq "C") {
      } elsif ($arrayOfFields[7] eq "G") {
      } elsif ($arrayOfFields[7] eq "T") {
      } else {
	return "SNPstructure.pm wrong reference base pair for ".$arrayOfFields[7]." for snpid = ".$arrayOfFields[4]."\n";
      }

    }elsif($arrayOfFields[6] eq "-"){
      foreach my $obsBP (split("/",$arrayOfFields[9])){
	if($obsBP     eq "A"){
	  push(@{$arrayOfObserved},"T");
	}elsif($obsBP eq "C"){
	  push(@{$arrayOfObserved},"G");
	}elsif($obsBP eq "G"){
	  push(@{$arrayOfObserved},"C");
	}elsif($obsBP eq "T"){
	  push(@{$arrayOfObserved},"A");
	}else{
	  return "SNPstructure.pm wrong base pair, field #9 =  ".$arrayOfFields[9]."  for snpid = ".$arrayOfFields[4]."\n";
	}
      }

      if ($arrayOfFields[7]      eq "A") {
	$bpRef="T";
      } elsif ($arrayOfFields[7] eq "C") {
	$bpRef="G";
      } elsif ($arrayOfFields[7] eq "G") {
	$bpRef="C";
      } elsif ($arrayOfFields[7] eq "T") {
	$bpRef="A";
      } else {
	return "SNPstructure.pm wrong reference base pair for ".$arrayOfFields[7]." for snpid = ".$arrayOfFields[4]."\n";
      }

    } else {
      print "SNPstructure.pm Strand ".$arrayOfFields[6]." is wrong for snpid = ".$arrayOfFields[4]."\n";
      die;
    }

    $hashToInsert={'id'    => $arrayOfFields[4] ,
		   'bpRef' => $arrayOfFields[7],
		   'bpObs' => $arrayOfObserved
		  };

    my $binIndex = int($arrayOfFields[3]/$self->{_genomebinsize} );

    if($self->{_binnedHashCoord}->[$binIndex]){ #bin has been initialized

      if(exists $self->{_binnedHashCoord}->[$binIndex]->{$arrayOfFields[3]}){ #cooord already exists
	$self->{_binnedHashCoord}->[$binIndex]->{$arrayOfFields[3]}->{'id'}.="#".$arrayOfFields[4];

	if($self->{_binnedHashCoord}->[$binIndex]->{$arrayOfFields[3]}->{'bpRef'} ne $hashToInsert->{'bpRef'}){
	  print "SNPstructure.pm Conflicting reference base pair for snpid = ".$arrayOfFields[4]."\n";
	  die;
	}


	foreach my $bp1 (@{$arrayOfObserved}) {
	  my $found=0;
	  foreach my $bp2 (@{$self->{_binnedHashCoord}->[$binIndex]->{$arrayOfFields[3]}->{'bpObs'}}) {
	    if ($bp1 eq $bp2) {
	      $found=1;
	    }
	  }
	  if (!$found) {
	    push(@{$self->{_binnedHashCoord}->[$binIndex]->{$arrayOfFields[3]}->{'bpObs'}},$bp1);
	  }
	}

      }else{ #new coord
	$self->{_binnedHashCoord}->[$binIndex]->{$arrayOfFields[3]} = $hashToInsert;
      }
    }else{ #new bin
      $self->{_binnedHashCoord}->[$binIndex]={$arrayOfFields[3] => $hashToInsert};
    }

  } #end if single


  return "";
}



sub annotateCoord{
  my ($self,$coord,$bpRef,$bpRead,$noRefBP)=@_;
  my $hashResult={};

  my $binIndex = int($coord/$self->{_genomebinsize} );

  if (exists $self->{_binnedHashCoord}->[$binIndex]->{$coord} ) {
    $hashResult->{'foundCoord'} = 1;


    my $snpRecord = $self->{_binnedHashCoord}->[$binIndex]->{$coord};

    my $hash={'snpid' => $snpRecord->{'id'}};
    my $found=0;
    if($noRefBP){
      #do not check ref base pair
    }else{
      if ($bpRef ne $snpRecord->{'bpRef'}) {
	die "The refence base pairs do not coincide\n";
      }
    }

    foreach my $obsBPSNP ( @{$snpRecord->{'bpObs'}} ) {
      if ($bpRead eq $obsBPSNP) {
	$found=1;
      }
    }

    $hash->{'alreadyFoundInDBSNP'} = $found;
    $hashResult->{'foundSNPs'} = $hash;

  } else {
    $hashResult->{'foundCoord'}  = 0;
  }
  return $hashResult;
}






1;
