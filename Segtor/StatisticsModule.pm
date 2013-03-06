package StatiticsModule;

use strict;
use warnings;
use Data::Dumper;

use POSIX qw(ceil floor);


sub new{
  my ($class,%arg)=(@_);
  my ($self) =bless{_numberOfSites    =>  0,
		    _upstream         =>  0,
		    _downstream       =>  0,
		    _exon             =>  0,
		    _intron           =>  0,
		    _exonIndices      =>  [],
		    _intronIndices    =>  [],
		    _exonIndicesMAX   =>  0,
		    _intronIndicesMAX =>  0,


		    _inTransSingleExon      =>  0,
		    _inTransExonicIntronic  =>  0,
		    _inTransIntronic        =>  0,
		    _rangeContainsTrans     =>  0,
		    _5PendFirstExon         =>  0,
		    _5PendExonsIntrons      =>  0,
		    _3PLastExon             =>  0,
		    _3PExonsIntrons         =>  0,

		    _SNP              =>  0,

		    _reliable         =>  0,
		    _unreliable       =>  0,

		    _exoncodingSNP    =>  0,
		    _exonnonCodingSNP =>  0,

		    _synSNP           =>  0,
		    _nonSynSNP        =>  0,

		    _nohits           =>  0,

		    _usedbsnp         =>  0,
		    _notindbsnp       =>  0,
		    _indbsnpdiff      =>  0,
		    _indbsnpsame      =>  0,

		    _useOfCustomDatabase     =>  0,
		    _customUpstream          =>  0,
		    _customDownstream        =>  0,
		    _customContained         =>  0,
		    _customOverlaps5p        =>  0,
		    _customOverlaps3p        =>  0,
		    _customSpans             =>  0,


		    _genomeUsed       =>  $arg{genomeUsed},
		    _databaseUsed     =>  $arg{databaseUsed},
		    _databaseFiles    =>  $arg{databaseFiles},
		    _range            =>  $arg{range},
		    _mode             =>  $arg{mode} }, $class;

  if($arg{mode} == 4){
    $self->{_bpInBins}   = $arg{bpInBins};
    $self->{_numberBins} = $arg{numberBins};

    $self->{_binsPositives}=[];
    $self->{_binsNegatives}=[];


    for(my $index=0;$index<$self->{_numberBins};$index++){
      $self->{_binsPositives}->[$index]=0;
      $self->{_binsNegatives}->[$index]=0;
    }

  }

  return $self;
}

sub useOfCustomDatabase{
  my ($self)=@_;
  $self->{_useOfCustomDatabase}=1;
}

sub addCustomUpstream{
  my ($self)=@_;
  $self->{_customUpstream}++;
}

sub addCustomDownstream{
  my ($self)=@_;
  $self->{_customDownstream}++;
}

sub addCustomContained{
  my ($self)=@_;
  $self->{_customContained}++;
}

sub addCustomOverlaps5p{
  my ($self)=@_;
  $self->{_customOverlaps5p}++;
}

sub addCustomOverlaps3p{
  my ($self)=@_;
  $self->{_customOverlaps3p}++;
}

sub addCustomSpans{
  my ($self)=@_;
  $self->{_customSpans}++;
}




sub addSite{
  my ($self)=@_;
  $self->{_numberOfSites}++;
}

sub addToNoHits{
  my ($self)=@_;
  $self->{_nohits}++;
}


sub addDBSNP{
  my ($self,$hashSNP)=@_;

  $self->{_usedbsnp}=1;
  if($hashSNP->{'foundCoord'} == 0){ #nothing in dbsnp
    $self->{_notindbsnp}++;
  }else{
    if($hashSNP->{'foundSNPs'}->{'alreadyFoundInDBSNP'} == 0){ #new bp in dbsnp
      $self->{_indbsnpdiff}++;
    }else{
      $self->{_indbsnpsame}++;
    }
  }

}


sub addClosestDistance{
  my ($self,$distance)=@_;

  if ($distance < 0) {
    $distance=-1*$distance;
    if ($distance > (( $self->{_numberBins}-1 )*$self->{_bpInBins} +1) ) {
      $self->{_binsNegatives}->[$self->{_numberBins}-1]++;
    } else {
      $self->{_binsNegatives}->[floor($distance/$self->{_bpInBins})]++;
    }
  } else {
    if ($distance > (( $self->{_numberBins}-1 )*$self->{_bpInBins} +1) ) {
      $self->{_binsPositives}->[$self->{_numberBins}-1]++;
    } else {
      $self->{_binsPositives}->[floor($distance/$self->{_bpInBins})]++;
    }
  }
}


sub addAnnotation{
  my ($self,$annotation)=@_;

  my $code=$annotation->{'code'};

  if (     $code  == 1) {
    $self->{_exon}++;

    if(exists $annotation->{"index"} ){
      if($self->{_exonIndicesMAX} <  $annotation->{"index"} ){
	$self->{_exonIndicesMAX} =  $annotation->{"index"};
      }

      if(defined($self->{_exonIndices}->[$annotation->{"index"}])){
	$self->{_exonIndices}->[$annotation->{"index"}]++;
      }else{
	$self->{_exonIndices}->[$annotation->{"index"}]=1;
      }
    }else{
      die "StatisticModules, addAnnotation cannot find exon index\n";
    }

  } elsif ( $code  == 2) {
    $self->{_intron}++;

    if(exists $annotation->{"index"} ){
      if($self->{_intronIndicesMAX} <  $annotation->{"index"} ){
	$self->{_intronIndicesMAX} =  $annotation->{"index"};
      }

      if(defined($self->{_intronIndices}->[$annotation->{"index"}])){
	$self->{_intronIndices}->[$annotation->{"index"}]++;
      }else{
	$self->{_intronIndices}->[$annotation->{"index"}]=1;
      }
    }else{
      die "StatisticModules, addAnnotation cannot find intron index\n";
    }

  } elsif ( $code  == 3) {
    $self->{_upstream}++;
  } elsif ( $code  == 4) {
    $self->{_downstream}++;
  } else {
    die "StatisticModules.pm in addAnnotation, wrong code\n";
  }
}

sub addAnnotationInterval{
  my ($self,$annotation)=@_;

  my $code=$annotation->{'code'};

  if (      $code  == 1) {
    $self->{_inTransSingleExon}++;

    if(exists $annotation->{"indicesEx"} ){
      foreach my $exonIndex (split(",",$annotation->{"indicesEx"})){
	if($self->{_exonIndicesMAX} <  $exonIndex ){
	  $self->{_exonIndicesMAX} =  $exonIndex;
	}

	if(defined($self->{_exonIndices}->[$exonIndex])){
	  $self->{_exonIndices}->[$exonIndex]++;
	}else{
	  $self->{_exonIndices}->[$exonIndex]=1;
	}
      }
    }else{
      die "StatisticModules, addAnnotationInterval cannot find exon index\n";
    }

  } elsif ( $code  == 2) {
    $self->{_inTransExonicIntronic}++;
    if(exists $annotation->{"indicesEx"} ){
      foreach my $exonIndex (split(",",$annotation->{"indicesEx"})){

	if($self->{_exonIndicesMAX} <  $exonIndex ){
	  $self->{_exonIndicesMAX} =  $exonIndex;
	}

	if(defined($self->{_exonIndices}->[$exonIndex])){
	  $self->{_exonIndices}->[$exonIndex]++;
	}else{
	  $self->{_exonIndices}->[$exonIndex]=1;
	}
      }
    }else{
      die "StatisticModules, addAnnotationInterval cannot find exon index\n";
    }

    if(exists $annotation->{"indicesIntron"} ){
      foreach my $intronIndex (split(",",$annotation->{"indicesIntron"})){
	if($self->{_intronIndicesMAX} <  $intronIndex ){
	  $self->{_intronIndicesMAX} =  $intronIndex;
	}

	if(defined($self->{_intronIndices}->[$intronIndex])){
	  $self->{_intronIndices}->[$intronIndex]++;
	}else{
	  $self->{_intronIndices}->[$intronIndex]=1;
	}
      }
    }else{
      die "StatisticModules, addAnnotationInterval cannot find intron index\n";
    }

  } elsif ( $code  == 3) {
    $self->{_inTransIntronic}++;
    if(exists $annotation->{"indicesIntron"} ){
      foreach my $intronIndex (split(",",$annotation->{"indicesIntron"})){
	if($self->{_intronIndicesMAX} <  $intronIndex ){
	  $self->{_intronIndicesMAX} =  $intronIndex;
	}

	if(defined($self->{_intronIndices}->[$intronIndex])){
	  $self->{_intronIndices}->[$intronIndex]++;
	}else{
	  $self->{_intronIndices}->[$intronIndex]=1;
	}
      }
    }else{
      die "StatisticModules, addAnnotationInterval cannot find intron index\n";
    }

  } elsif ( $code  == 4) {
    $self->{_rangeContainsTrans}++;
  } elsif ( $code  == 5) {
    $self->{_5PendFirstExon}++;

    if(exists $annotation->{"indicesEx"} ){
      foreach my $exonIndex (split(",",$annotation->{"indicesEx"})){
	if($self->{_exonIndicesMAX} <  $exonIndex ){
	  $self->{_exonIndicesMAX} =  $exonIndex;
	}

	if(defined($self->{_exonIndices}->[$exonIndex])){
	  $self->{_exonIndices}->[$exonIndex]++;
	}else{
	  $self->{_exonIndices}->[$exonIndex]=1;
	}
      }
    }else{
      die "StatisticModules, addAnnotationInterval cannot find exon index\n";
    }

  } elsif ( $code  == 6) {
    $self->{_5PendExonsIntrons}++;
    if(exists $annotation->{"indicesEx"} ){
      foreach my $exonIndex (split(",",$annotation->{"indicesEx"})){
	if($self->{_exonIndicesMAX} <  $exonIndex ){
	  $self->{_exonIndicesMAX} =  $exonIndex;
	}

	if(defined($self->{_exonIndices}->[$exonIndex])){
	  $self->{_exonIndices}->[$exonIndex]++;
	}else{
	  $self->{_exonIndices}->[$exonIndex]=1;
	}
      }
    }else{
      die "StatisticModules, addAnnotationInterval cannot find exon index\n";
    }

    if(exists $annotation->{"indicesIntron"} ){
      foreach my $intronIndex (split(",",$annotation->{"indicesIntron"})){
	if($self->{_intronIndicesMAX} <  $intronIndex ){
	  $self->{_intronIndicesMAX} =  $intronIndex;
	}

	if(defined($self->{_intronIndices}->[$intronIndex])){
	  $self->{_intronIndices}->[$intronIndex]++;
	}else{
	  $self->{_intronIndices}->[$intronIndex]=1;
	}
      }
    }else{
      die "StatisticModules, addAnnotationInterval cannot find intron index\n";
    }

  } elsif ( $code  == 7) {
    $self->{_3PLastExon}++;
    if(exists $annotation->{"indicesEx"} ){
      foreach my $exonIndex (split(",",$annotation->{"indicesEx"})){
	if($self->{_exonIndicesMAX} <  $exonIndex ){
	  $self->{_exonIndicesMAX} =  $exonIndex;
	}

	if(defined($self->{_exonIndices}->[$exonIndex])){
	  $self->{_exonIndices}->[$exonIndex]++;
	}else{
	  $self->{_exonIndices}->[$exonIndex]=1;
	}
      }
    }else{
      die "StatisticModules, addAnnotationInterval cannot find exon index\n";
    }

  } elsif ( $code  == 8) {
    $self->{_3PExonsIntrons}++;
    if(exists $annotation->{"indicesEx"} ){
      foreach my $exonIndex (split(",",$annotation->{"indicesEx"})){
	if($self->{_exonIndicesMAX} <  $exonIndex ){
	  $self->{_exonIndicesMAX} =  $exonIndex;
	}

	if(defined($self->{_exonIndices}->[$exonIndex])){
	  $self->{_exonIndices}->[$exonIndex]++;
	}else{
	  $self->{_exonIndices}->[$exonIndex]=1;
	}
      }
    }else{
      die "StatisticModules, addAnnotationInterval cannot find exon index\n";
    }

    if(exists $annotation->{"indicesIntron"} ){
      foreach my $intronIndex (split(",",$annotation->{"indicesIntron"})){
	if($self->{_intronIndicesMAX} <  $intronIndex ){
	  $self->{_intronIndicesMAX} =  $intronIndex;
	}

	if(defined($self->{_intronIndices}->[$intronIndex])){
	  $self->{_intronIndices}->[$intronIndex]++;
	}else{
	  $self->{_intronIndices}->[$intronIndex]=1;
	}
      }
    }else{
      die "StatisticModules, addAnnotationInterval cannot find intron index\n";
    }

  } elsif ( $code  == 9) {
    $self->{_upstream}++;
  } elsif ( $code  == 10) {
    $self->{_downstream}++;
  } else {
    die "StatisticModules.pm in addAnnotationInterval wrong code\n";
  }

}


sub addAnnotationSNP{
  my ($self,$annotation,$snpAnnotation)=@_;
  $self->{_SNP}=1;
  $self->addAnnotation($annotation);



  if ($snpAnnotation->{'reliable'} == 1) {
    $self->{_reliable}++;

    if ($snpAnnotation->{'code'} == 1 ) {
      if ($snpAnnotation->{'coding'} == 1 ) {
	$self->{_exoncodingSNP}++;


	if ( $snpAnnotation->{'synonymous'} ) {
	  $self->{_synSNP}++
	} else {
	  $self->{_nonSynSNP}++;
	}


      } else {
	$self->{_exonnonCodingSNP}++;
      }
    }

  } else {
    $self->{_unreliable}++;
  }

}

sub printBasicStats{
  my ($self,$message)=@_;
  my $stringToReturn="";
  $stringToReturn.="Annotation Statistics".$message."\n";
  $stringToReturn.="Genome Used      : ".$self->{_genomeUsed}."\n";
  $stringToReturn.="Database Used    : ".$self->{_databaseUsed}."\n";
  if($self->{_mode} != 4){
    $stringToReturn.="Range Used       : ".$self->{_range}."\n";
  }
  for(my $index=0;$index<$#{ $self->{_databaseFiles} };$index+=2){
    $stringToReturn.="Database file : ". $self->{_databaseFiles}->[$index]." downloaded on ". $self->{_databaseFiles}->[$index+1]."\n";
  }

  return $stringToReturn;
}


sub printStats{
  my ($self,$message)=@_;


  my $stringToReturn="";
  $stringToReturn.=$self->printBasicStats($message);

  if ($self->{_useOfCustomDatabase}) {
    if (    $self->{_mode} == 1) { #coordinates
      $stringToReturn.="\nTotal number of sites       : ".$self->{_numberOfSites}."\n";
      $stringToReturn.="Coordinates upstream          : ".$self->{_customUpstream}." ( ".100*$self->{_customUpstream}/$self->{_numberOfSites}."%)\n";
      $stringToReturn.="Coordinates downstream        : ".$self->{_customDownstream}." ( ".100*$self->{_customDownstream}/$self->{_numberOfSites}."%)\n";
      $stringToReturn.="Coordinates contained         : ".$self->{_customContained}." ( ".100*$self->{_customContained}/$self->{_numberOfSites}."%)\n";
      $stringToReturn.="Coordinates without hits      : ".$self->{_nohits}." ( ".100*$self->{_nohits}/$self->{_numberOfSites}."%)\n";


    }elsif ($self->{_mode} == 2) { #intervals



      $stringToReturn.="\nTotal number of sites       : ".$self->{_numberOfSites}."\n";
      $stringToReturn.="Intervals upstream          : ".$self->{_customUpstream}." ( ".100*$self->{_customUpstream}/$self->{_numberOfSites}."%)\n";
      $stringToReturn.="Intervals downstream        : ".$self->{_customDownstream}." ( ".100*$self->{_customDownstream}/$self->{_numberOfSites}."%)\n";
      $stringToReturn.="Intervals contained         : ".$self->{_customContained}." ( ".100*$self->{_customContained}/$self->{_numberOfSites}."%)\n";
      $stringToReturn.="Intervals spanning the 5' end         : ".$self->{_customOverlaps5p}." ( ".100*$self->{_customOverlaps5p}/$self->{_numberOfSites}."%)\n";
      $stringToReturn.="Intervals spanning the 3' end         : ".$self->{_customOverlaps3p}." ( ".100*$self->{_customOverlaps3p}/$self->{_numberOfSites}."%)\n";
      $stringToReturn.="Intervals spanning an entire record   : ".$self->{_customSpans}." ( ".100*$self->{_customSpans}/$self->{_numberOfSites}."%)\n";

      $stringToReturn.="Intervals without hits      : ".$self->{_nohits}." ( ".100*$self->{_nohits}/$self->{_numberOfSites}."%)\n";




    }else{
      die "StatisticModules, cannot use custom database for any other mode beside 1 and 2\n";
    }
  } else {

    if ($self->{_mode} != 4) {
      $stringToReturn.="\nTotal number of sites       : ".$self->{_numberOfSites}."\n";
      $stringToReturn.="Sites without hits           : ".$self->{_nohits}." ( ".100*$self->{_nohits}/$self->{_numberOfSites}."%)\n";

      if ($self->{_mode} == 2) { #intervals


	$stringToReturn.="Intervals inside a single exon                       : ".$self->{_inTransSingleExon}." ( ".100*$self->{_inTransSingleExon}/$self->{_numberOfSites}."%)\n";
	$stringToReturn.="Intervals overlapping both exons/introns             : ".$self->{_inTransExonicIntronic}." ( ".100*$self->{_inTransExonicIntronic}/$self->{_numberOfSites}."%)\n";
	$stringToReturn.="Intervals inside a single intron                     : ".$self->{_inTransIntronic}." ( ".100*$self->{_inTransIntronic}/$self->{_numberOfSites}."%)\n";
	$stringToReturn.="Intervals spanning an entire transcript              : ".$self->{_rangeContainsTrans}." ( ".100*$self->{_rangeContainsTrans}/$self->{_numberOfSites}."%)\n";
	$stringToReturn.="Intervals overlapping the 5' end and first exon      : ".$self->{_5PendFirstExon}." ( ".100*$self->{_5PendFirstExon}/$self->{_numberOfSites}."%)\n";
	$stringToReturn.="Intervals overlapping the 5' end and exons/introns   : ".$self->{_5PendExonsIntrons}." ( ".100*$self->{_5PendExonsIntrons}/$self->{_numberOfSites}."%)\n";
	$stringToReturn.="Intervals overlapping the 3' end and last exon       : ".$self->{_3PLastExon}." ( ".100*$self->{_3PLastExon}/$self->{_numberOfSites}."%)\n";
	$stringToReturn.="Intervals overlapping the 3' end and exons/introns   : ".$self->{_3PExonsIntrons}." ( ".100*$self->{_3PExonsIntrons}/$self->{_numberOfSites}."%)\n";

	$stringToReturn.="Intervals upstream of genes                          : ".$self->{_upstream}." ( ".100*$self->{_upstream}/$self->{_numberOfSites}."%)\n";
	$stringToReturn.="Intervals downstream of genes                        : ".$self->{_downstream}." ( ".100*$self->{_downstream}/$self->{_numberOfSites}."%)\n";
      } else {
	$stringToReturn.="Sites upstream of genes       : ".$self->{_upstream}." ( ".100*$self->{_upstream}/$self->{_numberOfSites}."%)\n";
	$stringToReturn.="Sites downstream of genes     : ".$self->{_downstream}." ( ".100*$self->{_downstream}/$self->{_numberOfSites}."%)\n";
	$stringToReturn.="Sites within exons of genes   : ".$self->{_exon}." ( ".100*$self->{_exon}/$self->{_numberOfSites}."%)\n";
	$stringToReturn.="Sites within introns of genes : ".$self->{_intron}." ( ".100*$self->{_intron}/$self->{_numberOfSites}."%)\n";
      }
      $stringToReturn.="Exons:\n";

      for (my $index=1;$index<= $self->{_exonIndicesMAX} ; $index++) {
	if (defined($self->{_exonIndices}->[$index])) {
	  $stringToReturn.="Exon#".($index).":".$self->{_exonIndices}->[$index]."\n";
	}
      }

      $stringToReturn.="Introns:\n";
      for (my $index=1;$index<=  $self->{_intronIndicesMAX} ; $index++) {
	if (defined($self->{_intronIndices}->[$index])) {
	  $stringToReturn.="Intron#".($index).":".$self->{_intronIndices}->[$index]."\n";
	}
      }

      if ( $self->{_SNP} ) {
	$stringToReturn.="SNVs hitting reliable genes                : ".$self->{_reliable}."\n";
	$stringToReturn.="SNVs hitting unreliable genes              : ".$self->{_unreliable}."\n";
	$stringToReturn.="SNVs hitting coding exons                  : ".$self->{_exoncodingSNP}."\n";
	$stringToReturn.="SNVs hitting non-coding exons              : ".$self->{_exonnonCodingSNP}."\n";
	$stringToReturn.="SNVs causing non-synonymous mutations      : ".$self->{_nonSynSNP}."\n";
	$stringToReturn.="SNVs hitting synonymous mutations          : ".$self->{_synSNP}."\n";

	if ($self->{_usedbsnp}) {
	  $stringToReturn.="SNVs not in dbSNP                          : ".$self->{_notindbsnp}."\n";
	  $stringToReturn.="SNVs in dbSNP with same BPs                : ".$self->{_indbsnpsame}."\n";
	  $stringToReturn.="SNVs in dbSNP with different BPs           : ".$self->{_indbsnpdiff}."\n";
	}
      }
    } else {
      $stringToReturn.="\nPositive distances\n";
      my $minBp=0;
      my $maxBp=$self->{_bpInBins}-1;

      for (my $index=0;$index<$self->{_numberBins};$index++) {
	if ($index == ($self->{_numberBins}-1) ) {
	  $stringToReturn.="$minBp and up\t".$self->{_binsPositives}->[$index]."\n";
	} else {
	  $stringToReturn.="$minBp to $maxBp\t".$self->{_binsPositives}->[$index]."\n";
	}

	$minBp+=$self->{_bpInBins};
	$maxBp+=$self->{_bpInBins};
      }

      $stringToReturn.="\nNegative distances\n";
      $minBp=0;
      $maxBp=$self->{_bpInBins}-1;

      for (my $index=0;$index<$self->{_numberBins};$index++) {
	if ($index ==0 ) {
	  $stringToReturn.="$minBp";
	} else {
	  $stringToReturn.="-$minBp";
	}
	if ($index == ($self->{_numberBins}-1) ) {
	  $stringToReturn.=" and up\t".$self->{_binsNegatives}->[$index]."\n";
	} else {
	  $stringToReturn.=" to -$maxBp\t".$self->{_binsNegatives}->[$index]."\n";
	}

	$minBp+=$self->{_bpInBins};
	$maxBp+=$self->{_bpInBins};
      }

    }
  }
  return $stringToReturn;
}





sub printBasicStatsXML{
  my ($self)=@_;
  my $stringToReturn="";

  $stringToReturn.="<GENOMEUSED>".$self->{_genomeUsed}."</GENOMEUSED>\n";
  $stringToReturn.="<DATABASEUSED>".$self->{_databaseUsed}."</DATABASEUSED>\n";
  if($self->{_mode} != 4){
    $stringToReturn.="<RANGEUSED>".$self->{_range}."</RANGEUSED>\n";
  }
  $stringToReturn.="<DATABASE>\n";
  for(my $index=0;$index<$#{ $self->{_databaseFiles} };$index+=2){
    $stringToReturn.="<DATABASEFILE name=\"". $self->{_databaseFiles}->[$index]." downloaded on ". $self->{_databaseFiles}->[$index+1]."\"/>"."\n";
  }
  $stringToReturn.="</DATABASE>\n";

  return $stringToReturn;
}

sub printStatsXML{
  my ($self)=@_;


  my $stringToReturn="<STATS>\n";

  $stringToReturn.=$self->printBasicStatsXML();
  if ($self->{_mode} != 4) {

    $stringToReturn.="<TOTALSITES>".$self->{_numberOfSites}."</TOTALSITES>\n";
    $stringToReturn.="<UPSTREAM>".$self->{_upstream}."</UPSTREAM>\n";
    $stringToReturn.="<DOWNSTREAM>".$self->{_downstream}."</DOWNSTREAM>\n";

    if ($self->{_mode} == 2) {	#intervals

      $stringToReturn.="<INTRANSSINGLEEXON>".$self->{_inTransSingleExon}."</INTRANSSINGLEEXON>\n";
      $stringToReturn.="<INTRANSEXONICINTRONIC>".$self->{_inTransExonicIntronic}."</INTRANSEXONICINTRONIC>\n";
      $stringToReturn.="<INTRANSINTRONIC>".$self->{_inTransIntronic}."</INTRANSINTRONIC>\n";
      $stringToReturn.="<RANGECONTAINSTRANS>".$self->{_rangeContainsTrans}."</RANGECONTAINSTRANS>\n";
      $stringToReturn.="<FIVEPENDFIRSTEXON>".$self->{_5PendFirstExon}."</FIVEPENDFIRSTEXON>\n";
      $stringToReturn.="<FIVEPENDEXONSINTRONS>".$self->{_5PendExonsIntrons}."</FIVEPENDEXONSINTRONS>\n";
      $stringToReturn.="<THREEPLASTEXON>".$self->{_3PLastExon}."</THREEPLASTEXON>\n";
      $stringToReturn.="<THREEPEXONSINTRONS>".$self->{_3PExonsIntrons}."</THREEPEXONSINTRONS>\n";
    }else{
      $stringToReturn.="<EXONS>".$self->{_exon}."</EXONS>\n";
      $stringToReturn.="<INTRONS>".$self->{_intron}."</INTRONS>\n";
    }

    $stringToReturn.="<EXONSINDEX>\n";

    for (my $index=1;$index<= $self->{_exonIndicesMAX} ; $index++) {
      if (defined($self->{_exonIndices}->[$index])) {
	$stringToReturn.="<EXON NUMBER=\"".($index)."\" VALUE=\"".$self->{_exonIndices}->[$index]."\"\/>"."\n";
      }
    }
    $stringToReturn.="</EXONSINDEX>\n";
    $stringToReturn.="<INTRONSINDEX>\n";

    for (my $index=1;$index<=  $self->{_intronIndicesMAX} ; $index++) {
      if (defined($self->{_intronIndices}->[$index])) {
	$stringToReturn.="<INTRON NUMBER=\"".($index)."\" VALUE=\"".$self->{_intronIndices}->[$index]."\"\/>"."\n";
      }
    }
    $stringToReturn.="</INTRONSINDEX>\n";

    if ( $self->{_SNP} ) {
      $stringToReturn.="<RELIABLE>".$self->{_reliable}."</RELIABLE>"."\n";
      $stringToReturn.="<UNRELIABLE>".$self->{_unreliable}."</UNRELIABLE>"."\n";
      $stringToReturn.="<CODINGEXONS>".$self->{_exoncodingSNP}."</CODINGEXONS>"."\n";
      $stringToReturn.="<NONCODINGEXONS>".$self->{_exonnonCodingSNP}."</NONCODINGEXONS>"."\n";
      $stringToReturn.="<NONSYN>".$self->{_nonSynSNP}."</NONSYN>"."\n";
      $stringToReturn.="<SYN>".$self->{_synSNP}."</SYN>"."\n";
    }

    $stringToReturn.="</STATS>\n";
  } else {

    $stringToReturn.="<DISTANCES>\n";
    my $minBp=0;
    my $maxBp=$self->{_bpInBins}-1;
    my @arraydistemp;
    for(my $index=0;$index<$self->{_numberBins};$index++){
      my $disstr="";

      $disstr.="<DISTANCE ";
      my $label;

      if($index ==0 ){
	$disstr.=" min=\"".$minBp."\" ";
	$label=$minBp;
      }else{
	$disstr.=" min=\"".-$minBp."\" ";
	$label=-$minBp;
      }

      if($index == ($self->{_numberBins}-1) ){
	$disstr.=" label=\"$label and up\" ";
	$disstr.=" max=\"up\" ";
	$disstr.=" bin=\"".$self->{_binsNegatives}->[$index]."\" ";
      }else{
	$disstr.=" label=\"$label to ".-$maxBp."\" ";
	$disstr.=" max=\"".-$maxBp."\" ";
	$disstr.=" bin=\"".$self->{_binsNegatives}->[$index]."\" ";
      }

      $minBp+=$self->{_bpInBins};
      $maxBp+=$self->{_bpInBins};
      $disstr.="/>";
      push(@arraydistemp,$disstr);
    }
    $stringToReturn.=join("\n",reverse(@arraydistemp))."\n";
    $stringToReturn.="<DISTANCE min=\"TSS\" label=\"TSS\" max=\"\" />\n";

    $minBp=0;
    $maxBp=$self->{_bpInBins}-1;

    for(my $index=0;$index<$self->{_numberBins};$index++){
      $stringToReturn.="<DISTANCE ";
      if($index == ($self->{_numberBins}-1) ){
	$stringToReturn.=" min=\"".$minBp."\" ";
	$stringToReturn.=" label=\"$minBp and up\" ";
	$stringToReturn.=" max=\"up\" ";
	$stringToReturn.=" bin=\"".$self->{_binsPositives}->[$index]."\" ";
      }else{
	$stringToReturn.=" min=\"".$minBp."\" ";
	$stringToReturn.=" label=\"$minBp to $maxBp\" ";
	$stringToReturn.=" max=\"".$maxBp."\" ";
	$stringToReturn.=" bin=\"".$self->{_binsPositives}->[$index]."\" ";
      }

      $minBp+=$self->{_bpInBins};
      $maxBp+=$self->{_bpInBins};
      $stringToReturn.="/>\n";
    }
    $stringToReturn.="</DISTANCES>\n";

    $stringToReturn.="</STATS>\n";
  }

  return $stringToReturn;

}


1;
