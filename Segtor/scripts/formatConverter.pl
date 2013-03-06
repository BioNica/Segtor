#!/usr/bin/perl

=head1 NAME

   snvmix2Segtor.pl


=head1 Disclaimer

This software and all the data therein were developed by the Bioinformatics and Computational Biology Lab (Laboratorio de Bioinformatica e Biologia Computational (LBBC)) of the Brazilian National Cancer Institute (Instituto Nacional de Cancer (INCA)) in Rio de Janeiro, RJ, Brazil. The INCA is part of the Health Ministry (Ministerio da Saude) of the Federal Government of Brazil. This software was developed by the staff of the LBBC as part of their official duties. According to federal regulations, this software is therefore public domain.

This software is provided as is and the INCA does make any guarantee regarding the accuracy/quality/reliability of the results. The INCA cannot be held accountable for any hardware or software damage that might be caused through the use or misuse of this software.

The LBBC encourages researchers using our software to acknowledge us. The software can be distributed freely and modified given that any derivative bears an acknowledgment that it was derived from our original software.

=head1 Author

Gabriel Renaud grenaud@inca.gov.br


=cut

local $| = 1;

sub isDNA{
  my ($bp)=@_;
  $bp=uc($bp);

  if ( $bp eq "A" ||
       $bp eq "C" ||
       $bp eq "G" ||
       $bp eq "T" ) {
    return 1;
  } else {
    return 0;
  }
}


sub convertIUPAC{
  my ($bp)=@_;
  $bp=uc($bp);

  if ( $bp eq "A" ){
    return ("A");
  }elsif ( $bp eq "C" ){
    return ("C");
  }elsif ( $bp eq "G" ){
    return ("G");
  }elsif ( $bp eq "T" ){
    return ("T");

  }elsif ( $bp eq "R" ){
    return ("A","G");
  }elsif ( $bp eq "Y" ){
    return ("C","T");
  }elsif ( $bp eq "S" ){
    return ("G","C");
  }elsif ( $bp eq "W" ){
    return ("A","T");
  }elsif ( $bp eq "K" ){
    return ("G","T");
  }elsif ( $bp eq "M" ){
    return ("A","C");

  }elsif ( $bp eq "B" ){
    return ("C","G","T");
  }elsif ( $bp eq "D" ){
    return ("A","G","T");
  }elsif ( $bp eq "H" ){
    return ("A","C","T");
  }elsif ( $bp eq "V" ){
    return ("A","C","G");

  }elsif ( $bp eq "N" ){
    return ("A","C","G","T");


  } else {
    die "Invalid IUPAC code $bp\n";
  }
}


use strict;
use warnings;
use Getopt::Std;
use Data::Dumper;
use List::Util qw[min max];

my %opts;
getopts('i:o:f:t:h', \%opts);
#                   0     1            2     3        4
my @inputFormats=("BAM","BAMpairend","BED","snvmix","samtools-cv","varscan");
#                  0            1          2     3
my @outputFormats=("coordinate","interval","SNV","indel");
my $limitPEdistance=2000;
my $usage="$0\n".
  "-i [format] input format options\n".
  "\t".join("\n\t",@inputFormats)."\n".
  "\nBAM and BAMpairend require to have bamToBed from BEDTools in your path\n".
  "-o [format] output Segtor format\n".
  "\t".join("\n\t",@outputFormats)."\n".
  "-f [file] (optional) input file to read from default: STDIN \n".
  "-t [file] (optional) input file to write to default: STDOUT \n".
  "\nWarning:\n".
  "\tif you select BAMpairend, the file must be sorted by read name\n\tand only the paired reads will be reported\n";


if ( scalar(keys %opts) == 0 ) {
  print "usage:\n".$usage."";
  exit(0);
}

if ($opts{'h'}) {
  print "usage:\n".$usage."";
  exit(0);
}

my $inputFormat;
my $outputFormat;
my $inputFormatCode=-1;
my $outputFormatCode=-1;

if ($opts{'i'}) {
  my $found=0;

  for (my $i=0;$i<=$#inputFormats;$i++) {
    if ($opts{'i'} eq $inputFormats[$i]) {
      $found=1;
      $inputFormatCode=$i;
    }
  }

  if ($found == 0) {
    print "usage:\n".$usage."";
    exit(0);
  } else {
    $inputFormat = $opts{'i'};
  }
} else {
  print "usage:\n".$usage."";
  exit(0);
}






if ($opts{'o'}) {
  my $found=0;
  #foreach my $format (@outputFormats){
  for (my $i=0;$i<=$#outputFormats;$i++) {
    if ($opts{'o'} eq $outputFormats[$i]) {
      $found=1;
      $outputFormatCode=$i;
    }
  }
  if ($found == 0) {
    print "usage:\n".$usage."";
    exit(0);
  } else {
    $outputFormat = $opts{'o'};
  }

} else {
  print "usage:\n".$usage."";
  exit(0);
}




open(INPUTSTREAM, '-') or die "Can't open STDIN ";

if ($opts{'f'}) {
  my $file=$opts{'f'};
  if ($inputFormatCode ==  0 || $inputFormatCode ==  1) {
    my $command = "bamToBed -i $file";
    open(INPUTSTREAM, "$command |") or die ("Can't run command \"$command\"\n");
  } else {
    open(INPUTSTREAM, $file) or die ("Can't open $file ");
  }

} else {
  #STDIN
  if ($inputFormatCode ==  0 || $inputFormatCode ==  1) {
    my $command = "bamToBed -i /dev/stdin";
    open(INPUTSTREAM, "$command |") or die ("Can't run command \"$command\"\n");
  } else {
    #nothing to do
  }

}


open(OUTPUTSTREAM,'>-') or die "Cannot open STDOUT\n";


if ($opts{'t'}) {
  my $file=">".$opts{'t'};
  open(OUTPUTSTREAM, $file) or die ("Can't open $file ");
} else {

}


my $id=0;


my $counterINDEL=1;

my $previousChr    ;
my $previousCoord1 ;
my $previousCoord2 ;
my $previousName   ;
my $previousId     ;
my $previousScore  ;
my $previousStrand ;


while (my $line = <INPUTSTREAM>) {

  # bam paired end
  if ($inputFormatCode ==  1) {
    chomp($line);
    my @temparray=split("\t",$line);
    if ($#temparray < 5) {
      die "At least the first 5 fields must be define in the BED format\n";
    }

    my $chr    = $temparray[0];
    my $coord1 = $temparray[1];
    my $coord2 = $temparray[2];
    my $name   = $temparray[3];
    my $score  = $temparray[4];
    my $strand = $temparray[5];

    if($name =~ /^(\S+)\/1$/){
      $previousChr    = $chr ;
      $previousCoord1 = $coord1;
      $previousCoord2 = $coord2;
      $previousName   = $name;
      $previousId     = $1;
      $previousScore  = $score;
      $previousStrand = $strand;
    }elsif($name =~ /(\S+)\/2$/){
      if($1 eq $previousId){
	if ($outputFormatCode == 0) { #coordinate
	  if ($strand eq "+") {
	    print OUTPUTSTREAM $chr."\t".$previousCoord2."\t".$previousId."\n";
	  } elsif ($strand eq "-") {
	    print OUTPUTSTREAM $chr."\t".$previousCoord1."\t".$previousId."\n";
	  } else {
	    die "Wrong input/output pair\n";
	  }
	} elsif ($outputFormatCode == 1) { #interval
	  my $min=min($previousCoord1,$coord1);
	  my $max=max($previousCoord2,$coord2);
	  if($min > $max ){
	   warn "Skipping line ".$chr."\t".$min."\t".$max."\t".$previousId."\n";
	  }else{
	    print OUTPUTSTREAM $chr."\t".$min."\t".$max."\t".$previousId."\n";
	  }
	} else {
	  die "Wrong input/output pair\n";
	}
      }

    }else{
      die "The line $line did not parse or it not paired-end\n";
    }

  }



  # BAM and BED
  if ($inputFormatCode ==  0   || $inputFormatCode ==  2) {
    chomp($line);
    my @temparray=split("\t",$line);
    if ($#temparray < 5) {
      die "At least the first 5 fields must be define in the BED format\n";
    }

    my $chr    = $temparray[0];
    my $coord1 = $temparray[1];
    my $coord2 = $temparray[2];
    my $name   = $temparray[3];
    my $score  = $temparray[4];
    my $strand = $temparray[5];

    if ($outputFormatCode == 0) { #coordinate
      if ($strand eq "+") {
	print OUTPUTSTREAM $chr."\t".$coord1."\t".$name."\n";
      } elsif ($strand eq "-") {
	print OUTPUTSTREAM $chr."\t".$coord2."\t".$name."\n";
      } else {
	die "Strand $strand is wrong in line $line\n";
      }
    } elsif ($outputFormatCode == 1) { #interval
      my $min=min($coord1,$coord2);
      my $max=max($coord1,$coord2);
      print OUTPUTSTREAM $chr."\t".$min."\t".$max."\t".$name."\n";
    } else {
      die "Wrong input/output pair\n";
    }
  }


  #SNVMIX
  #chr1:567062	C	T	C:2,T:43,0.0000000000,0.0000000143,0.9999999857,3
  if ($inputFormatCode ==  3) {
    if ($outputFormatCode != 2) { #SNV
      die "Wrong input/output pair\n";
    }
    if ($line =~ /^(\S+):(\d+)\s+(\w)\s+(\w)\s+(\w):(\d+),(\w):(\d+),([\d\.]+),([\d\.]+),([\d\.]+),(\d)$/) { 
      if ( uc($3) eq 'N' ) {
	warn "Skipping line ".$line."\n";
	next;
      }
      if ( ! isDNA($3) ) {
	die "Base pair $3 in line $line did not parse\n";
      }
      if ( ! isDNA($4) ) {
	die "Base pair $4 in line $line did not parse\n";
      }
      if ( ! isDNA($5) ) {
	die "Base pair $5 in line $line did not parse\n";
      }
      if ( ! isDNA($7) ) {
	die "Base pair $7 in line $line did not parse\n";
      }

      my $tag="";
      $tag.=$5."_".$6;
      $tag.=$7."_".$8."_";

      if ($12 == 1) {
	$tag.="homoRef";
      } elsif ($12 == 2) {
	$tag.="hetero";
      } elsif ($12 == 3) {
	$tag.="homoNonRef";
      } else {
	die "Wrong zygosity in line $line did not parse\n";
      }
      $tag.=$id;
      $id++;
      if ($outputFormatCode == 0) { #coordinate
	print OUTPUTSTREAM $1."\t".$2."\t".$tag."\n";
      } elsif ($outputFormatCode == 1) { #interval
	die "Wrong input/output pair\n";
      } elsif ($outputFormatCode == 2) { #snv
	print OUTPUTSTREAM $1."\t".$2."\t".uc($3)."\t".uc($4)."\t".$tag."\n";
      } else {
	die "Wrong input/output pair\n";

      }


    } else {
      die "Line $line did not parse\n";
    }
  }				#end $inputFormatCode = 3

  #samtools cv
  if ($inputFormatCode ==  4) {
    #                 1       2      3     4
    if ($line =~ /^(\S+)\s+(\d+)\s+(\w)\s+(\w)\s+/) {
      if ($outputFormatCode == 2) { #indel
	my @arrayBPs=convertIUPAC($4);
	foreach my $bp (@arrayBPs){
	  my $tag="snv#".$id;
	  $id++;
	  print OUTPUTSTREAM $1."\t".$2."\t".uc($3)."\t".$bp."\t".$tag."\n";
	}
      }else {
	die "Wrong input/output pair\n";
      }
    }else{
      warn "Skipping ".$line."\n";
    }
  }

  #var scan
  if ($inputFormatCode ==  5) {
    #chrUn_gl000241	11828	A	-GTACCA
    if ($line =~ /^Chrom\s+/) {
      next;
    }

    if ($line =~ /^(\S+)\s+(\d+)\s+(\w)\s+([+-])(\w+)\s+/) {
      if ($outputFormatCode == 3) { #indel
	my $chr    = $1;
	my $coord  = $2;
	my $ref    = $3;
	my $type   = $4;
	my $seq    = $5;
	if ($type eq "+") {
	  print OUTPUTSTREAM "INS\t".$chr."\t".$coord."\t".$seq."\tins".$counterINDEL++."\n";
	} elsif ($type eq "-") {
	  print OUTPUTSTREAM "DEL\t".$chr."\t".($coord+1)."\t".($coord+length($seq))."\tdel".$counterINDEL++."\n";
	} else {
	  die "Line $line did not parse\n";
	}
      } else {
	die "Wrong input/output pair\n";
      }

    } else {
      die "Line $line did not parse\n";
    }
  }

}


close(INPUTSTREAM);
close(OUTPUTSTREAM);



