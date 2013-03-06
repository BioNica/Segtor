#!/usr/bin/perl

=head1 NAME

   splitFasta.pl


=head1 Disclaimer

This software and all the data therein were developed by the Bioinformatics and Computational Biology Lab (Laboratorio de Bioinformatica e Biologia Computational (LBBC)) of the Brazilian National Cancer Institute (Instituto Nacional de Cancer (INCA)) in Rio de Janeiro, RJ, Brazil. The INCA is part of the Health Ministry (Ministerio da Saude) of the Federal Government of Brazil. This software was developed by the staff of the LBBC as part of their official duties. According to federal regulations, this software is therefore public domain.

This software is provided as is and the INCA does make any guarantee regarding the accuracy/quality/reliability of the results. The INCA cannot be held accountable for any hardware or software damage that might be caused through the use or misuse of this software.

The LBBC encourages researchers using our software to acknowledge us. The software can be distributed freely and modified given that any derivative bears an acknowledgment that it was derived from our original software.

=head1 Author

Gabriel Renaud grenaud@inca.gov.br



=cut

use strict;
use warnings;


use Getopt::Std;

use Cwd;
use Cwd 'abs_path';
use File::Basename;
use Data::Dumper;

my %opts;
getopts('f:e', \%opts);

my $usage="Usage:\n".
  "./splitFasta.pl -f [inputFile] [options]\n".
  "-f [fasta file] The input fasta file\t".
  "Options:\n".
  "-e errase the original file\n";


my $inputFastaFile;
my $eraseOriginalFile=0;

if( !(exists $opts{'f'}) || !$opts{'f'}){
  print "The -f argument is mandatory\n";
  die $usage;
}else{
  $inputFastaFile=$opts{'f'};
}


if( (exists $opts{'e'}) ){
  $eraseOriginalFile=1;
}

$inputFastaFile=abs_path($inputFastaFile);

my $abs_path          = dirname(abs_path($inputFastaFile));
#my (undef,undef,$suffix) = fileparse($inputFastaFile,qr{\..*});
my $suffix;
if($inputFastaFile =~ /(\.\w+?)$/){
  $suffix=$1;
}else{
  die "Wrong pattern for file $inputFastaFile\n";
}


my $counterOfDeflines = 0;
my $filePtrOpen       = 0;
my $FILEPTR           = undef;
my %checkUniqueness;
my $foundADefWithFileName =0;
my $foundADefWithFileNameTEMPNAME =0;


open(INFO, $inputFastaFile) or die "splitFasta.pl:  Cannot open file ".$inputFastaFile."\n";
while(my $line = <INFO>){
  chomp($line);
  if(substr($line,0,1) eq ">"){
    if(index($line, " ") != -1){
      die "splitFasta.pl: The defline cannot contain white spaces\n";
    }

    my $fileNameOUT=substr($line,1);

    if (exists $checkUniqueness{ $fileNameOUT }) {
      die "splitFasta.pl: Duplicate defline ".$fileNameOUT."\n";
    } else {
      $checkUniqueness{$fileNameOUT}=1;
    }

    if(  $inputFastaFile eq $abs_path ."/". $fileNameOUT.$suffix ){
      $foundADefWithFileName=1;
      $fileNameOUT.="TEMP";
      $foundADefWithFileNameTEMPNAME=$fileNameOUT.$suffix;
    }

    if($counterOfDeflines == 0){
      open($FILEPTR,">".$abs_path ."/".$fileNameOUT.$suffix) or die "splitFasta.pl:  Cannot open file ".$abs_path ."/".$fileNameOUT.$suffix."\n";
      print $FILEPTR $line."\n";
    }else{
      close($FILEPTR);
      open($FILEPTR,">".$abs_path ."/".$fileNameOUT.$suffix) or die "splitFasta.pl: Cannot open file ".$abs_path ."/".$fileNameOUT.$suffix."\n";
      print $FILEPTR $line."\n";
    }
    $filePtrOpen=1;
    $counterOfDeflines++;
  }else{

    if(!$filePtrOpen){
      die "splitFasta.pl: waiting for a defline and found line=$line\n";
    }
    print $FILEPTR $line."\n";
  }
}
close(INFO);
close($FILEPTR);

if(!$filePtrOpen){
  die "splitFasta.pl: The file $inputFastaFile appears to be empty\n";
}



if($foundADefWithFileName){

  if($counterOfDeflines != 1){
    die "splitFasta.pl: The file $inputFastaFile contains a fasta record with the same name as the file but has more than 1 fasta record\n";
  }

  unlink($abs_path ."/". $foundADefWithFileNameTEMPNAME) or die "Cannot delete ".$abs_path ."/". $foundADefWithFileNameTEMPNAME."\n";
}else{

  if($eraseOriginalFile){
    unlink($inputFastaFile) or die "Cannot delete ".$inputFastaFile."\n";
  }
}

