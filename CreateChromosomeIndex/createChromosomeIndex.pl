#!/usr/bin/perl

=head1 NAME

   CreateChromosomeIndex::createChromosomeIndex.pl


=head1 Disclaimer

This software and all the data therein were developed by the Bioinformatics and Computational Biology Lab (Laboratorio de Bioinformatica e Biologia Computational (LBBC)) of the Brazilian National Cancer Institute (Instituto Nacional de Cancer (INCA)) in Rio de Janeiro, RJ, Brazil. The INCA is part of the Health Ministry (Ministerio da Saude) of the Federal Government of Brazil. This software was developed by the staff of the LBBC as part of their official duties. According to federal regulations, this software is therefore public domain.

This software is provided as is and the INCA does make any guarantee regarding the accuracy/quality/reliability of the results. The INCA cannot be held accountable for any hardware or software damage that might be caused through the use or misuse of this software.

The LBBC encourages researchers using our software to acknowledge us. The software can be distributed freely and modified given that any derivative bears an acknowledgment that it was derived from our original software.

=head1 Author

Gabriel Renaud grenaud@inca.gov.br



=cut



use strict;
use warnings;

use Data::Dumper;


use Cwd         qw(abs_path);

my $path;

BEGIN {
  my @pathtemp = split("/", abs_path($0) );
  delete($pathtemp[$#pathtemp]);
  delete($pathtemp[$#pathtemp]);
  $path=join("/",@pathtemp);
}

use lib $path."/fasta/";


use FastaParser;



my @hashDeflinepattern2Value=("chr",
			      "chrun.",
			      "scaffold_",
			      "contig",
			      "ultra",
			      "zv._scaffold",
			      "zv._na");



sub trimAlpha{
  my ($stringToTrim)=@_;
  my $stringToReturn="";

  for(my $index=0;$index<length($stringToTrim);$index++){
    my $char=substr($stringToTrim,$index,1);
    if($char =~ /[\d\.]/){
      if($char =~ /\d/){
	$stringToReturn.=$char;
      }elsif($char =~ /\./){
	#do nothing
      }else{
      	die "Problem#1203\n";
      }
    }else{
      if(length($stringToReturn) != 0){
	return $stringToReturn;
      }
    }
  }
  return $stringToReturn;
}





# a is lesser  than b return -1
# a is greater than b return 1
sub compareChr{
  my $aName=lc($a);
  my $bName=lc($b);
  my $aIndex=-1;
  my $bIndex=-1;
  my $aRest;
  my $bRest;

  #  print "a = $aName b $bName\n";

  if($aName eq $bName){
    die "Duplicate chromsome names (Name = $aName found twice)\n";
  }

  for(my $index=0;$index<=$#hashDeflinepattern2Value;$index++){
    my $pattern="^".$hashDeflinepattern2Value[$index]."(.+)\$";
    if($aName =~ /$pattern/){
      $aIndex = $index;
      $aRest = $1;
    }
  }


  for(my $index=0;$index<=$#hashDeflinepattern2Value;$index++){
    my $pattern="^".$hashDeflinepattern2Value[$index]."(.+)\$";
    if($bName =~ /$pattern/){
      $bIndex = $index;
      $bRest = $1;
    }
  }


  if($aIndex == -1 && $bIndex == -1 ){
    return ($aName cmp $bName);
  }


  if( $aIndex < $bIndex ){
    return -1;
  }elsif( $aIndex > $bIndex ){
    return 1;
  }else{
    my $aNameNum=trimAlpha($aName);
    my $bNameNum=trimAlpha($bName);

    if(length($aNameNum) == 0 && length($bNameNum) == 0 ){
      return ($aName cmp $bName);
    }

    #a<b
    if(length($aNameNum) != 0 && length($bNameNum) == 0 ){
      return -1;
    }

    #a>b
    if(length($aNameNum) == 0 && length($bNameNum) != 0 ){
      return 1;
    }

    #a>b
    if(length($aNameNum) != 0 && length($bNameNum) != 0 ){
      if($aNameNum eq $bNameNum){
	return ($aName cmp $bName);
      }else{
	return ($aNameNum <=> $bNameNum);
      }
    }

  }

  die "createChrosomeIndex.pl : Wrong Statement reached\n";
}

sub sortChromosome{
  my (@listOfChr)=@_;
  return (sort compareChr @listOfChr);
}

my $usage="Please enter the directory where the chromsomes are located\nUsage:\n$0 [directory where the fasta files for each chromosome are]\n";

if(! $ARGV[0]){
  die $usage;
}else{
  my $directory=$ARGV[0];
  if(!(-d $directory)){
    print "Directory : ".$directory." does not exists\n";
    die $usage;
  }else{
    if(substr($directory,-1) ne "/"){
      $directory.="/";
    }

    #    my $command="/bin/ls -1 $directory*fa";
    my $command="/usr/bin/find $directory |grep \"fa\$\"";
    my @result=split("\n",`$command`);
    #my @arrayOfDeflines;
    my @arrayOfFiles;

    my $hashFile2Chr;



    foreach my $file (@result){
      warn "Reading file $file\n";
      my $filename;

      if($file =~ /^$directory(.+)$/){
	$filename=$1;
	#push(@{$arrayOfFiles},$filename);
	push(@arrayOfFiles,$filename);
      }else{
	die "Wrong file pattern : ".$file."\n";
      }

      #print "file $file\n";
      my $fp = FastaParser->new(filename => $file);
      while($fp->hasRecord()){
	my $deflineToAdd=$fp->getNextRecord()->getDefline();
	if($deflineToAdd =~ />(.*)/){
	  $deflineToAdd=$1;
	}else{
	  die "Wrong defline pattern";
	}

	if(exists $hashFile2Chr->{$filename}){
	  push(@{$hashFile2Chr->{$filename}}, $deflineToAdd);
	}else{
	  $hashFile2Chr->{$filename} = [$deflineToAdd];
	}

      }

    }

    my @arrayOfFilesSort = sort compareChr  @arrayOfFiles;

    foreach my $file (@arrayOfFilesSort) {
      my @arrayOfChrTemp1=@{$hashFile2Chr->{$file}};
      my @arrayOfChrTemp2 = sort compareChr  @arrayOfChrTemp1;
      $hashFile2Chr->{$file} = [@arrayOfChrTemp2];
    }



    foreach my $file (@arrayOfFilesSort) {
      print $file."\n";
      foreach my $chr (@{$hashFile2Chr->{$file}}) {
	print ">".$chr."\n";
      }
    }


  }


}
