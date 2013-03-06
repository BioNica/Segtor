#!/usr/bin/perl



=head1 NAME

   Segtor::main.pl


=head1 Disclaimer

This software and all the data therein were developed by the Bioinformatics and Computational Biology Lab (Laboratorio de Bioinformatica e Biologia Computational (LBBC)) of the Brazilian National Cancer Institute (Instituto Nacional de Cancer (INCA)) in Rio de Janeiro, RJ, Brazil. The INCA is part of the Health Ministry (Ministerio da Saude) of the Federal Government of Brazil. This software was developed by the staff of the LBBC as part of their official duties. According to federal regulations, this software is therefore public domain.

This software is provided as is and the INCA does make any guarantee regarding the accuracy/quality/reliability of the results. The INCA cannot be held accountable for any hardware or software damage that might be caused through the use or misuse of this software.

The LBBC encourages researchers using our software to acknowledge us. The software can be distributed freely and modified given that any derivative bears an acknowledgment that it was derived from our original software.

=head1 Author

Gabriel Renaud grenaud@inca.gov.br


=cut

local $| = 1;

our $VERSION = '1.3';


use strict;
use warnings;

#Perl libraries
use Cwd         qw(abs_path);
use List::Util  qw(min max);
use Time::HiRes qw(time);
use POSIX       qw(ceil floor);
use Storable    qw(store retrieve freeze thaw dclone);

use Getopt::Std;
use File::stat;

use Data::Dumper;

my $path;

BEGIN {
  my @pathtemp = split("/", abs_path($0) );
  delete($pathtemp[$#pathtemp]);
  delete($pathtemp[$#pathtemp]);
  $path=join("/",@pathtemp);
}




use lib $path."/SegmentTreeArray/";
use SegmentTreeRange;
use SegmentTree;
use SegmentTree;

use lib $path."/RetrieveDataFromUCSC/";
use RetrieveDataUCSC;

use lib $path."/GeneDataStructure/";
use CreateRefSeqGenes;
use CreateKnownGenes;
use CreateEnsemblGenes;
use CreateEsts;

use lib $path."/ChromosomeIndexing/";
use ChromosomeIndexing;

use lib $path."/Segtor/";
use SNPstructure;
use StatisticsModule;


my $createChrIndexPath = $path."/CreateChromosomeIndex/";
my $splitFastaPath     = $path."/SplitFasta/";

my $dataDirectory;

my @speciesWithSNPData=("hg17","hg18","hg19");



my $colorBest   = "0xc6d2ff";
my $colorSingle = "0x92a9ff";
my $colorNormal = "0x5b7eff";

my $directoryInHome=".segtor";
my $genomeDirectory="chromosomes";
my $databaseDirectory="database";
my $structureDirectory="structures";

my $numberOfSecondsToPromptUpdate="15552000"; # 6months x 30 days x 24 hours x 60 minutes x 60 seconds
my $genomebinsize="100000";

my $availableDatabases={"refseq"      => ["refGene","refLink","refFlat","refSeqAli"],
			"ensembl"     => ["ensGene"],
			"knowngene"   => ["knownGene"],
			"ests"        => ["all_est"]  };

my $singleTranscriptDatabase=1;
my $MAXBPFORINTERVALS=2000;
my $MAXBPFORRANGE=100000000;
my $arrayOfTrees=[];
my $arrayOfTreesSNP=[];
my $arrayOfTSS=[];
my $arrayOfSNPStructures=[];
my $suffixForStructure;
my $promptBeforeDownload=1;

my $onlyBuildStruct=0;

my $useExternalFile=0;
my $listOfExternalFilesPath;


my $genomeIsNeededForBuild=1;
my $doNotPrintGenesOut=1;
my $doNotPrintXML=0;
my $printFullFiles=1;
my $printStatsOnly=0;

my $noRefBP=0;

my @arrayOfDatabases;

my @arrayOfRanges;
my $defaultRange=0;


my $databaseFiles=[];


my $usage="$0 -s [species code] -m [mode] -f [file] -d [database]\n".

  "\nMandatory parameters\n".
  "\t-s [species code] Code for species used by UCSC\n".

  "\n\tFor the annotation mode, please specify either:\n".
  "\t\t-m 1 Annotate coordinate mode\n".
  "\t\t-m 2 Annotate interval mode\n".
  "\t\t-m 3 Annotate SNPs\n".
  "\t\t-m 4 Find closest 5' end\n".
  "\t\t-m 5 Insertion/Deletion/Translocation \n".

  "\n\tFor the input, please specify either:\n".
  "\t\t-f [file] File to annotate\n".
  "\t\tor\n".
  "\t\t-l [file of list of files] File the full path of the files to annotate\n".

  "\n\tFor the database, please specify either:\n".
  "\t\t-d [ests]           The annotation will be based on ESTs (mode 1,2 & 4 only)\n".
  "\t\t-d [knowngene]      The annotation will be based on UCSC Known Genes (mode 1,2 & 4 only)\n".
  "\t\t-d [refseq]         The annotation will be based on refSeqs\n".
  "\t\t-d [ensembl]        The annotation will be based on Ensembl Genes\n".
  "\t\t-d [other database] Another database that has been constructed using the -e for mode 1,2\n".


  "\n\nParameters with default values\n".
  "\t\t-c [directory] Use this directory to store indices and download the chromosomes/database files\n".
  "\t\t               (default: \$HOME/".$directoryInHome."/\n".
  "\t\t-r [range] Range(s) (for -m 1,2,3 only) (default: $defaultRange) in base pairs to use to annotate coordinates\n".
  "\t\t           (ex: to use 1kb and 5kb : -r \"1000,5000\")\n".


  "\n\nOptional parameters\n".
  "\tAnnotation options\n".
  "\t\t-b Check dbSNP for existing snps (for -m 3 only, requires a large RAM)\n".
  "\t\t-a Produce amino acid sequence (for -m 3 only, triggered by default for -m 5)\n".
  "\t\t-i Do not check the reference base pair\n".
  "\t\t-n Report inputs in terms of genes (genes.out file) \n".
  "\t\t-g n:m for mode 4, report n bins of m base pairs each\n".
  "\t\t       (ex: -g 10:1000 will create 10 bins of 1kb each side of a tss)\n".
  "\t\t-x do not produce XML output files, used to save space\n".
  "\t\t-t produce only a part of the total output (all.out), triggers -x\n".
  "\t\t-z only produce stats\n".



  "\n\tIndex building options\n".
  "\t\t-o Build the index files for the given species/database/mode (no annotation)\n".
  "\t\t\t-e [files] When using -o, specify a comma-separated list of files in psl or bed (with strand) format \n".
  "\t\t\t            (from UCSC Table Browser for example) to create the segment tree\n".
  "\t\t-b Download dbSNP and build index for it (for -m 3 only)\n".
  "\t\t-y Do not prompt prior to downloading the chromosomes/database files\n".
  "\t\t-p [proxy server] FTP proxy server to use (ex: http://127.0.0.1:21)\n".

  "\n\n";



my $ftpTarget  = "ftp://hgdownload.cse.ucsc.edu";


my $speciesCode;
my $databaseName;
my $ranges;
my $annotationMode;
my $customDatabase;

my $proxyServer;
my $retrieveDataUCSC;
my $produceAASeq=0;
my $compareWithDBSNP=0;

my $dbSNPFile;

my $numberBins=10;
my $bpInBins  =1000;

my @arrayOfFilesToAnnotate;








































################################################
################################################
##                                            ##
##                                            ##
##           BEGIN SUBROUTINES                ##
##                                            ##
##                                            ##
################################################
################################################


=item terminate

   Kills the current program with the
   kill -[code]
   with code being given as argument

=cut
sub terminate{
  my ($code)=@_;
  print "Terminating program\n";
  kill($code, $$);
}







=item selectSingleBestTranscriptSNP

    Subroutine to select the best result among an array
    annotations for a SNP.

=cut
sub selectSingleBestTranscriptSNP{
  my (@arrayOfResults)=@_;
  my $minScore=10;
  my $coding=0;
  my $nonSyn=0;
  my $reliable=0;
  my $lengthTranscript=-1;
  my $singleBestResult;
  my $spliceForBest=0;


  foreach my $resultHash (@arrayOfResults){

    if($minScore > $resultHash->{'annoResult'}->{'code'} ){ #level 1, if one with a better code was found take it
	$singleBestResult = $resultHash;
	$lengthTranscript = $resultHash->{'transcript'}->getTranscriptionRange()->getLength();
	$minScore         = $resultHash->{'annoResult'}->{'code'} ;

	if( ($resultHash->{'annoResult'}->{'code'} == 2) ){ #new record is in intron
	  $spliceForBest=$resultHash->{'snpResult'}->{'splice'};
	}

	if( ($resultHash->{'annoResult'}->{'code'} == 1) && #new record is in exon
	    ($resultHash->{'snpResult'}->{'coding'}) ){     #and in a coding region

	  if($resultHash->{'snpResult'}->{'reliable'}){ #the new one is reliable
	    my $resultNonSyn;
	    if ($resultHash->{'snpResult'}->{'aaRef'} eq $resultHash->{'snpResult'}->{'aaRead'}) {
	      $resultNonSyn   = 0;
	    } else {
	      $resultNonSyn   = 1;
	    }
	    $reliable = $resultHash->{'snpResult'}->{'reliable'};
	    $coding   = $resultHash->{'snpResult'}->{'coding'};
	    $nonSyn   = $resultNonSyn;
	  }else{      #the new one is reliable
	    $reliable = $resultHash->{'snpResult'}->{'reliable'};
	    $coding   = $resultHash->{'snpResult'}->{'coding'};
	  }
	}

    }elsif($minScore == $resultHash->{'annoResult'}->{'code'} ){ #level 1, tie in code, look to see if we can find a coding one

      if($resultHash->{'annoResult'}->{'code'} != 1){ #level2, not an exon, just pick the longest

	if($resultHash->{'annoResult'}->{'code'} == 2){ #both the best one and the new one are in intron

	  if( ($spliceForBest == 0 ) &&                        # the case where the new one causes a splice disruption
	      ($resultHash->{'snpResult'}->{'splice'} == 1)){ # but not the current best one
	    #                                                    the best one becomes the current one
	    $singleBestResult = $resultHash;
	    $lengthTranscript = $resultHash->{'transcript'}->getTranscriptionRange()->getLength();
	    $minScore         = $resultHash->{'annoResult'}->{'code'} ;
	    $spliceForBest    = $resultHash->{'snpResult'}->{'splice'};
	  }elsif( ($spliceForBest == 1 ) &&                   # where the new one does not cause a splice disruption
		  ($resultHash->{'snpResult'}->{'splice'} == 0) ){ # but the current best one does
	    #do nothing, we already have the best
	  }else{ #we have a tie in splice disruption, we pick the longest
	    if ($lengthTranscript < $resultHash->{'transcript'}->getTranscriptionRange()->getLength() ) {
	      $singleBestResult = $resultHash;
	      $lengthTranscript = $resultHash->{'transcript'}->getTranscriptionRange()->getLength();
	      $minScore         = $resultHash->{'annoResult'}->{'code'} ;
	      $spliceForBest    = $resultHash->{'snpResult'}->{'splice'};
	    }
	  }

	}else{  #not in intron or exon, we pick the longest
	  if ($lengthTranscript < $resultHash->{'transcript'}->getTranscriptionRange()->getLength() ) {
	    $singleBestResult = $resultHash;
	    $lengthTranscript = $resultHash->{'transcript'}->getTranscriptionRange()->getLength();
	    $minScore         = $resultHash->{'annoResult'}->{'code'} ;
	  }
	}

      }else{ #level2, in an exon

	if(    $coding == 0 && $resultHash->{'snpResult'}->{'coding'} == 0){ #level3, both are non coding

	  if(     $reliable==0 && $resultHash->{'snpResult'}->{'reliable'} == 0){ #level4, both unreliable

	    if ($lengthTranscript < $resultHash->{'transcript'}->getTranscriptionRange()->getLength() ) {
	      $singleBestResult = $resultHash;
	      $lengthTranscript = $resultHash->{'transcript'}->getTranscriptionRange()->getLength();
	      $minScore         = $resultHash->{'annoResult'}->{'code'} ;
	    }

	  }elsif( $reliable==0 && $resultHash->{'snpResult'}->{'reliable'} == 1){ #level4, new one is reliable, take it
	    $reliable = $resultHash->{'snpResult'}->{'reliable'};
	    $singleBestResult = $resultHash;
	    $lengthTranscript = $resultHash->{'transcript'}->getTranscriptionRange()->getLength();
	    $minScore         = $resultHash->{'annoResult'}->{'code'} ;
	  }elsif( $reliable==1 && $resultHash->{'snpResult'}->{'reliable'} == 0){ #level4, new one is unreliable, ignore it
	    #do nothing
	  }elsif( $reliable==1 && $resultHash->{'snpResult'}->{'reliable'} == 1){ #level4, both reliable

	    if ($lengthTranscript < $resultHash->{'transcript'}->getTranscriptionRange()->getLength() ) {
	      $singleBestResult = $resultHash;
	      $lengthTranscript = $resultHash->{'transcript'}->getTranscriptionRange()->getLength();
	      $minScore         = $resultHash->{'annoResult'}->{'code'} ;
	    }

	  }




	}elsif($coding == 0 && $resultHash->{'snpResult'}->{'coding'} == 1){ #level3, the new one is coding, take it


	  if(    $resultHash->{'snpResult'}->{'reliable'} == 0){
	    $singleBestResult = $resultHash;
	    $lengthTranscript = $resultHash->{'transcript'}->getTranscriptionRange()->getLength();
	    $minScore         = $resultHash->{'annoResult'}->{'code'} ;

	  } elsif ($resultHash->{'snpResult'}->{'reliable'} == 1) {
	    my $resultNonSyn;
	    if ($resultHash->{'snpResult'}->{'aaRef'} eq $resultHash->{'snpResult'}->{'aaRead'}) {
	      $resultNonSyn   = 0;
	    } else {
	      $resultNonSyn   = 1;
	    }

	    $reliable = $resultHash->{'snpResult'}->{'reliable'};
	    $coding   = $resultHash->{'snpResult'}->{'coding'};
	    $nonSyn   = $resultNonSyn;
	    $singleBestResult = $resultHash;
	    $lengthTranscript = $resultHash->{'transcript'}->getTranscriptionRange()->getLength();
	    $minScore         = $resultHash->{'annoResult'}->{'code'} ;

	  } else {
	    print STDERR "Invalid state in selectSingleBestTranscriptSNP()\n";
	    terminate(9);
	  }

	}elsif($coding == 1 && $resultHash->{'snpResult'}->{'coding'} == 0){ #level3, the new one is non-coding
	  #do nothing
	}elsif ($coding == 1 && $resultHash->{'snpResult'}->{'coding'} == 1) { #level3, both are coding



	  if(     $reliable==0 && $resultHash->{'snpResult'}->{'reliable'} == 0){ #level4, both unreliable

	    if ($lengthTranscript < $resultHash->{'transcript'}->getTranscriptionRange()->getLength() ) {
	      $singleBestResult = $resultHash;
	      $lengthTranscript = $resultHash->{'transcript'}->getTranscriptionRange()->getLength();
	      $minScore         = $resultHash->{'annoResult'}->{'code'} ;
	    }

	  }elsif( $reliable==0 && $resultHash->{'snpResult'}->{'reliable'} == 1){ #level4, new one is reliable, take it
	    $reliable = $resultHash->{'snpResult'}->{'reliable'};
	    $singleBestResult = $resultHash;
	    $lengthTranscript = $resultHash->{'transcript'}->getTranscriptionRange()->getLength();
	    $minScore         = $resultHash->{'annoResult'}->{'code'} ;
	  }elsif( $reliable==1 && $resultHash->{'snpResult'}->{'reliable'} == 0){ #level4, new one is unreliable, ignore it
	    #do nothing
	  }elsif( $reliable==1 && $resultHash->{'snpResult'}->{'reliable'} == 1){ #level4, both reliable

	    my $resultNonSyn;
	    if ($resultHash->{'snpResult'}->{'aaRef'} eq $resultHash->{'snpResult'}->{'aaRead'}) {
	      $resultNonSyn   = 0;
	    } else {
	      $resultNonSyn   = 1;
	    }

	    if(    $nonSyn == 0 && $resultNonSyn == 0){ #level5, both non-syn
	      if ($lengthTranscript < $resultHash->{'transcript'}->getTranscriptionRange()->getLength() ) {
	        $reliable = $resultHash->{'snpResult'}->{'reliable'};
		$coding   = $resultHash->{'snpResult'}->{'coding'};
		$nonSyn   = $resultNonSyn;
		$singleBestResult = $resultHash;
		$lengthTranscript = $resultHash->{'transcript'}->getTranscriptionRange()->getLength();
		$minScore         = $resultHash->{'annoResult'}->{'code'} ;
	      }
	    }elsif($nonSyn == 0 && $resultNonSyn == 1){ #level5, new one if syn, take it
	      $reliable = $resultHash->{'snpResult'}->{'reliable'};
	      $coding   = $resultHash->{'snpResult'}->{'coding'};
	      $nonSyn   = $resultNonSyn;
	      $singleBestResult = $resultHash;
	      $lengthTranscript = $resultHash->{'transcript'}->getTranscriptionRange()->getLength();
	      $minScore         = $resultHash->{'annoResult'}->{'code'} ;
	    }elsif($nonSyn == 1 && $resultNonSyn == 0){ #level5, new one if non-syn
	      #do nothing
	    }elsif($nonSyn == 1 && $resultNonSyn == 1){ #level5, both are syn
	      if ($lengthTranscript < $resultHash->{'transcript'}->getTranscriptionRange()->getLength() ) {
	        $reliable = $resultHash->{'snpResult'}->{'reliable'};
		$coding   = $resultHash->{'snpResult'}->{'coding'};
		$nonSyn   = $resultNonSyn;
		$singleBestResult = $resultHash;
		$lengthTranscript = $resultHash->{'transcript'}->getTranscriptionRange()->getLength();
		$minScore         = $resultHash->{'annoResult'}->{'code'} ;
	      }
	    }else{ #end level5
	      print STDERR "Invalid state in selectSingleBestTranscriptSNP()\n";
	      terminate(9);

	    }
	  } #end level4, reliable



	} else { #end level3, coding
	  print STDERR "Invalid state in selectSingleBestTranscriptSNP()\n";
	  terminate(9);
	}

      }#end level2, in an exon


    }else{ #level 1, new one has a worse code
      #do nothing
    } #end level 1, new one has a worse code

  }


  return $singleBestResult;
}


=item printXML

    Subroutine to print the results as
    XML rather than the standard .out

=cut
sub printXML{
  my ($snp,$annotationMode,$INPUTALLXML,$coordInput,$chrInput,$idInput,$bestResult,$snpString,$singleGenes,$hashOfGenes,$arrayOfFieldHeader,$singleTranscriptDatabase)=@_;

  print $INPUTALLXML "<INPUTID name=\"".$idInput." $chrInput $coordInput\" chr=\"".$chrInput."\" coord=\"".$coordInput."\">"."\n";

  my $arrayXML=[];
  if($singleTranscriptDatabase){
    foreach my $geneCounter (keys %{$hashOfGenes}) {
	my $arrayXMLTemp=[];
	foreach my $record ( @{$hashOfGenes->{$geneCounter}} ) {
	  push(@{$arrayXMLTemp},$record);
	}
	push(@{$arrayXML},$arrayXMLTemp);
    }

    $colorSingle = $colorNormal;
    $colorBest   = $colorNormal;

  } else {
    foreach my $geneCounter (keys %{$hashOfGenes}) {

      if ($bestResult->{'transcript'}->getGene()->getGeneCounter() == $geneCounter) {

	my $arrayXMLTemp=[];
	foreach my $record ( @{$hashOfGenes->{$geneCounter}} ) {
	  if ($record->{'transcript'}->getId() eq $bestResult->{'transcript'}->getId()) {
	    unshift(@{$arrayXMLTemp},$record);
	  } else {
	    push(@{$arrayXMLTemp},$record);
	  }
	}

	unshift(@{$arrayXML},$arrayXMLTemp);
      } else {
	my $arrayXMLTemp=[];
	foreach my $record ( @{$hashOfGenes->{$geneCounter}} ) {
	  if ( $record->{'transcript'}->getId() eq $singleGenes->{$geneCounter} ) {
	    unshift(@{$arrayXMLTemp},$record);
	  } else {
	    push(@{$arrayXMLTemp},$record);
	  }
	}

	push(@{$arrayXML},$arrayXMLTemp);
      }
    }
  }


  my $bestGeneRecord=1;
  foreach my $geneRECORD (@{$arrayXML}) {
    if(!$singleTranscriptDatabase){
      print $INPUTALLXML "<GENE name=\"".$geneRECORD->[0]->{'transcript'}->getGeneSymbol()."\" link=\"NCBI\">"."\n";
    }
    my $singleGeneRecord=1;

    foreach my $transRECORD (@{$geneRECORD}) {
      print $INPUTALLXML "<TRANS name=\"".$transRECORD->{'transcript'}->getId()."\" ";






      if($snp){
	#                         0        1      2        3          4       5          6        7          8              9        10        11         12          13      14       15
	# @arrayOfFieldHeader=("#INPUTID","CHR","COORD","TRANSCRIPT","GENE","POSITION","INDEX","DISTANCE5p","DISTANCE3p","PARTIAL","CODONREF","CODONREAD","COORDSNP","AAREF","AAREAD","DBSNP");
	my @arraytemp1=split("\t",formatAnnotation(1, $transRECORD->{'annoResult'} ));
	my @arraytemp2=split("\t",formatSNPAnno($transRECORD->{'snpResult'}));

	my $indexStart1=5;
	my $indexStart2=10;
	for (my $indexFIELD=$indexStart1;$indexFIELD<$indexStart2;$indexFIELD++) {
	  print $INPUTALLXML "".lc($arrayOfFieldHeader->[$indexFIELD])."=\"". $arraytemp1[$indexFIELD-$indexStart1] ."\" ";
	}


	if($#arraytemp2 == -1){
	  for (my $indexFIELD=$indexStart2;$indexFIELD<=$#{$arrayOfFieldHeader};$indexFIELD++) {
	    print $INPUTALLXML "".lc($arrayOfFieldHeader->[$indexFIELD])."=\" \" ";
	  }

	  print $INPUTALLXML "comment=\"  \" ";

	}elsif($#arraytemp2 == 0){

	  for (my $indexFIELD=$indexStart2;$indexFIELD<=$#{$arrayOfFieldHeader};$indexFIELD++) {
	    print $INPUTALLXML "".lc($arrayOfFieldHeader->[$indexFIELD])."=\" \" ";
	  }

	  print $INPUTALLXML "comment=\"". $arraytemp2[0] ."\" ";
	}else{

	  if($snpString ne "-1"){
	   push(@arraytemp2,$snpString);
	  }

	  for (my $indexFIELD=$indexStart2;$indexFIELD<=$#{$arrayOfFieldHeader};$indexFIELD++) {
	    print $INPUTALLXML "".lc($arrayOfFieldHeader->[$indexFIELD])."=\"". $arraytemp2[$indexFIELD-$indexStart2] ."\" ";
	  }

	  print $INPUTALLXML "comment=\"  \" ";
	}

      }else{ #not snp
	#                         0        1      2        3          4       5          6        7          8              9
	#@arrayOfFieldHeader=("#INPUTID","CHR","COORD","TRANSCRIPT","GENE","POSITION","INDEX","DISTANCE5p","DISTANCE3p","PARTIAL");
	my @arraytemp;
	my $indexStart;
	if($annotationMode == 1){
	  @arraytemp=split("\t",formatAnnotation(1, $transRECORD->{'annoResult'} ));
	  $indexStart=5;
	}else{
	  @arraytemp=split("\t",formatIntervalAnnotation($transRECORD->{'annoResult'} ));
	  $indexStart=6;
	}


	for (my $indexFIELD=$indexStart;$indexFIELD<=$#{$arrayOfFieldHeader};$indexFIELD++) {
	  print $INPUTALLXML "".lc($arrayOfFieldHeader->[$indexFIELD])."=\"". $arraytemp[$indexFIELD-$indexStart] ."\" ";
	}

      }


      if ($bestGeneRecord) {
	print $INPUTALLXML "color=\"". $colorBest ."\" ";
	$bestGeneRecord=0;
	$singleGeneRecord=0;
      } elsif ($singleGeneRecord) {
	print $INPUTALLXML "color=\"". $colorSingle ."\" ";
	$singleGeneRecord=0;
      } else {
	print $INPUTALLXML "color=\"". $colorNormal ."\" ";
      }

      print $INPUTALLXML "link=\"NCBI\"/>"."\n";
    }
    if(!$singleTranscriptDatabase){
      print $INPUTALLXML "</GENE>"."\n";
    }
  }
  print $INPUTALLXML "</INPUTID>"."\n";



}



=item selectSingleBestTranscript

    Subroutine to select and return
    the best result among an array of results
    for mode 1 and 2

=cut
sub selectSingleBestTranscript{
  my (@arrayOfResults)=@_;
  my $minScore=10;
  my $lengthTranscript=-1;
  my $singleBestResult;

  foreach my $resultHash (@arrayOfResults){

    if($minScore > $resultHash->{'annoResult'}->{'code'} ){
      $singleBestResult = $resultHash;
      $lengthTranscript = -1;
      $minScore         = $resultHash->{'annoResult'}->{'code'} ;
    }elsif($minScore == $resultHash->{'annoResult'}->{'code'} ){
      if($lengthTranscript < $resultHash->{'transcript'}->getTranscriptionRange()->getLength() ){
	$singleBestResult = $resultHash;
	$lengthTranscript = $resultHash->{'transcript'}->getTranscriptionRange()->getLength() ;
      }
    }
  }

  return $singleBestResult;
}


=item findNewestSNPFile

    Subroutine to select the most recent dbSNP
    file on the hard disk.

=cut
sub findNewestSNPFile{
  my @filesInDatabase = (@_);
  my $mostRecentFile = -1;
  my $fileForMostRecentFile = -1;
  my @arrayOfOldSNPFiles;
  my $fileSuffix="";

  foreach my $dataFile (@filesInDatabase){
    if($dataFile =~ /^snp(\d+)(\.txt(.gz)?)?$/){
      if($1 > $mostRecentFile){
	if($mostRecentFile != -1){
	  push(@arrayOfOldSNPFiles,$fileForMostRecentFile);
	}
	$fileForMostRecentFile=$dataFile;
	$mostRecentFile=$1;
	if($2){
	  if($3){
	    $fileSuffix=$2.$3;
	  }else{
	    $fileSuffix=$2;
	  }
	}
      }else{
	push(@arrayOfOldSNPFiles,$dataFile);
      }
    }
  }

  if($mostRecentFile == -1){
    return (-1,-1);
  }else{
    return ("snp".$mostRecentFile.$fileSuffix,@arrayOfOldSNPFiles);
  }

}


=item eliminateDuplicateTranscripts

    Subroutine which received an reference to array of transcripts
    and returns a reference to an new array of transcript but
    without transcript duplicate (using transcript counter)

=cut
sub eliminateDuplicateTranscripts{
  my ($initialArrayTranscripts)=@_;
  my $arrayTranscriptsToReturn;
  my %hashIDtrans;

  foreach my $trans (@{$initialArrayTranscripts}){
    if(exists $hashIDtrans{$trans->{'transcript'}->getTranscriptCounter()}){
    }else{
      $hashIDtrans{$trans->{'transcript'}->getTranscriptCounter()}=$trans;
    }
  }

  $arrayTranscriptsToReturn=[(values %hashIDtrans)];

  return $arrayTranscriptsToReturn;
}


=item eliminateDuplicateRecords

    Subroutine which received an reference to array of records
    built from a custom database and returns a reference to an
    new array of transcript but without records duplicate

=cut
sub eliminateDuplicateRecords{
  my ($initialArrayRecords)=@_;
  my $arrayRecordsToReturn;
  my %hashIDrecord;

  foreach my $recordInArray (@{$initialArrayRecords}){
    if(exists $hashIDrecord{$recordInArray->{'record'}->{'id'}."#".$recordInArray->{'record'}->{'range'}->getStart()."#".$recordInArray->{'record'}->{'range'}->getEnd()."#".$recordInArray->{'record'}->{'strand'}}){
    }else{
      $hashIDrecord{$recordInArray->{'record'}->{'id'}."#".$recordInArray->{'record'}->{'range'}->getStart()."#".$recordInArray->{'record'}->{'range'}->getEnd()."#".$recordInArray->{'record'}->{'strand'}}=$recordInArray;
    }
  }

  $arrayRecordsToReturn=[(values %hashIDrecord)];

  return $arrayRecordsToReturn;
}



=item formatDBSNP

   Subroutine to format the output
   for a dbSNP annotation

=cut
sub formatDBSNP{
  my ($hashSNP)=@_;
  if($hashSNP->{'foundCoord'} == 0){ #nothing in dbsnp
    return "no records found in dbsnp";
  }else{
    if($hashSNP->{'foundSNPs'}->{'alreadyFoundInDBSNP'} == 0){ #new bp in dbsnp
      return "Existing SNPs : ".$hashSNP->{'foundSNPs'}->{'snpid'}." but different bp";
    }else{
      return "Existing SNPs : ".$hashSNP->{'foundSNPs'}->{'snpid'}." with same bp";
    }
  }
}

=item formatSNPAnno

   Subroutine to format the output
   line for a SNP annotation

=cut
sub formatSNPAnno{
  my ($hashAnnotation)=@_;

  if ($hashAnnotation->{'reliable'} == 1) {

    if ($hashAnnotation->{'code'} == 1 ){ #in exon
      if ($hashAnnotation->{'coding'} == 1 ) {

	if ($hashAnnotation->{'seqRef'}  &&  $hashAnnotation->{'seqRead'}) {
	  return $hashAnnotation->{'codonRef'}."\t".$hashAnnotation->{'codonRead'}."\t".$hashAnnotation->{'aaRef'}."\t".$hashAnnotation->{'aaRead'}."\t".$hashAnnotation->{'coordInAASeq'}."\t".$hashAnnotation->{'seqRef'}."\t".$hashAnnotation->{'seqRead'};
	} else {
	  return $hashAnnotation->{'codonRef'}."\t".$hashAnnotation->{'codonRead'}."\t".$hashAnnotation->{'aaRef'}."\t".$hashAnnotation->{'aaRead'};
	}
      } else {
	if ($hashAnnotation->{'utr'} == 5) {
	  return "untranslated 5'UTR";		#untranslated
	}elsif ($hashAnnotation->{'utr'} == 3) {
	  return "untranslated 3'UTR";		#untranslated
	}else{
	  return "untranslated";		#untranslated
	}
      }
    } else {
      if (  $hashAnnotation->{'splice'} ) {
	return "potential splice disruption";		#untranslated
      }else{
	return " ";
      }
    }
  } else {
    if ($hashAnnotation->{'code'} == 1 ){
      return "unreliable frames";
    }
  }
  return " ";
}

=item formatInsertAnno

   Subroutine to format the output
   line for an insert annotation

=cut
sub formatInsertAnno{
  my ($hashAnnotation,$hashAnnotationCoord)=@_;

  if ($hashAnnotation->{'reliable'} == 1) {

    if ($hashAnnotation->{'code'} == 1 ){
      if ($hashAnnotation->{'coding'} == 1 ) {
	if($hashAnnotationCoord->{'code'} != 2 ){
	  return $hashAnnotation->{'coordInAASeq'}."\t".$hashAnnotation->{'seqRef'}."\t".$hashAnnotation->{'seqRead'};
	}else{
	  return "intronic"; #within intron
	}
      } else {
	return "untranslated";		#untranslated
      }
    } else {
      return " ";		#non-exon
    }
  } else {
    if ($hashAnnotation->{'code'} == 1 ){
      return "unreliable frames";
    }
  }
  return " ";
}


=item formatDeletionAnno

   Subroutine to format the output
   line for a deletion annotation

=cut
sub formatDeletionAnno{
  my ($hashAnnotation,$hashAnnotationInterval)=@_;

  if ($hashAnnotation->{'reliable'} == 1) {
    if ($hashAnnotation->{'coding'} == 1 ) {
      if($hashAnnotationInterval->{'code'} != 3 ){
	return $hashAnnotation->{'seqRef'}."\t".$hashAnnotation->{'seqRead'};
      }else{
	return "intronic"; #within intron
      }
    } else {
      return "untranslated";	#untranslated
    }
  } else {
    if ($hashAnnotation->{'coding'} == 1 ) {
      return "unreliable frames";
    }
  }
  return " ";
}



=item formatTransAnno

   Subroutine to format the output
   line for a translocation annotation

=cut
sub formatTransAnno{
  my ($hashAnnotation)=@_;

  if ($hashAnnotation->{'reliable'} == 1) {

    if ($hashAnnotation->{'coding'} == 1 ) {
      return $hashAnnotation->{'seqRead'};
    } else {
      return "untranslated";	#untranslated
    }
  } else {
    if ($hashAnnotation->{'coding'} == 1 ) {
      return "unreliable frames";
    }
  }
  return " ";
}



=item formatAnnotation

   Subroutine to format the output
   line for a coordinate annotation

=cut
sub formatAnnotation{
  my ($putIndex,$hashAnnotation)=@_;

  my $code       =$hashAnnotation->{'code'};
  my $distance5P =$hashAnnotation->{'distance5P'};
  my $distance3P =$hashAnnotation->{'distance3P'};
  my $partial ;
  my $index      =-1;
  if(exists $hashAnnotation->{"index"} ){
    $index      = $hashAnnotation->{"index"};
  }


  if (     $code  == 1) {
    $code="EXN";
    $partial= "###";
  } elsif ( $code  == 2) {
    $code="INT";
    $partial= "###";
  } elsif ( $code  == 3) {
    $code="UPS";
    $partial= $hashAnnotation->{"partial"};
  } elsif ( $code  == 4) {
    $code="DWS";
    $partial= $hashAnnotation->{"partial"};
  } else {
    print STDERR "Problem in formatAnnotation, wrong code\n";
    terminate(9);
  }

  if ($putIndex) {
    return "".$code."\t".$index."\t".$distance5P."\t".$distance3P."\t".$partial;
  } else {
    return "".$code."\t".$distance5P."\t".$distance3P."\t".$partial;
  }

}


=item formatIntervalAnnotation

   Subroutine to format the output
   line for a interval annotation

=cut
sub formatIntervalAnnotation{
  my ($hashAnnotation)=@_;

  my $code       =$hashAnnotation->{'code'};
  my $indicesEx;
  my $indicesIn;
  my $partial ;
  my $stringCoding="";

  if ( $hashAnnotation->{'coding'} ) {
    if (      $hashAnnotation->{'coding'} == 1 ) {
      $stringCoding="not in coding range";
    } elsif ( $hashAnnotation->{'coding'} == 2 ) {
      $stringCoding="encompassed in coding range";
    } elsif ( $hashAnnotation->{'coding'} == 3 ) {
      $stringCoding="overlaps a coding start/end";
    } else {
      print STDERR "Wrong code for coding annotation\n";
      terminate(9);
    }

  }


  if (          $code  == 1     ) {
    $code="contained in transcript and contained in a single exon ".$stringCoding;
    $indicesEx  =$hashAnnotation->{'indicesEx'};
  } elsif (     $code  == 2     ) {
    $code="contained in transcript and overlaps exonic and intronic regions ".$stringCoding;
    $indicesEx  =$hashAnnotation->{'indicesEx'};
    $indicesIn  =$hashAnnotation->{'indicesIntron'};
  } elsif (     $code  == 3     ) {
    $code="contained in transcript and contained in single intron ".$stringCoding;
    $indicesIn  =$hashAnnotation->{'indicesIntron'};
  } elsif (     $code  == 4     ) {
    $code="interval completely contains the transcript ".$stringCoding;
  } elsif (     $code  == 5     ) {
    $code="overlaps 5' end and first exon only ".$stringCoding;
    $partial = $hashAnnotation->{'partial'};
    $indicesEx  =$hashAnnotation->{'indicesEx'};
  } elsif (     $code  == 6     ) {
    $code="overlaps 5' end and exons, introns ".$stringCoding;
    $indicesEx  =$hashAnnotation->{'indicesEx'};
    $indicesIn  =$hashAnnotation->{'indicesIntron'};
    $partial    =$hashAnnotation->{'partial'};
  } elsif (     $code  == 7     ) {
    $code="overlaps 3' end and last exon only ".$stringCoding;
    $partial    =$hashAnnotation->{'partial'};
    $indicesEx  =$hashAnnotation->{'indicesEx'};
  } elsif (     $code  == 8     ) {
    $code="overlaps 3' end and exons, introns ".$stringCoding;
    $indicesEx  =$hashAnnotation->{'indicesEx'};
    $indicesIn  =$hashAnnotation->{'indicesIntron'};
    $partial    =$hashAnnotation->{'partial'};
  } elsif (     $code  == 9     ) {
    $code="upstream";
    $partial    =$hashAnnotation->{'partial'};
  } elsif (     $code  == 10     ) {
    $code="downstream";
    $partial    =$hashAnnotation->{'partial'};
  }else{
    print STDERR "Problem in formatIntervalAnnotation, wrong code\n";
    terminate(9);
  }

  my $stringToReturn="".$code;

  if ($indicesEx) {
    $stringToReturn.="\t".$indicesEx;
  }else{
    $stringToReturn.="\t"."-";
  }

  if ($indicesIn) {
    $stringToReturn.="\t".$indicesIn;
  }else{
    $stringToReturn.="\t"."-";
  }

  if ($partial) {
    $stringToReturn.="\t".$partial;
  }else{
    $stringToReturn.="\t"."###";
  }


  return $stringToReturn;
}



=item testInternetConnection

   Subroutine to check if the current
   machine can get access to the internet.
   First called by findClosestTSS()

=cut
sub testInternetConnection{

  while (1) {
    print "\ntesting internet connection...";
    my $response=$retrieveDataUCSC->testInternet();
    if ($response) {
      print "success!\n";
      last;
    } else {
      print "failed, enter your proxy ip (ex: http://127.0.0.1:21) :";
      while(1){
	$proxyServer=<>;

	if($proxyServer =~ /^http:\/\/\d+\.\d+\.\d+\.\d+:\d+\/?$/){
	  last;
	}else{
	  print "\nwrong format, please enter a proxy in the following format: http://127.0.0.1:21 :";
	}
      }
      chomp($proxyServer);
      $retrieveDataUCSC = RetrieveDataUCSC->new(ftpTarget  => $ftpTarget,
						httpProxy  => $proxyServer,
						ftpProxy   => $proxyServer,
						waisProxy  => $proxyServer);
    }
  }


}



=item returnClosestTSS

   Recursive subroutine to search the $arrayToSearch
   using the $numberToLookFor

=cut
sub returnClosestTSS{
  my ($arrayToSearch,$numberToLookFor,$index1,$index2)=@_;

  if( ($index1 + 1) == $index2){
    if(  abs($arrayToSearch->[$index1]->{'tss'} - $numberToLookFor) <
	 abs($arrayToSearch->[$index2]->{'tss'} - $numberToLookFor) ){
      return $arrayToSearch->[$index1];
    }else{
      return $arrayToSearch->[$index2];
    }
  }else{
    my $middleIndex=floor( ($index2+$index1)/2);

    if (      $arrayToSearch->[$middleIndex]->{'tss'} == $numberToLookFor) {
      return $arrayToSearch->[$middleIndex];
    } elsif ( $arrayToSearch->[$middleIndex]->{'tss'}  < $numberToLookFor) {
      returnClosestTSS($arrayToSearch,$numberToLookFor,$middleIndex,$index2)
    } elsif ( $arrayToSearch->[$middleIndex]->{'tss'}  > $numberToLookFor) {
      returnClosestTSS($arrayToSearch,$numberToLookFor,$index1,$middleIndex)
    } else {
      print STDERR "Major error\n";
      terminate(9);
    }
  }
}


=item returnClosestTSS

   Subroutine to call returnClosestTSS
   on $arrayRef and $genomicCoordinate

=cut
sub findClosestTSS{
  my ($arrayRef,$genomicCoordinate)=@_;

  if($#{$arrayRef} == 0){
    return $arrayRef->[0];
  }else{
    return returnClosestTSS($arrayRef,$genomicCoordinate,0,$#{$arrayRef});
  }

}



################################################
################################################
##                                            ##
##                                            ##
##             END SUBROUTINES                ##
##                                            ##
##                                            ##
################################################
################################################







































################################################
################################################
##                                            ##
##                                            ##
##           BEGIN MAIN PROGRAM               ##
##                                            ##
##                                            ##
################################################
################################################





############################################
##                                        ##
##                                        ##
##  BEGIN READING COMMAND LINE ARGUMENT   ##
##                                        ##
##                                        ##
############################################



my %opts;
getopts('abc:d:e:f:g:hil:m:nop:r:s:tuyxz', \%opts);


if( scalar(keys %opts) == 0 || (exists $opts{'h'}) ){
  print "usage:\n".$usage."";
  exit(0);
}


####### MANDATORY PARAMETERS #######
#species
if( !(exists $opts{'s'}) || !$opts{'s'}){
  print STDERR "The -s argument is mandatory\n"."usage:\n".$usage."";
  terminate(9);
}else{
  $speciesCode=$opts{'s'};
}



#mode
if( !(exists $opts{'m'}) || !$opts{'m'}){
  print STDERR "The -m argument is mandatory\n". "usage:\n".$usage."";
  terminate(9);
}else{
  if($opts{'m'} !~ /\d/){
    print STDERR "The -m argument must be numerical\n". "usage:\n".$usage."";
    terminate(9);
  }

  $annotationMode=$opts{'m'};
}

#database
if( !(exists $opts{'d'}) || !$opts{'d'}){
  print STDERR "The -d argument is mandatory\n". "usage:\n".$usage."";
  terminate(9);
}else{
  $databaseName=$opts{'d'};

  if( exists $availableDatabases->{$databaseName} ){
    @arrayOfDatabases=@{ $availableDatabases->{$databaseName} };
    $customDatabase=0;
  }else{
    if($opts{'d'} !~ /^\w+$/){
      print STDERR "Invalid database name: $databaseName\n"." the name must be only alphanumerical characters and underscores ";
      terminate(9);

    }
    $customDatabase=1;

    if ( $annotationMode == 1 || $annotationMode == 2 ) {
      #fine
    }else{
      print STDERR "The use of a custom database is only available for mode 1 and 2, otherwise please use one of the following databases ".join(" ",(keys %{$availableDatabases}))."\n";
      terminate(9);
    }

  }
}


#checking mode again
if( $annotationMode == 1){ #coord annotation mode
  $suffixForStructure="treedat";
  $genomeIsNeededForBuild=0;
}elsif( $annotationMode == 2){ #interval annotation mode
  $suffixForStructure="treedat";
  $genomeIsNeededForBuild=0
}elsif( $annotationMode == 3){ #SNP annotation mode
  $genomeIsNeededForBuild=1;
  if($databaseName eq "ests" || $databaseName eq "knowngene"){
    print STDERR "Cannot use the ESTs or Knowngene database while using mode 2\n". "usage:\n".$usage."";
    terminate(9);
  }
  $suffixForStructure="treesnpdat";
}elsif( $annotationMode == 4){ #closest TSS
  $genomeIsNeededForBuild=0;
  $suffixForStructure="tssdat";
  if($databaseName eq "ests" || $databaseName eq "knowngene"){
    print STDERR "Cannot use the ESTs or Knowngene database while using mode 2\n". "usage:\n".$usage."";
    terminate(9);
  }
}elsif( $annotationMode == 5){ #insdeltrans
  $genomeIsNeededForBuild=1;
  if($databaseName eq "ests" || $databaseName eq "knowngene"){
    print STDERR "Cannot use the ESTs or Knowngene database while using mode 2\n". "usage:\n".$usage."";
    terminate(9);
  }
  $suffixForStructure="treesnpdat";
  $produceAASeq=1;
}else{
  print STDERR "The -m should either be 1,2,3,4 or 5\n". "usage:\n".$usage."";
  terminate(9);
}

if($databaseName eq "refseq" || $databaseName eq "ensembl"){
  $singleTranscriptDatabase=0;
}else{
  $singleTranscriptDatabase=1;
}
#end mode














####### PARAMETERS WITH DEFAULT VALUES #######
# chromosome file
if( (exists $opts{'c'}) ){
  if($opts{'c'}){
    $dataDirectory = $opts{'c'};
  }else{
    print STDERR "Please specify the directory for storing the chromosome/database files\n";
    terminate(9);
  }
}else{
  $dataDirectory = $ENV{ HOME };
}
# range
if( !(exists $opts{'r'}) ){
  push(@arrayOfRanges,$defaultRange);
}else{
  if( $annotationMode == 4){
    print STDERR "Cannot specify range for mode 4 Find closest 5' end\n";
    terminate(9);
  }
  @arrayOfRanges=split("," , $opts{'r'});
  foreach my $range (@arrayOfRanges){
    if($range !~ /\d+/){
      print STDERR "The ranges must be numerical (-r option)\n";
      terminate(9);
    }elsif($range > $MAXBPFORRANGE){
	print STDERR "The ranges cannot be greater than $MAXBPFORRANGE bp\n";
	terminate(9);
    }
  }
}





####### OPTIONAL PARAMETERS #######
#Annotation options
if($opts{'a'}){
  if( $annotationMode != 3){
    print STDERR "Cannot produce amino acid for any other mode besides 3\n";
    terminate(9);
  }
  $produceAASeq=1;
}

if($opts{'n'}){
  $doNotPrintGenesOut=0;
}else{
  $doNotPrintGenesOut=1;
}

#dbSNP
if($opts{'b'}){
  if( $annotationMode != 3){
    print STDERR "Cannot compare with dbSNP for any other mode besides 3\n";
    terminate(9);
  }

  my $speciesFoundSNP=0;
  foreach my $speciesWithDBSNP (@speciesWithSNPData){
    if($speciesCode eq $speciesWithDBSNP){
      $speciesFoundSNP=1;
    }
  }

  if(!$speciesFoundSNP){
    print STDERR "The species that you have specified does not have snp available\n";
    terminate(9);
  }

  $compareWithDBSNP=1;
}

#do not check ref bp
if($opts{'i'}){
  $noRefBP=1;
}

if($opts{'x'}){
  $doNotPrintXML=1;
}


if($opts{'t'}){
 $printFullFiles=0;
 $doNotPrintXML=1;
}
if($opts{'z'}){
  if($opts{'t'}){
    print STDERR "Cannot use option -z and -t at once\n";
    terminate(9);
  }

  if($annotationMode == 5){
    print STDERR "Cannot use option -z and mode 5 since no files will be printed\n";
    terminate(9);
  }

  $printStatsOnly=1;
}
#bins
if($opts{'g'}){
  if( $annotationMode != 4){
    print STDERR "Cannot specify any other mode than mode 4\n";
    terminate(9);
  }

  if($opts{'g'} =~ /^(\d+)\:(\d+)$/){
    $numberBins=$1;
    $bpInBins  =$2;
  }else{
    print STDERR "The bins must be specified as n:m to report n bins of m base pairs each (ex: -g 10:1000 will create 10 bins of 1kb each side of a tss)\n";
    terminate(9);
  }
}




#Index building options
if($opts{'o'}){
  $onlyBuildStruct=1;
  if($opts{'e'}){
    if(!$customDatabase){
      print STDERR "Please specify a database name not among the following : ".join(" ",(keys %{$availableDatabases}))." using only alphanumerical characters and underscores\n";
      terminate(9);
    }

    if ( $annotationMode == 1 || $annotationMode == 2 ) {
      #fine
    }else{
      print STDERR "The use of the -e option is only available for mode 1 and 2\n";
      terminate(9);
    }
    $useExternalFile=1;
    $listOfExternalFilesPath=$opts{'e'};
  }
}else{
  if($opts{'e'}){
    print STDERR "The -e option is only used upon using -o\n";
    terminate(9);
  }
}

#do not prompt to download
if( $opts{'y'}  ){
  $promptBeforeDownload=0;
}

#do not promt update
#if($opts{'u'}){
#  $promptUpdate=0;
#}

#proxy
if( !(exists $opts{'p'}) || !$opts{'p'}){
  $retrieveDataUCSC = RetrieveDataUCSC->new(ftpTarget  => $ftpTarget);
}else{
  $proxyServer=$opts{'p'};
  $retrieveDataUCSC = RetrieveDataUCSC->new(ftpTarget  => $ftpTarget,
					    httpProxy  => $proxyServer,
					    ftpProxy   => $proxyServer,
					    waisProxy  => $proxyServer);
}





#####################################################

#input file(s)
if (!  $onlyBuildStruct) {

  if ( ( !(exists $opts{'f'}) || !$opts{'f'})  &&
       ( !(exists $opts{'l'}) || !$opts{'l'})  ) {
    print STDERR "Either specify the -f or -l option\n". "usage:\n".$usage."";
    terminate(9);
  } else {
    if ( ( (exists $opts{'f'}) && $opts{'f'} ) &&
	 ( (exists $opts{'l'}) && $opts{'l'} )) {
      print STDERR "Either specify the -f or -l option\n". "usage:\n".$usage."";
      terminate(9);
    }

    if ( (exists $opts{'f'}) && $opts{'f'} ) {
      push(@arrayOfFilesToAnnotate,$opts{'f'});
    }
    if ( (exists $opts{'l'}) && $opts{'l'} ) {
      if(!open(FILEOFFILES,$opts{'l'}) ) { 
	print STDERR "Cannot open file ".$opts{'l'}."\n";
	terminate(9);
      };
      while (my $line = <FILEOFFILES>) {
	chomp($line);
	if(length($line)>0){
	  push(@arrayOfFilesToAnnotate,$line);
	}
      }
      close(FILEOFFILES);
    }

  }
}else{
  if ( ((exists $opts{'f'}) || $opts{'f'}) ||
       ((exists $opts{'l'}) || $opts{'l'})  ){
    print STDERR "Cannot specify the -f or the -l option when building indices\n". "usage:\n".$usage."";
    terminate(9);

  }
}

##########################################
##                                      ##
##                                      ##
##  END READING COMMAND LINE ARGUMENT   ##
##                                      ##
##                                      ##
##########################################

















































#######################################################
#######################################################
####                                               ####
####                                               ####
####          BEGIN INDEX BUILDING MODE            ####
####                                               ####
####                                               ####
#######################################################
#######################################################

if (  $onlyBuildStruct) {

  ############################################
  ##                                        ##
  ##                                        ##
  ##       BEGIN DOWNLOADING GENOME         ##
  ##                                        ##
  ##                                        ##
  ############################################


  if (-d $dataDirectory."/".$directoryInHome."/") {
  } else {
    print "Creating directory ".$dataDirectory."/".$directoryInHome."/\n";
    if (! mkdir($dataDirectory."/".$directoryInHome)) {
      print STDERR "Unable to create directory ".$dataDirectory."/".$directoryInHome." error = ".$!."\n";
      terminate(9);
    }
  }


  if (-d $dataDirectory."/".$directoryInHome."/".$speciesCode) {
    #directory exists, do nothing
  } else {
    print "Creating directory ".$dataDirectory."/".$directoryInHome."/".$speciesCode."\n";
    if (! mkdir($dataDirectory."/".$directoryInHome."/".$speciesCode)) {
      print STDERR "Unable to create directory ".$dataDirectory."/".$directoryInHome."/".$speciesCode." error = ".$!."\n";
      terminate(9);
    }
  }




  my $genomeWasDownloaded=0;
  if ($genomeIsNeededForBuild) {
    print "Genome is needed to build the index for the mode you selected, make sure you do not simply have the chromosome index\n";
  }

  if ( (-e $dataDirectory."/".$directoryInHome."/".$speciesCode."/".$genomeDirectory."/index.txt" ) ) {
    #index exists
  } else {
    if (! (-d $dataDirectory."/".$directoryInHome."/".$speciesCode."/".$genomeDirectory ) ) {
      print "Creating directory ".$dataDirectory."/".$directoryInHome."/".$speciesCode."/".$genomeDirectory."\n";
      if (! mkdir($dataDirectory."/".$directoryInHome."/".$speciesCode."/".$genomeDirectory)) {
	print STDERR "Unable to create directory ".$dataDirectory."/".$directoryInHome."/".$speciesCode."/".$genomeDirectory." error = ".$!."\n";
	terminate(9);
      }
    }


    print "Downloading genome for species:$speciesCode\n";

    #Downloading genome
    testInternetConnection();

    my ($log,$warnLog,$dieLog,@arrayOfChr)=$retrieveDataUCSC->listChromosomes($speciesCode);
    my $sizeTotal=0;
    foreach my $record (@arrayOfChr) {
      $sizeTotal+=$record->{'size'};
    }

    if ($promptBeforeDownload) {
      while (1) {
	print "We will download ".($#arrayOfChr+1)." files for a total of $sizeTotal bytes of data\nIs this ok ? (y/n)\n";
	my $userAnswer=<>;
	chomp($userAnswer);
	if ($userAnswer eq "y") {
	  last;
	} elsif ($userAnswer eq "n") {
	  print "\nexiting\n";
	  exit;
	} else {
	  print "Please answer y or n\n";
	}
      }
    }

    ($log,$warnLog,$dieLog,@arrayOfChr)=$retrieveDataUCSC->downloadGenome($speciesCode);
    $genomeWasDownloaded=1;
    if ($dieLog) {
      print "Downloading genome for species:$speciesCode failed\n";
    }

    foreach my $chrToMove (@arrayOfChr) {
      my $command="/bin/mv -f $chrToMove ".$dataDirectory."/".$directoryInHome."/".$speciesCode."/".$genomeDirectory."/";
      my $output=`$command`;
      if ($output) {
	print STDERR "Command $command failed\n";
	terminate(9);
      }

      $command="$splitFastaPath"."splitFasta.pl -f ".$dataDirectory."/".$directoryInHome."/".$speciesCode."/".$genomeDirectory."/".$chrToMove." -e";
      print "Examing $chrToMove\n";
      $output=`$command`;
      if ($output) {
	print STDERR "Command $command failed\n";
	terminate(9);
      }


    }

  }




  if (-e $dataDirectory."/".$directoryInHome."/".$speciesCode."/".$genomeDirectory."/index.txt") {
    if ($genomeWasDownloaded) {
      print STDERR "Error, the file ".$dataDirectory."/".$directoryInHome."/".$speciesCode."/".$genomeDirectory."/index.txt should not be there, please delete the partially genome and try again\n";
      terminate(9);
    }
  } else {
    print "Creating chromosome index, please wait, this may take a few minutes ...\n";
    my $commandChrIndex=$createChrIndexPath."createChromosomeIndex.pl  ".$dataDirectory."/".$directoryInHome."/".$speciesCode."/".$genomeDirectory."/  >  ".$dataDirectory."/".$directoryInHome."/".$speciesCode."/".$genomeDirectory."/index.txt";
    my $outputChrIndex=`$commandChrIndex`;
    if ($outputChrIndex) {
      print STDERR "Command $commandChrIndex failed\n";
      terminate(9);
    }
    print "\n...done\n";
  }




  ############################################
  ##                                        ##
  ##                                        ##
  ##         END DOWNLOADING GENOME         ##
  ##                                        ##
  ##                                        ##
  ############################################









  ############################################
  ##                                        ##
  ##                                        ##
  ##    BEGIN  DOWNLOADING   DATABASES      ##
  ##                                        ##
  ##                                        ##
  ############################################

  my $stringIdForFreeze;
  my $stringIdForFreezeSNP;

  ############################################
  ##    DOWNLOAD DATABASES FROM UCSC        ##
  ############################################

  if (!$useExternalFile) {

    my @arrayDatabasesToDownload;

    $stringIdForFreeze=$databaseName;
    $stringIdForFreezeSNP="snp";


    if (-d $dataDirectory."/".$directoryInHome."/".$speciesCode."/".$databaseDirectory) {
      #directory exists
    } else {
      print "Creating directory ".$dataDirectory."/".$directoryInHome."/".$speciesCode."/".$databaseDirectory."\n";
      if (! mkdir($dataDirectory."/".$directoryInHome."/".$speciesCode."/".$databaseDirectory)) {
	print STDERR "Unable to create directory ".$dataDirectory."/".$directoryInHome."/".$speciesCode."/".$databaseDirectory." error = ".$!."\n";
	terminate(9);
      }
    }


    if (! (-d $dataDirectory."/".$directoryInHome."/".$speciesCode."/".$structureDirectory) ) {
      print "Creating directory ".$dataDirectory."/".$directoryInHome."/".$speciesCode."/".$structureDirectory."\n";
      if (! mkdir($dataDirectory."/".$directoryInHome."/".$speciesCode."/".$structureDirectory)) {
	print STDERR "Unable to create directory ".$dataDirectory."/".$directoryInHome."/".$speciesCode."/".$structureDirectory." error = ".$!."\n";
	terminate(9);
      }

      #    $createStructure=1;
    } else {

      if (  (-e $dataDirectory."/".$directoryInHome."/".$speciesCode."/".$structureDirectory."/".$stringIdForFreeze.".".$suffixForStructure)  ) {

	if ($promptBeforeDownload) {
	  while (1) {
	    print "An index file ".$dataDirectory."/".$directoryInHome."/".$speciesCode."/".$structureDirectory."/".$stringIdForFreeze.".".$suffixForStructure." already exists, are you sure you want to overwrite the file ? (y/n)\n";
	    my $userAnswer=<>;
	    chomp($userAnswer);

	    if ($userAnswer eq "y") {
	      last;
	    } elsif ($userAnswer eq "n") {
	      print "exiting\n";
	      terminate(9);
	    } else {
	      print "Please answer y or n\n";
	    }
	  }
	}

      }

    }




    foreach my $databaseToDownload (@arrayOfDatabases) {
      push(@arrayDatabasesToDownload,$databaseToDownload);
    }


    if ($#arrayDatabasesToDownload >= 0) {

      if ( $#arrayDatabasesToDownload != $#arrayOfDatabases) {
	print "Cannot download only part of the databases\n";
	exit;
      }

      #Downloading database
      testInternetConnection();

      my $sizeTotal=0;
      my $numberOfFiles=0;

      foreach my $databaseToDownload (@arrayDatabasesToDownload) {
	my ($log,$warnLog,$dieLog,$sizeFile);
	($log,$warnLog,$dieLog,$sizeFile)=$retrieveDataUCSC->getSizeDatabase($speciesCode,$databaseToDownload);
	$sizeTotal=$sizeFile;
	$numberOfFiles++;
      }

      if ($promptBeforeDownload) {
	while (1) {
	  print "We will download $numberOfFiles files for a total of $sizeTotal bytes of data\nIs this ok ? (y/n)\n";
	  my $userAnswer=<>;
	  chomp($userAnswer);
	  if ($userAnswer eq "y") {
	    last;
	  } elsif ($userAnswer eq "n") {
	    print "\nexiting\n";
	    exit;
	  } else {
	    print "Please answer y or n\n";
	  }
	}
      }



      foreach my $databaseToDownload (@arrayDatabasesToDownload) {


	my ($log,$warnLog,$dieLog)=$retrieveDataUCSC->downloadDatabase($speciesCode,$databaseToDownload);



	##Delete previous database files
	opendir(TARGETDIR, $dataDirectory."/".$directoryInHome."/".$speciesCode."/".$databaseDirectory);

	foreach my $file (readdir(TARGETDIR)) {
	  if ($file =~ /^$databaseToDownload/) {
	    print "removing previous file ".$dataDirectory."/".$directoryInHome."/".$speciesCode."/".$databaseDirectory."/".$file."\n";
	    if(!unlink($dataDirectory."/".$directoryInHome."/".$speciesCode."/".$databaseDirectory."/".$file)){
	      print STDERR "Could not delete ".$dataDirectory."/".$directoryInHome."/".$speciesCode."/".$databaseDirectory."/".$file."\n";
	    }
	  }

	}
	closedir(TARGETDIR);

	my $command="/bin/mv -f ".$databaseToDownload.".txt  ".$dataDirectory."/".$directoryInHome."/".$speciesCode."/".$databaseDirectory."/";
	my $output=`$command`;
	if ($output) {
	  print STDERR "Command $command failed\n";
	  terminate(9);
	}
      }

    }


    foreach my $databaseFileName (@arrayOfDatabases) {

      if (-e $dataDirectory."/".$directoryInHome."/".$speciesCode."/".$databaseDirectory."/".$databaseFileName.".txt") {

	my $st;
	if (!($st = stat($dataDirectory."/".$directoryInHome."/".$speciesCode."/".$databaseDirectory."/".$databaseFileName.".txt") ) ) {
	  print STDERR "Could not stat ".$dataDirectory."/".$directoryInHome."/".$speciesCode."/".$databaseDirectory."/".$databaseFileName.".txt: $!";
	  terminate(9);
	}
	my $mtime = $st->mtime;

	#$stringIdForFreeze.=$mtime;
	push(@{$databaseFiles},$databaseFileName.".txt");
	push(@{$databaseFiles}, scalar(localtime($st->mtime)) );

      } else {
	print "Error, the database $databaseFileName does not exist\n";
      }

    }




    #### BEGIN DBSNP ####
    if ($compareWithDBSNP) {


      my $command="/bin/ls -1 ".$dataDirectory."/".$directoryInHome."/".$speciesCode."/".$databaseDirectory."/";
      my $output=`$command`;

      if ($?) {
	print STDERR "Command $command failed\n";
	terminate(9);
      }

      my @filesInDatabase = split("\n",$output);
      my ($mostRecentFile,@arrayOfOlddbSNP) = findNewestSNPFile(@filesInDatabase);
      my $downloadNewDBSNP=1;

      if ($mostRecentFile ne "-1" ) {
	my $st;
	if ( !($st= stat($dataDirectory."/".$directoryInHome."/".$speciesCode."/".$databaseDirectory."/".$mostRecentFile) )) {
	  print STDERR "Could not stat ".$dataDirectory."/".$directoryInHome."/".$speciesCode."/".$databaseDirectory."/".$mostRecentFile.": $!";
	  terminate(9);
	}
	my $mtime = $st->mtime;
	$dbSNPFile=$dataDirectory."/".$directoryInHome."/".$speciesCode."/".$databaseDirectory."/".$mostRecentFile;

	if ((time()-$mtime) >  $numberOfSecondsToPromptUpdate) {
	  print "The existing file $dbSNPFile is more than 6 months old\n";
	}else{
	  print "A file $dbSNPFile already exists\n";

	}

	if ($promptBeforeDownload) {
	  while (1) {
	    print "Do you want to download a new version of the file ? (y/n)\n";
	    my $userAnswer=<>;
	    chomp($userAnswer);
	    if ($userAnswer eq "y") {
	      $downloadNewDBSNP=1;
	      last;
	    } elsif ($userAnswer eq "n") {
	      $downloadNewDBSNP=0;
	      last;
	    } else {
	      print "Please answer y or n\n";
	    }
	  }
	}

	foreach my $oldFileToDelete (@arrayOfOlddbSNP) {
	  if ($promptBeforeDownload) {
	    while (1) {
	      print "An older file ".$dataDirectory."/".$directoryInHome."/".$speciesCode."/".$databaseDirectory."/".$oldFileToDelete." already exists, do you wish to delete it ? (y/n)\n";
	      my $userAnswer=<>;
	      chomp($userAnswer);
	      if ($userAnswer eq "y") {
		if(!unlink($dataDirectory."/".$directoryInHome."/".$speciesCode."/".$databaseDirectory."/".$oldFileToDelete)){
		  print STDERR "Could not delete ".$dataDirectory."/".$directoryInHome."/".$speciesCode."/".$databaseDirectory."/".$oldFileToDelete."\n";
		}
		last;
	      } elsif ($userAnswer eq "n") {
		last;
	      } else {
		print "Please answer y or n\n";
	      }
	    }
	  }
	}




	#	if ($promptUpdate) {
	#	  if ((time()-$mtime) >  $numberOfSecondsToPromptUpdate) {
	#
	#	    if ( $onlyBuildStruct ) {
	#	      $downloadNewDBSNP=1;
	#	    } else {
	#	      while (1) {
	#		print "The file $mostRecentFile is more than 6 months old\nDo you want to download a new version of the file ? (y/n)\n";
	#		my $userAnswer=<>;
	#		chomp($userAnswer);
	#		if ($userAnswer eq "y") {
	#		  $downloadNewDBSNP=1;
	#		  last;
	#		} elsif ($userAnswer eq "n") {
	#		  last;
	#		} else {
	#		  print "Please answer y or n\n";
	#		}
	#	      }
	#	    }
	#	  }
	#	}

      } else {			#if no file, download
	$downloadNewDBSNP=1;
      }


      if ( $downloadNewDBSNP ) {
	testInternetConnection();

	my ($log,$warnLog,$dieLog,@arrayOfDatabases)=$retrieveDataUCSC->listDatabases($speciesCode);
	my @arrayOfDBSNPfiles;
	foreach my $databaseHash (@arrayOfDatabases) {
	  if ($databaseHash->{'name'} =~ /^snp\d+$/) {
	    push(@arrayOfDBSNPfiles,$databaseHash->{'name'});
	  }
	}
	my ($mostRecentFileToDownload,@arrayUCSCdbSNPold) = findNewestSNPFile(@arrayOfDBSNPfiles);
	my $mostRecentFileToDownloadSize=-1;

	foreach my $databaseHash (@arrayOfDatabases) {
	  if ($databaseHash->{'name'} eq $mostRecentFileToDownload) {
	    $mostRecentFileToDownloadSize = $databaseHash->{'size'};
	  }
	}

	if ( $mostRecentFileToDownloadSize == -1 ) {
	  print STDERR "Unable to find dbSNP in the UCSC database\n";
	  terminate(9);
	}

	if ($promptBeforeDownload) {
	  while (1) {
	    print "We will download 1 file for a total of $mostRecentFileToDownloadSize bytes of data\nIs this ok ? (y/n)\n";
	    my $userAnswer=<>;
	    chomp($userAnswer);
	    if ($userAnswer eq "y") {
	      last;
	    } elsif ($userAnswer eq "n") {
	      print "\nexiting\n";
	      exit;
	    } else {
	      print "Please answer y or n\n";
	    }
	  }
	}

	($log,$warnLog,$dieLog) = $retrieveDataUCSC->downloadDatabase($speciesCode,$mostRecentFileToDownload.".txt.gz");
	my $command="/bin/mv -f ".$mostRecentFileToDownload.".txt  ".$dataDirectory."/".$directoryInHome."/".$speciesCode."/".$databaseDirectory."/";
	my $output=`$command`;
	if ($output) {
	  print STDERR "Command $command failed\n";
	  terminate(9);
	}

	my $st;
	if (!($st = stat($dataDirectory."/".$directoryInHome."/".$speciesCode."/".$databaseDirectory."/".$mostRecentFileToDownload.".txt") )) {
	  print STDERR "Could not stat ".$dataDirectory."/".$directoryInHome."/".$speciesCode."/".$databaseDirectory."/".$mostRecentFileToDownload.".txt: $!";
	  terminate(9);
	}
	my $mtime = $st->mtime;

	$dbSNPFile=$dataDirectory."/".$directoryInHome."/".$speciesCode."/".$databaseDirectory."/".$mostRecentFileToDownload.".txt";


	push(@{$databaseFiles},$mostRecentFileToDownload.".txt");
	push(@{$databaseFiles}, scalar(localtime($st->mtime)) );

      } else {			# we do not download
	my $st;
	if (!($st = stat($dataDirectory."/".$directoryInHome."/".$speciesCode."/".$databaseDirectory."/".$mostRecentFile) )) {
	  print STDERR "Could not stat ".$dataDirectory."/".$directoryInHome."/".$speciesCode."/".$databaseDirectory."/".$mostRecentFile.": $!";
	  terminate(9);
	}
	my $mtime = $st->mtime;
	#$stringIdForFreezeSNP.=$mtime;

	push(@{$databaseFiles},$mostRecentFile);
	push(@{$databaseFiles}, scalar(localtime($st->mtime)) );
      }

    } #end compareWithDBSNP
    ####  END DBSNP ####




    ############################################
    ##                                        ##
    ##                                        ##
    ##     END   DOWNLOADING   DATABASES      ##
    ##                                        ##
    ##                                        ##
    ############################################






    my $startTIME2 = time;















    ############################################
    ##                                        ##
    ##                                        ##
    ##     BEGIN READING GENE DATABASES       ##
    ##                                        ##
    ##                                        ##
    ############################################

    ChromosomeIndexing::setBaseDir($dataDirectory."/".$directoryInHome."/");






    my $logfile = $dataDirectory."/".$directoryInHome."/".$speciesCode."/".$structureDirectory."/".$stringIdForFreeze.".log";

    open(STDERR, ">$logfile");

    print "Reading gene files...\n";
    my $genesObject;

    if ($databaseName eq "refseq") {
      $genesObject = CreateRefSeqGenes->new(species => $speciesCode,
					    baseDir => $dataDirectory."/".$directoryInHome."/");
    } elsif ($databaseName eq "ensembl") {
      $genesObject = CreateEnsemblGenes->new(species => $speciesCode,
					     baseDir => $dataDirectory."/".$directoryInHome."/");
    } elsif ($databaseName eq "knowngene") {
      $genesObject = CreateKnownGenes->new(species => $speciesCode,
					   baseDir => $dataDirectory."/".$directoryInHome."/");
    } elsif ($databaseName eq "ests") {
      $genesObject = CreateEsts->new(species => $speciesCode,
				     baseDir => $dataDirectory."/".$directoryInHome."/");
    } else {
      print STDERR "Invalid database name\n";
      terminate(9);
    }



    print "...done\n";



    ############################################
    ##                                        ##
    ##                                        ##
    ##       END READING GENE DATABASES       ##
    ##                                        ##
    ##                                        ##
    ############################################













    ############################################
    ##                                        ##
    ##                                        ##
    ##     BEGIN BUILDING TREE STRUCTURE      ##
    ##                                        ##
    ##                                        ##
    ############################################


    my $indexOfChrBeingBuilt=1;
    my $indexOfChrMAX=scalar(ChromosomeIndexing::listChr($speciesCode))+1;

    print "Building tree structure ...\n";
    #For each chromosome

    foreach my $indexTree (ChromosomeIndexing::listChr($speciesCode)) {

      #Prepare array of transcripts
      my $arrayOfR   =[];
      my $arrayOfTSSForChr =[];

      print "\nBuilding for $indexTree ".$indexOfChrBeingBuilt." of ".($indexOfChrMAX-1)."\n";

      my @arrayOfTranscriptsForChr=$genesObject->getArrayOfTranscriptsForGivenChr($indexTree);

      if ($#arrayOfTranscriptsForChr == -1) {
	#printing full progress bar
	print "\r";
	print "\t";
	print  RetrieveDataUCSC::printBar(100,100);
      }

      my $minBinIndexFound;
      #print "chr ".$arrayOfTranscriptsForChr[0]->getChromosome()."\n";

      #For each transcript for that chromosome
      for (my $i=0;$i<=$#arrayOfTranscriptsForChr;$i++) {

	#printing progress bar
	print "\r";
	print "\t";
	print  RetrieveDataUCSC::printBar($i,$#arrayOfTranscriptsForChr);


	if ( $annotationMode == 3 || $annotationMode == 5 ) { #SNP & indel/trans annotation mode
	  $arrayOfTranscriptsForChr[$i]->buildExonCoding();
	}

	my $indexbinst=int($arrayOfTranscriptsForChr[$i]->getTranscriptionRange()->getStart()/$genomebinsize);
	my $indexbinen=int($arrayOfTranscriptsForChr[$i]->getTranscriptionRange()->getEnd()/$genomebinsize);
	#print "id     ".$arrayOfTranscriptsForChr[$i]->getId()."\n";
	#print "range  ".$arrayOfTranscriptsForChr[$i]->getTranscriptionRange()->print()."\n";
	#print "indexbinst $indexbinst\n";
	#print "indexbinen $indexbinen\n";

	if($i==0){
	  $minBinIndexFound=$indexbinst;
	}else{
	  if($indexbinst<$minBinIndexFound){
	    $minBinIndexFound=$indexbinst;
	  }
	}

	for (my $indexbin=$indexbinst;$indexbin<=$indexbinen;$indexbin++) {
	  #perhaps modify to only include the limits of the bins
	  my $segmentTreeRange=new SegmentTreeRange(start => $arrayOfTranscriptsForChr[$i]->getTranscriptionRange()->getStart(),
						    end   => $arrayOfTranscriptsForChr[$i]->getTranscriptionRange()->getEnd()    );
	  my $hash={transcript    =>  $arrayOfTranscriptsForChr[$i],
		    range         =>  $segmentTreeRange};
	  if($arrayOfR->[$indexbin]){
	    push(@{$arrayOfR->[$indexbin]},$hash);
	  }else{
	    $arrayOfR->[$indexbin]=[$hash];
	  }
	}

	my $hashTSS={transcript    =>  $arrayOfTranscriptsForChr[$i],
		     tss           =>  $arrayOfTranscriptsForChr[$i]->get5Prime() };
	push(@{$arrayOfTSSForChr},$hashTSS);

      } #for each transcript

      if($#{$arrayOfR} != -1){ #if we found something
	$arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$indexTree)]=[];
	my $segTree=new SegmentTree(arrayOfInputs => $arrayOfR->[$minBinIndexFound]);
	$arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$indexTree)]->[$minBinIndexFound]=$segTree;

	for (my $indexbin=($minBinIndexFound+1);$indexbin<=$#{$arrayOfR};$indexbin++) {
	  if($arrayOfR->[$indexbin]){
	    my $segTree=new SegmentTree(arrayOfInputs => $arrayOfR->[$indexbin]);
	    $arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$indexTree)]->[$indexbin]=$segTree;
	  }
	}
      }

      my @arrayTempTSS=sort {$a->{'tss'} <=> $b->{'tss'} } @{$arrayOfTSSForChr};
      $arrayOfTSSForChr=\@arrayTempTSS;

      $arrayOfTSS->[ChromosomeIndexing::chr2index($speciesCode,$indexTree)]=$arrayOfTSSForChr;
      $indexOfChrBeingBuilt++;

    } #for each chromosome
    print "...done\n";

    close(STDERR);

    ##Delete previous files with $suffixForStructure and .tssdat
    opendir(TARGETDIR, $dataDirectory."/".$directoryInHome."/".$speciesCode."/".$structureDirectory."/");
    foreach my $file (readdir(TARGETDIR)) {
      if ($file =~ /^$databaseName\d+\.$suffixForStructure$/) {
	print "removing previous file ".$dataDirectory."/".$directoryInHome."/".$speciesCode."/".$structureDirectory."/".$file."\n";
	if(!unlink($dataDirectory."/".$directoryInHome."/".$speciesCode."/".$structureDirectory."/".$file)){
	  print STDERR "Could not delete ".$dataDirectory."/".$directoryInHome."/".$speciesCode."/".$structureDirectory."/".$file."\n";
	}
      }

      if ($file =~ /^$databaseName\d+\.log$/ && $file !~ /$stringIdForFreeze/) {
	print "removing previous file ".$dataDirectory."/".$directoryInHome."/".$speciesCode."/".$structureDirectory."/".$file."\n";
	if(!unlink($dataDirectory."/".$directoryInHome."/".$speciesCode."/".$structureDirectory."/".$file)){
	  print STDERR "Could not delete ".$dataDirectory."/".$directoryInHome."/".$speciesCode."/".$structureDirectory."/".$file."\n";
	}
      }



      if ($file =~ /^$databaseName\d+\.tssdat$/) {
	print "removing previous file ".$dataDirectory."/".$directoryInHome."/".$speciesCode."/".$structureDirectory."/".$file."\n";
	if(!unlink($dataDirectory."/".$directoryInHome."/".$speciesCode."/".$structureDirectory."/".$file)){
	  print STDERR "Could not delete ".$dataDirectory."/".$directoryInHome."/".$speciesCode."/".$structureDirectory."/".$file."\n";
	}
      }
    }
    closedir(TARGETDIR);

    my $fileOutFreeze ;
    $fileOutFreeze        = $dataDirectory."/".$directoryInHome."/".$speciesCode."/".$structureDirectory."/".$stringIdForFreeze.".".$suffixForStructure;
    my $fileOutFreezeTSS  = $dataDirectory."/".$directoryInHome."/".$speciesCode."/".$structureDirectory."/".$stringIdForFreeze.".tssdat";

    if (!store($arrayOfTrees, $fileOutFreeze) ) {
      print STDERR "Can't store tree in $fileOutFreeze!\n";
      terminate(9);
    }

    if (!store($arrayOfTSS, $fileOutFreezeTSS) ) {
      print STDERR "Can't store tree in $fileOutFreezeTSS!\n";
      terminate(9);
    }


    print "done saving to file\n";
    ############################################
    ##                                        ##
    ##                                        ##
    ##      END BUILDING TREE STRUCTURE       ##
    ##                                        ##
    ##                                        ##
    ############################################






    if ($compareWithDBSNP) {
      my $createdbsnpDatastructure=1;

      if (  (-e $dataDirectory."/".$directoryInHome."/".$speciesCode."/".$structureDirectory."/".$stringIdForFreezeSNP.".dbsnp")  ) {

	if ($promptBeforeDownload) {
	  while (1) {
	    print "An index file ".$dataDirectory."/".$directoryInHome."/".$speciesCode."/".$structureDirectory."/".$stringIdForFreezeSNP.".dbsnp already exists, are you sure you want to overwrite the file ? (y/n)\n";
	    my $userAnswer=<>;
	    chomp($userAnswer);
	    if ($userAnswer eq "y") {
	      $createdbsnpDatastructure=1;
	      last;
	    } elsif ($userAnswer eq "n") {
	      $createdbsnpDatastructure=0;
	      last;
	    } else {
	      print "Please answer y or n\n";
	    }
	  }
	}

      }

      if ($createdbsnpDatastructure) {
	print "Building dbsnp datastructure ...\n";
	#build data
	foreach my $chrInSpecies (ChromosomeIndexing::listChr($speciesCode)) {
	  my $dbSNPstruct=SNPstructure->new(genomebinsize => $genomebinsize);
	  $arrayOfSNPStructures->[ChromosomeIndexing::chr2index($speciesCode,$chrInSpecies)] = $dbSNPstruct;
	}

	my $fileOutFreezeSNPlog  = $dataDirectory."/".$directoryInHome."/".$speciesCode."/".$structureDirectory."/".$stringIdForFreezeSNP.".dbsnp.log";
	if (!open(SNPFILElog, ">".$fileOutFreezeSNPlog) ) {
	  print STDERR "Can't open log file for dbSNP : $fileOutFreezeSNPlog\n";
	  terminate(9);
	}

	if (!open(SNPFILE, $dbSNPFile)) {
	  print STDERR "Can't open file for dbSNP : $dbSNPFile\n";
	  terminate(9);
	}

	my $numberOfTotalLinesSNP =0;
	my $currentNumberOfTotalLinesSNP =0;

	while (my $line = <SNPFILE>) {
	  $numberOfTotalLinesSNP++;
	}
	close(SNPFILE);

	if (!open(SNPFILE, $dbSNPFile)) {
	  print STDERR "Can't open file for dbSNP : $dbSNPFile\n";
	  terminate(9);
	}

	while (my $line = <SNPFILE>) {
	  if( (($currentNumberOfTotalLinesSNP%250000) == 0) ||
	      ($numberOfTotalLinesSNP-$currentNumberOfTotalLinesSNP)<100 ){
	    print "\r";
	    print "\t";
	    print  RetrieveDataUCSC::printBar($currentNumberOfTotalLinesSNP,$numberOfTotalLinesSNP);
	  }
	  $currentNumberOfTotalLinesSNP++;

	  my @arrayOfFields=split("\t",$line);
	  my $returnLine = $arrayOfSNPStructures->[ChromosomeIndexing::chr2index($speciesCode,$arrayOfFields[1])]->addLine($line);
	  print SNPFILElog $returnLine;
	}
	close(SNPFILE);

	close(SNPFILElog);

	my $fileOutFreezeSNP  = $dataDirectory."/".$directoryInHome."/".$speciesCode."/".$structureDirectory."/".$stringIdForFreezeSNP.".dbsnp";

	if (!store($arrayOfSNPStructures, $fileOutFreezeSNP)) {
	  print STDERR "Can't store tree in $fileOutFreezeSNP!\n";
	  terminate(9);
	}
	print "\ndone\n";
      }


    }

    ###################################################################
    ##  DO NOT  DOWNLOAD DATABASES FROM UCSC AND USE LOCAL FILE      ##
    ###################################################################
  } else {			#if $useExternalFile
    $stringIdForFreeze=$databaseName;

    ChromosomeIndexing::setBaseDir($dataDirectory."/".$directoryInHome."/");

    if ($compareWithDBSNP) {
      print STDERR "Cannot specify dbSNP when using the file to build the segment tree";
      terminate(9);
    }
    #store data
    my $arrayOfInitialData=[];
    my %hashOfIDs;

    foreach my $indexTree (ChromosomeIndexing::listIndexChr($speciesCode)) {
      $arrayOfInitialData->[$indexTree]=[];
    }


    #read psl file

    foreach my $externalFilePath ( split(",",$listOfExternalFilesPath) ) {
      if ( !open(PSLFILE,$externalFilePath) ) {
	print STDERR "Cannot open file ".$externalFilePath."\n";
	terminate(9);
      }


      while (my $line = <PSLFILE>) {
	chomp($line);
	my @arrayOfFields=split("\t",$line);
	if (substr($line,0,1) eq "#") { #skip header or comments
	  next;
	}

	my $idForRecord;
	my $chrForRecord;
	my $startForRecord;
	my $endForRecord;
	my $strandForRecord;

	if ( ($arrayOfFields[5] eq "+" ) ||
	     ($arrayOfFields[5] eq "-" )  ) { #bed no bin field

	  $strandForRecord=  $arrayOfFields[5];
	  $idForRecord    =  $arrayOfFields[3];
	  $chrForRecord   =  $arrayOfFields[0];
	  $startForRecord =  $arrayOfFields[1];
	  $endForRecord   =  $arrayOfFields[2];
	} elsif ( ($arrayOfFields[6] eq "+" ) ||
		  ($arrayOfFields[6] eq "-" )  ) { #bed with bin field
	  $strandForRecord=  $arrayOfFields[6];
	  $idForRecord    =  $arrayOfFields[4];
	  $chrForRecord   =  $arrayOfFields[1];
	  $startForRecord =  $arrayOfFields[2];
	  $endForRecord   =  $arrayOfFields[3];

	} elsif ( ($arrayOfFields[8] eq "+" ) ||
		  ($arrayOfFields[8] eq "-" )  ) { #psl no bin field
	  $strandForRecord=  $arrayOfFields[8];
	  $idForRecord    =  $arrayOfFields[9];
	  $chrForRecord   =  $arrayOfFields[13];
	  $startForRecord =  $arrayOfFields[15];
	  $endForRecord   =  $arrayOfFields[16];
	}elsif ( ($arrayOfFields[9] eq "+" ) ||
		 ($arrayOfFields[9] eq "-" )  ) { #psl with bin field
	  $strandForRecord=  $arrayOfFields[9];
	  $idForRecord    =  $arrayOfFields[10];
	  $chrForRecord   =  $arrayOfFields[14];
	  $startForRecord =  $arrayOfFields[16];
	  $endForRecord   =  $arrayOfFields[17];
	} else {
	  print STDERR "Format error with line $line in file $externalFilePath\ncheck to make sure the file is in psl or bed (with strand) format\n";
	  terminate(9);
	}

	if (exists $hashOfIDs{$idForRecord."#".$chrForRecord."#".$startForRecord."#".$endForRecord."#".$strandForRecord}) {
	  print STDERR "Warning for line $line in file $externalFilePath, id: ".$idForRecord."#".$chrForRecord."#".$startForRecord."#".$endForRecord."#".$strandForRecord." was previously found, skipping\n";
	  #terminate(9);
	} else {
	  $hashOfIDs{$idForRecord."#".$chrForRecord."#".$startForRecord."#".$endForRecord."#".$strandForRecord} = 0;
	}

	my $hashRange=new Range(start => $startForRecord,
				end   => $endForRecord);
	my $hashOfData={'id'     => $idForRecord,
			'range'  => $hashRange,
			'strand' => $strandForRecord};
	push(@{$arrayOfInitialData->[ChromosomeIndexing::chr2index($speciesCode,$chrForRecord)]},$hashOfData);

      }# end for each line
      close(PSLFILE);
    }				#end for each external file


    #done reading files, build tree
    my $indexOfChrBeingBuilt=1;
    my $indexOfChrMAX=scalar(ChromosomeIndexing::listChr($speciesCode))+1;

    print "Building tree structure ...\n";
    #For each chromosome
    foreach my $indexTree (ChromosomeIndexing::listIndexChr($speciesCode)) {

      #Prepare array of transcripts
      my $arrayOfR   =[];

      print "\nBuilding for ".ChromosomeIndexing::index2chr($speciesCode,$indexTree)." ".$indexOfChrBeingBuilt." of ".($indexOfChrMAX-1)."\n";

      my @arrayOfRecordsForChr= @{$arrayOfInitialData->[$indexTree]};

      if ($#arrayOfRecordsForChr == -1) {
	#printing full progress bar
	print "\r";
	print "\t";
	print  RetrieveDataUCSC::printBar(100,100);
      } else {
	my $minBinIndexFound;
	#For each transcript for that chromosome
	for (my $i=0;$i<=$#arrayOfRecordsForChr;$i++) {
	  #printing progress bar
	  print "\r";
	  print "\t";
	  print  RetrieveDataUCSC::printBar($i,$#arrayOfRecordsForChr);

	  my $indexbinst=int($arrayOfRecordsForChr[$i]->{'range'}->getStart()/$genomebinsize);
	  my $indexbinen=int($arrayOfRecordsForChr[$i]->{'range'}->getEnd()/$genomebinsize);

	  if ($i==0) {
	    $minBinIndexFound=$indexbinst;
	  } else {
	    if ($indexbinst<$minBinIndexFound) {
	      $minBinIndexFound=$indexbinst;
	    }
	  }

	  for (my $indexbin=$indexbinst;$indexbin<=$indexbinen;$indexbin++) {

	    my $segmentTreeRange=new SegmentTreeRange(start => $arrayOfRecordsForChr[$i]->{'range'}->getStart(),
						      end   => $arrayOfRecordsForChr[$i]->{'range'}->getEnd());
	    my $hash={record        =>  $arrayOfRecordsForChr[$i],
		      range         =>  $segmentTreeRange};


	    if ($arrayOfR->[$indexbin]) {
	      push(@{$arrayOfR->[$indexbin]},$hash);
	    } else {
	      $arrayOfR->[$indexbin]=[$hash];
	    }
	  }
	}
	#end for each transcripts

	if ($#{$arrayOfR} != -1) { #if we found something
	  $arrayOfTrees->[$indexTree]=[];
	  my $segTree=new SegmentTree(arrayOfInputs => $arrayOfR->[$minBinIndexFound]);
	  $arrayOfTrees->[$indexTree]->[$minBinIndexFound]=$segTree;

	  for (my $indexbin=($minBinIndexFound+1);$indexbin<=$#{$arrayOfR};$indexbin++) {
	    if ($arrayOfR->[$indexbin]) {
	      my $segTree=new SegmentTree(arrayOfInputs => $arrayOfR->[$indexbin]);
	      $arrayOfTrees->[$indexTree]->[$indexbin]=$segTree;
	    }
	  }
	}

	$indexOfChrBeingBuilt++;
      }
    }#end for each index chromosome
    print "...done\n";

    #storing
    my $fileOutFreeze ;
    $fileOutFreeze        = $dataDirectory."/".$directoryInHome."/".$speciesCode."/".$structureDirectory."/".$stringIdForFreeze.".".$suffixForStructure;

    if (!store($arrayOfTrees, $fileOutFreeze) ) {
      print STDERR "Can't store tree in $fileOutFreeze!\n"; terminate(9);
    }




  }				#end not ucsc but local file


  terminate(9);








  #######################################################
  #######################################################
  ####                                               ####
  ####                                               ####
  ####          END  INDEX BUILDING MODE             ####
  ####                                               ####
  ####                                               ####
  #######################################################
  #######################################################

} else {





  my $startTIME = time;
  #######################################################
  #######################################################
  ####                                               ####
  ####                                               ####
  ####           BEGIN ANNOTATION MODE               ####
  ####                                               ####
  ####                                               ####
  #######################################################
  #######################################################

  if (-d $dataDirectory."/".$directoryInHome."/") {
  } else {
    print STDERR "Missing directory ".$dataDirectory."/".$directoryInHome." \n";
    terminate(9);
  }


  if (-d $dataDirectory."/".$directoryInHome."/".$speciesCode) {
  } else {
    print STDERR "Error: missing directory ".$dataDirectory."/".$directoryInHome."/".$speciesCode." \n";
    terminate(9);
  }

  if (-d $dataDirectory."/".$directoryInHome."/".$speciesCode."/".$genomeDirectory) {
  } else {
    print STDERR "Error: missing directory ".$dataDirectory."/".$directoryInHome."/".$speciesCode."/".$genomeDirectory." \n";
    terminate(9);
  }


  if (-e $dataDirectory."/".$directoryInHome."/".$speciesCode."/".$genomeDirectory."/index.txt" ) {
  } else {
    print STDERR "Error: missing chromosome index file for $speciesCode ".$dataDirectory."/".$directoryInHome."/".$speciesCode." \n";
    terminate(9);
  }
  ChromosomeIndexing::setBaseDir($dataDirectory."/".$directoryInHome."/");


  my $stringIdForFreeze=$databaseName;
  my $stringIdForFreezeSNP="snp";

  my $fileOutFreeze        = $dataDirectory."/".$directoryInHome."/".$speciesCode."/".$structureDirectory."/".$stringIdForFreeze.".".$suffixForStructure;
  #my $fileOutFreezeTSS  = $dataDirectory."/".$directoryInHome."/".$speciesCode."/".$structureDirectory."/".$stringIdForFreeze.".tssdat";

  if (-d $dataDirectory."/".$directoryInHome."/".$speciesCode."/".$structureDirectory) {
  } else {
    print STDERR "Error: missing directory ".$dataDirectory."/".$directoryInHome."/".$speciesCode."/".$structureDirectory." \n";
    terminate(9);
  }

  if (-e $fileOutFreeze) {
  } else {
    print STDERR "The required index file ".$fileOutFreeze." does not exist, please download it from our website or build a new one\n";
    terminate(9);
  }



  if (  (-e $fileOutFreeze) ) {
    print "Retrieving tree file ...";
    $arrayOfTrees = retrieve($fileOutFreeze);
    print "done\n";
  }


  if ($compareWithDBSNP) {
    if (  (-e $dataDirectory."/".$directoryInHome."/".$speciesCode."/".$structureDirectory."/".$stringIdForFreezeSNP.".dbsnp") ) {
      #retrieve data
      print "Retrieving dbsnp data ...";
      my $fileOutFreezeSNP  = $dataDirectory."/".$directoryInHome."/".$speciesCode."/".$structureDirectory."/".$stringIdForFreezeSNP.".dbsnp";
      $arrayOfSNPStructures=retrieve($fileOutFreezeSNP);
      print "done\n";
    }
  }














  ############################################
  ##                                        ##
  ##                                        ##
  ##           BEGIN ANNOTATION             ##
  ##                                        ##
  ##                                        ##
  ############################################

  #files to print to
  #mode 1&2 annotation of coordinates and intervals
  # File name           1) Input or just genes ?         2) Filter for input ?             3) Filter for genes
  # input.all.out       input                            all inputs                        none (all the transcripts, all the genes)
  # input.single.out    input                            all inputs                        pick one transcript per gene but report all the different genes
  # input.best.out      input                            all inputs                        pick the most representative hit (in exon, in intron,

  # genes.all.out       genes                            all inputs                        none (all the transcripts, all the genes)
  # genes.single.out    genes                            all inputs                        pick one transcript per gene but report all the different genes
  # genes.best.out      genes                            all inputs                        pick the most representative hit (in exon, in intron,

  #                                                                                        upstream, downstream) in that order
  # input.none.out      input                            those with no hit                 (no gene will be reported)
  # input.exon.out      input                            All those in exon                 pick one transcript per gene but report all the different genes
  # input.intron.out    input                            All those in introns              pick one transcript per gene but report all the different genes
  # input.upstrm.out    input                            All those upstream of a gene      pick one transcript per gene but report all the different genes
  # input.dwnstrm.out   input                            All those downstream of a gene    pick one transcript per gene but report all the different genes

  #
  #mode 3 snp
  #
  # input.all.out       input                            all inputs                        none (all the transcripts, all the genes)
  # input.single.out    input                            all inputs                        pick one transcript per gene but report all the different genes
  # input.best.out      input                            all inputs                        pick the most representative hit (in exon, in intron,
  #                                                                                        upstream, downstream) in that order

  # genes.all.out       genes                            all inputs                        none (all the transcripts, all the genes)
  # genes.single.out    genes                            all inputs                        pick one transcript per gene but report all the different genes
  # genes.best.out      genes                            all inputs                        pick the most representative hit (in exon, in intron,
  #                                                                                        upstream, downstream) in that order



  # input.none.out      input                            those with no hit                 (no gene will be reported)
  # input.cde.out       input                            All those in coding exons         pick one transcript per gene but only genes with a coding exon hit
  # input.ncde.out      input                            All those in non-coding exons     pick one transcript per gene but only genes with a coding exon hit
  # input.nsyn.out      input                            Those causing non-syn mutations   pick one transcript per gene but only genes with non-syn mutation
  # input.syn.out       input                            Those causing synon mutations     pick one transcript per gene but only genes with syn mutation
  # input.intron.out    input                            All those in introns              pick one transcript per gene but report all the different genes
  # input.upstrm.out    input                            All those upstream of a gene      pick one transcript per gene but report all the different genes
  # input.dwnstrm.out   input                            All those downstream of a gene    pick one transcript per gene but report all the different genes

  #mode 4 closest TSS
  # input.closest.out   input                            all inputs on chr with genes      the closest transcript
  # input.none.out      input                            all inputs on chr without genes   (no gene will be reported)

  print "Beginning annotation ...";

  foreach my $fileName (@arrayOfFilesToAnnotate) { #For each file to annotate



    my $INPUTALLXML;

    ############################################
    ##                                        ##
    ##    BEGIN ANNOTATION MODE 1,2  and 3    ##
    ##                                        ##
    ############################################
    if ( $annotationMode == 1 || $annotationMode == 2 || $annotationMode == 3 ) { #coord,range,snp annotation mode

      foreach my $range (@arrayOfRanges) { #For each range to consider


	my $allStats=StatiticsModule->new(genomeUsed      => $speciesCode,
					  databaseUsed    => $databaseName,
					  databaseFiles   => $databaseFiles,
					  range           => $range,
					  mode            => $annotationMode);
	my $singleStats=StatiticsModule->new(genomeUsed      => $speciesCode,
					     databaseUsed    => $databaseName,
					     databaseFiles   => $databaseFiles,
					     range           => $range,
					     mode            => $annotationMode);
	my $bestStats=StatiticsModule->new(genomeUsed      => $speciesCode,
					   databaseUsed    => $databaseName,
					   databaseFiles   => $databaseFiles,
					   range           => $range,
					   mode            => $annotationMode);



	if ( $annotationMode == 1 || $annotationMode == 2 ) {
	  if (!$printStatsOnly) {
	    open(INPUTNONEOUT,    ">".$fileName.".".$range.".".$databaseName.".input.none.out");
	    open(INPUTALLOUT,     ">".$fileName.".".$range.".".$databaseName.".input.all.out");
	  }

	  if ($printFullFiles) {
	    open(INPUTALLSTAT,    ">".$fileName.".".$range.".".$databaseName.".input.all.stat");

	    if (!$customDatabase) {

	      if (!$doNotPrintXML) {
		if (!$printStatsOnly) {
		  open($INPUTALLXML,    ">".$fileName.".".$range.".".$databaseName.".input.all.xml");
		}
		open(INPUTALLSTATXML,    ">".$fileName.".".$range.".".$databaseName.".input.all.stat.xml");
	      }

	      if (!$singleTranscriptDatabase) {

		if (!$printStatsOnly) {
		  open(INPUTSINGLEOUT,  ">".$fileName.".".$range.".".$databaseName.".input.single.out");
		  open(INPUTBESTOUT,    ">".$fileName.".".$range.".".$databaseName.".input.best.out");
		}
		open(INPUTSINGLESTAT, ">".$fileName.".".$range.".".$databaseName.".input.single.stat");
		open(INPUTBESTSTAT,   ">".$fileName.".".$range.".".$databaseName.".input.best.stat");
		if (!$doNotPrintXML) {
		  open(INPUTSINGLESTATXML, ">".$fileName.".".$range.".".$databaseName.".input.single.stat.xml");
		  open(INPUTBESTSTATXML,   ">".$fileName.".".$range.".".$databaseName.".input.best.stat.xml");
		}

	      }else{ #single transcript

		if (!$printStatsOnly) {
		  open(INPUTBESTOUT,    ">".$fileName.".".$range.".".$databaseName.".input.best.out");
		}
		open(INPUTBESTSTAT,   ">".$fileName.".".$range.".".$databaseName.".input.best.stat");
		if (!$doNotPrintXML) {
		  open(INPUTBESTSTATXML,   ">".$fileName.".".$range.".".$databaseName.".input.best.stat.xml");
		}

	      }

	      if (!$printStatsOnly) {
		if (!$doNotPrintGenesOut) {
		  open(GENESALLOUT,     ">".$fileName.".".$range.".".$databaseName.".genes.all.out");

		  if (!$singleTranscriptDatabase) {
		    open(GENESSINGLEOUT,  ">".$fileName.".".$range.".".$databaseName.".genes.single.out");
		    open(GENESBESTOUT,    ">".$fileName.".".$range.".".$databaseName.".genes.best.out");
		  }else{
		    open(GENESBESTOUT,    ">".$fileName.".".$range.".".$databaseName.".genes.best.out");
		  }
		}


		open(INPUTEXONOUT,    ">".$fileName.".".$range.".".$databaseName.".input.exon.out");
		open(INPUTINTRONOUT,  ">".$fileName.".".$range.".".$databaseName.".input.intron.out");
		open(INPUTUPSTRMOUT,  ">".$fileName.".".$range.".".$databaseName.".input.upstrm.out");
		open(INPUTDWNSTRMOUT, ">".$fileName.".".$range.".".$databaseName.".input.dwnstrm.out");
	      }
	    }else{
	      if ($printFullFiles) {
		$allStats->useOfCustomDatabase();
		$singleStats->useOfCustomDatabase();
		$bestStats->useOfCustomDatabase();
	      }
	    }
	  }
	} elsif ($annotationMode == 3 ) {

	  if ($compareWithDBSNP) {
	    if (!$printStatsOnly) {
	      open(INPUTDBSNPOUT,     ">".$fileName.".".$range.".".$databaseName.".input.dbsnp.out");
	    }
	  }
	  if (!$printStatsOnly) {
	    open(INPUTALLOUT,     ">".$fileName.".".$range.".".$databaseName.".input.all.out");
	  }
	  if ($printFullFiles) {
	    open(INPUTALLSTAT,    ">".$fileName.".".$range.".".$databaseName.".input.all.stat");

	    if (!$printStatsOnly) {
	      open(INPUTSINGLEOUT,  ">".$fileName.".".$range.".".$databaseName.".input.single.out");
	      open(INPUTBESTOUT,    ">".$fileName.".".$range.".".$databaseName.".input.best.out");

	      open(INPUTSINGLESTAT, ">".$fileName.".".$range.".".$databaseName.".input.single.stat");
	      open(INPUTBESTSTAT,   ">".$fileName.".".$range.".".$databaseName.".input.best.stat");

	      if (!$doNotPrintXML) {
		open($INPUTALLXML,    ">".$fileName.".".$range.".".$databaseName.".input.all.xml");
		open(INPUTALLSTATXML, ">".$fileName.".".$range.".".$databaseName.".input.all.stat.xml");
		open(INPUTSINGLESTATXML, ">".$fileName.".".$range.".".$databaseName.".input.single.stat.xml");
		open(INPUTBESTSTATXML,   ">".$fileName.".".$range.".".$databaseName.".input.best.stat.xml");
	      }

	      if (!$doNotPrintGenesOut) {
		open(GENESALLOUT,     ">".$fileName.".".$range.".".$databaseName.".genes.all.out");
		open(GENESSINGLEOUT,  ">".$fileName.".".$range.".".$databaseName.".genes.single.out");
		open(GENESBESTOUT,    ">".$fileName.".".$range.".".$databaseName.".genes.best.out");
	      }

	      open(INPUTCDEOUT,     ">".$fileName.".".$range.".".$databaseName.".input.cde.out");
	      open(INPUTNCDEOUT,    ">".$fileName.".".$range.".".$databaseName.".input.ncde.out");
	      open(INPUTNSYNOUT,    ">".$fileName.".".$range.".".$databaseName.".input.nsyn.out");
	      open(INPUTSYNOUT,     ">".$fileName.".".$range.".".$databaseName.".input.syn.out");
	      open(INPUTEXONOUT,    ">".$fileName.".".$range.".".$databaseName.".input.exon.out");
	      open(INPUTINTRONOUT,  ">".$fileName.".".$range.".".$databaseName.".input.intron.out");
	      open(INPUTUPSTRMOUT,  ">".$fileName.".".$range.".".$databaseName.".input.upstrm.out");
	      open(INPUTDWNSTRMOUT, ">".$fileName.".".$range.".".$databaseName.".input.dwnstrm.out");
	    }

	  }
	  if (!$printStatsOnly) {
	    open(INPUTNONEOUT,    ">".$fileName.".".$range.".".$databaseName.".input.none.out");
	  }
	} else {
	  print STDERR "invalid mode\n";
	  terminate(9);
	}

	my $resultOpen=open(INPUTFILE,$fileName) or die "Unable to open $fileName\n";

	my %hashOfIdInput;


	my @arrayOfFieldHeader;
	if ( $annotationMode == 1) { #coord annotation mode

	  if (!$customDatabase) {
	    @arrayOfFieldHeader=("#INPUTID","CHR","COORD","TRANSCRIPT","GENE","POSITION","INDEX","DISTANCE5p","DISTANCE3p","PARTIAL");
	  } else {
	    @arrayOfFieldHeader=("#INPUTID","CHR","COORD","DBID","POSITION","DISTANCE5p","DISTANCE3p");
	  }

	} elsif ( $annotationMode == 2) { #range annotation mode

	  if (!$customDatabase) {
	    @arrayOfFieldHeader=("#INPUTID","CHR","COORD1","COORD2","TRANSCRIPT","GENE","POSITION","INDEXexons","INDEXintrons","PARTIAL");
	  } else {
	    @arrayOfFieldHeader=("#INPUTID","CHR","COORD1","COORD2","DBID","POSITION");
	  }

	} elsif ( $annotationMode == 3) { #coord annotation mode
	  @arrayOfFieldHeader=("#INPUTID","CHR","COORD","TRANSCRIPT","GENE","POSITION","INDEX","DISTANCE5p","DISTANCE3p","PARTIAL","CODONREF","CODONREAD","AAREF","AAREAD");
	  if ($produceAASeq) {
	    push(@arrayOfFieldHeader,"INDEXAA");
	    push(@arrayOfFieldHeader,"SEQREF");
	    push(@arrayOfFieldHeader,"SEQREAD");
	  }

	  if ($compareWithDBSNP) {
	    if (!$printStatsOnly) {
	      print INPUTDBSNPOUT    join("\t",( $arrayOfFieldHeader[0],$arrayOfFieldHeader[1],$arrayOfFieldHeader[2],"DBSNP" ))."\n";
	    }
	  }

	} else {
	  print STDERR "Wrong mode for field header\n";
	  terminate(9);
	}
	if (!$printStatsOnly) {
	  print INPUTALLOUT     join("\t",@arrayOfFieldHeader)."\n";
	  if ($printFullFiles) {
	    if (!$customDatabase) {
	      if (!$singleTranscriptDatabase) {
		print INPUTSINGLEOUT  join("\t",@arrayOfFieldHeader)."\n";
		print INPUTBESTOUT    join("\t",@arrayOfFieldHeader)."\n";
	      }else{
		print INPUTBESTOUT    join("\t",@arrayOfFieldHeader)."\n";
	      }

	      print INPUTNONEOUT    join("\t",( $arrayOfFieldHeader[0],$arrayOfFieldHeader[1],$arrayOfFieldHeader[2] ))."\n";
	      print INPUTEXONOUT    join("\t",@arrayOfFieldHeader)."\n";
	      print INPUTINTRONOUT  join("\t",@arrayOfFieldHeader)."\n";
	      print INPUTUPSTRMOUT  join("\t",@arrayOfFieldHeader)."\n";
	      print INPUTDWNSTRMOUT join("\t",@arrayOfFieldHeader)."\n";
	    }
	  }
	}
	my %hashOfGenesAll;
	my %hashOfGenesSingle;
	my %hashOfGenesBest;

	my $blankLines=0;


	while (my $line =<INPUTFILE>) { #for each line
	  chomp($line);

	  if (length($line) == 0) {
	    $blankLines++;
	    next;
	  } else {
	    if ($blankLines > 0) {
	      print STDERR "Blank lines found prior to line $line in $fileName \n";
	      terminate(9);
	    }
	  }




	  if ( $annotationMode == 1) { #coord annotation mode

	    # First input format, for coordinate annotation
	    #chr15	67532214	myCoordinate1293
	    my $chrInput     ;
	    my $coordInput  ;
	    my $idInput      ;


	    if ($line =~ /^(\S+)\s+(\d+)\s+(\S+)\s*$/) {
	      $chrInput     = $1;
	      $coordInput   = $2;
	      $idInput      = $3;
	    } else {		#end for each line
	      print STDERR "Line $line in $fileName did not parse\n";
	      terminate(9);
	    }



	    if (exists $hashOfIdInput{$idInput}) {
	      print STDERR "The input id ".$idInput." was found twice\n";
	      terminate(9);
	    } else {
	      $hashOfIdInput{$idInput}=0;
	    }

	    if ( !ChromosomeIndexing::chrEXISTS($speciesCode,$chrInput) ) {
	      print STDERR "Chromosome in file $fileName contains chromosome $chrInput which is not recognized for the current build\n";
	      terminate(9);
	    }

	    #put if
	    if (!$customDatabase) {
	      if ($printFullFiles) {
		$allStats->addSite();
		$singleStats->addSite();
		$bestStats->addSite();
	      }
	      #find coordinate in tree
	      my $array=[];


	      if ($range == 0) {
		if ($arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$chrInput)]) {
		  my $indexbinst=int($coordInput/$genomebinsize);
		  my $indexbinen=int($coordInput/$genomebinsize);

		  for (my $indexbin=$indexbinst;$indexbin<=$indexbinen;$indexbin++) {
		    if ($arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$chrInput)]->[$indexbin]) {
		      push( @{$array} , @{$arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$chrInput)]->[$indexbin]->seekCoord($coordInput)} );
		    }
		  }

		}
	      } else {
		if ($arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$chrInput)]) {
		  my $indexbinst=int(($coordInput-$range)/$genomebinsize);
		  my $indexbinen=int(($coordInput+$range)/$genomebinsize);

		  for (my $indexbin=$indexbinst;$indexbin<=$indexbinen;$indexbin++) {
		    if ($arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$chrInput)]->[$indexbin]) {
		      push( @{$array},  @{$arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$chrInput)]->[$indexbin]->seekRange($coordInput-$range,$coordInput+$range)} );
		    }
		  }
		}
	      }



	      if ($#{$array} >= 0 ) {
		$array=eliminateDuplicateTranscripts($array);


		my %hashOfGenes;

		for (my $indexArray=0;$indexArray<=$#{$array};$indexArray++   ) {
		  if ($printFullFiles) {
		    if (exists $hashOfGenes{$array->[$indexArray]->{'transcript'}->getGeneCounter()} ) {
		      push( @{$hashOfGenes{$array->[$indexArray]->{'transcript'}->getGeneCounter()}} , $array->[$indexArray]);
		    } else {
		      $hashOfGenes{$array->[$indexArray]->{'transcript'}->getGeneCounter()}=[$array->[$indexArray]];
		    }

		    my $hashForGeneFile ;

		    if ($singleTranscriptDatabase) {
		      $hashForGeneFile = {'idInput' => $idInput,
					  'trans'   => $array->[$indexArray]->{'transcript'}->getId()};
		    } else {
		      $hashForGeneFile = {'idInput' => $idInput,
					  'symbol'  => $array->[$indexArray]->{'transcript'}->getGeneSymbol(),
					  'trans'   => $array->[$indexArray]->{'transcript'}->getId()};
		    }

		    if (!$doNotPrintGenesOut) {
		      if (exists $hashOfGenesAll{$array->[$indexArray]->{'transcript'}->getGeneCounter()} ) {
			push( @{$hashOfGenesAll{$array->[$indexArray]->{'transcript'}->getGeneCounter()}} , $hashForGeneFile); #$array->[$indexArray] );
		      } else {
			$hashOfGenesAll{$array->[$indexArray]->{'transcript'}->getGeneCounter()}= [$hashForGeneFile]; #[$array->[$indexArray]];
		      }
		    }
		  }
		  $array->[$indexArray]->{'annoResult'} = $array->[$indexArray]->{'transcript'}->annotateCoordinate($coordInput);

		  my $stringToPrint;
		  $stringToPrint=join("\t",($idInput,$chrInput,$coordInput))."\t".$array->[$indexArray]->{'transcript'}->getId()."\t";

		  if (!$singleTranscriptDatabase) { #not single database
		    $stringToPrint.=$array->[$indexArray]->{'transcript'}->getGeneSymbol()."\t";
		  }

		  $stringToPrint.=formatAnnotation(1,$array->[$indexArray]->{'annoResult'} )."\n";

		  if (!$printStatsOnly) {
		    print INPUTALLOUT $stringToPrint;
		  }
		  if ($printFullFiles) {
		    $allStats->addAnnotation($array->[$indexArray]->{'annoResult'});
		    if (     $array->[$indexArray]->{'annoResult'}->{'code'} == 1) {
		      if (!$printStatsOnly) {
			print INPUTEXONOUT    $stringToPrint;
		      }
		    } elsif ( $array->[$indexArray]->{'annoResult'}->{'code'} == 2) {
		      if (!$printStatsOnly) {
			print INPUTINTRONOUT  $stringToPrint;
		      }
		    } elsif ( $array->[$indexArray]->{'annoResult'}->{'code'} == 3) {
		      if (!$printStatsOnly) {
			print INPUTUPSTRMOUT  $stringToPrint;
		      }
		    } elsif ( $array->[$indexArray]->{'annoResult'}->{'code'} == 4) {
		      if (!$printStatsOnly) {
			print INPUTDWNSTRMOUT $stringToPrint;
		      }
		    } else {
		      print STDERR "Wrong code\n";
		      terminate(9);
		    }

		  }
		} #for each result


		if ($printFullFiles) {
		  my @allGenes;
		  my %singleGenes;
		  my $bestResult;

		  if (!$singleTranscriptDatabase) { #not single database

		    foreach my $geneCounter (keys %hashOfGenes) {
		      my $singleBestResult=selectSingleBestTranscript( @{$hashOfGenes{$geneCounter}} );
		      push(@allGenes,@{$hashOfGenes{$geneCounter}});

		      if (exists $singleGenes{$geneCounter} ) {
			print STDERR "Duplicate gene counters\n";
			terminate(9);
		      } else {
			$singleGenes{$geneCounter} = $singleBestResult->{'transcript'}->getId() ;
		      }

		      $singleStats->addAnnotation($singleBestResult->{'annoResult'});

		      if (!$printStatsOnly) {
			print INPUTSINGLEOUT join("\t",($idInput,$chrInput,$coordInput))."\t".$singleBestResult->{'transcript'}->getId()."\t".$singleBestResult->{'transcript'}->getGeneSymbol()."\t".formatAnnotation(1,$singleBestResult->{'annoResult'})."\n";
		      }
		      my $hashForGeneFile = {'idInput' => $idInput,
					     'symbol'  => $singleBestResult->{'transcript'}->getGeneSymbol() ,
					     'trans'   => $singleBestResult->{'transcript'}->getId() };
		      if (!$doNotPrintGenesOut) {
			if (exists $hashOfGenesSingle{$singleBestResult->{'transcript'}->getGeneCounter()} ) {
			  push( @{$hashOfGenesSingle{$singleBestResult->{'transcript'}->getGeneCounter()}} , $hashForGeneFile); #$singleBestResult );
			} else {
			  $hashOfGenesSingle{$singleBestResult->{'transcript'}->getGeneCounter()}= [$hashForGeneFile]; #[$singleBestResult];
			}
		      }
		    }		#for each gene counter

		    $bestResult=selectSingleBestTranscript( @allGenes );
		    $bestStats->addAnnotation($bestResult->{'annoResult'});

		    if (!$printStatsOnly) {
		      print INPUTBESTOUT join("\t",($idInput,$chrInput,$coordInput))."\t".$bestResult->{'transcript'}->getId()."\t".$bestResult->{'transcript'}->getGeneSymbol()."\t".formatAnnotation(1, $bestResult->{'annoResult'} )."\n";
		    }
		    my $hashForGeneFile = {'idInput' => $idInput,
					   'symbol'  => $bestResult->{'transcript'}->getGeneSymbol(),
					   'trans'   => $bestResult->{'transcript'}->getId()  };

		    if (!$doNotPrintGenesOut) {
		      if (exists $hashOfGenesBest{$bestResult->{'transcript'}->getGeneCounter()} ) {
			push( @{$hashOfGenesBest{$bestResult->{'transcript'}->getGeneCounter()}} , $hashForGeneFile); #$bestResult );
		      } else {
			$hashOfGenesBest{$bestResult->{'transcript'}->getGeneCounter()}= [$hashForGeneFile]; #[$bestResult];
		      }
		    }

		  } else {	# single database
		    foreach my $geneCounter (keys %hashOfGenes) {
		      push(@allGenes,@{$hashOfGenes{$geneCounter}});
		    }

		    $bestResult=selectSingleBestTranscript( @allGenes );
		    $bestStats->addAnnotation($bestResult->{'annoResult'});

		    if (!$printStatsOnly) {
		      print INPUTBESTOUT join("\t",($idInput,$chrInput,$coordInput))."\t".$bestResult->{'transcript'}->getId()."\t".formatAnnotation(1, $bestResult->{'annoResult'} )."\n";
		    }
		    my $hashForGeneFile = {'idInput' => $idInput,
					   'trans'   => $bestResult->{'transcript'}->getId()  };

		    if (!$doNotPrintGenesOut) {
		      if (exists $hashOfGenesBest{$bestResult->{'transcript'}->getTranscriptCounter()} ) {
			push( @{$hashOfGenesBest{$bestResult->{'transcript'}->getTranscriptCounter()}} , $hashForGeneFile); #$bestResult );
		      } else {
			$hashOfGenesBest{$bestResult->{'transcript'}->getTranscriptCounter()}= [$hashForGeneFile]; #[$bestResult];
		      }
		    }

		  }




		  if (!$doNotPrintXML) {
		    if (!$printStatsOnly) {
		      printXML(0,$annotationMode,$INPUTALLXML,$coordInput,$chrInput,$idInput,$bestResult,-1,\%singleGenes,\%hashOfGenes,\@arrayOfFieldHeader,$singleTranscriptDatabase);
		    }
		  }
		}
	      } else {
		if ($printFullFiles) {
		  $allStats->addToNoHits();
		  $bestStats->addToNoHits();
		  $singleStats->addToNoHits();
		}
		if (!$printStatsOnly) {
		  print INPUTNONEOUT join("\t",($idInput,$chrInput,$coordInput))."\n";
		}
	      }

	    } else {		#for custom databases

	      $allStats->addSite();

	      #find coordinate in tree
	      my $array=[];



	      if ($range == 0) {
		if ($arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$chrInput)]) {
		  my $indexbinst=int($coordInput/$genomebinsize);
		  my $indexbinen=int($coordInput/$genomebinsize);

		  for (my $indexbin=$indexbinst;$indexbin<=$indexbinen;$indexbin++) {
		    if ($arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$chrInput)]->[$indexbin]) {
		      push( @{$array} , @{$arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$chrInput)]->[$indexbin]->seekCoord($coordInput)} );
		    }
		  }

		}
	      } else {
		if ($arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$chrInput)]) {
		  my $indexbinst=int(($coordInput-$range)/$genomebinsize);
		  my $indexbinen=int(($coordInput+$range)/$genomebinsize);

		  for (my $indexbin=$indexbinst;$indexbin<=$indexbinen;$indexbin++) {
		    if ($arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$chrInput)]->[$indexbin]) {
		      push( @{$array},  @{$arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$chrInput)]->[$indexbin]->seekRange($coordInput-$range,$coordInput+$range)} );
		    }
		  }
		}
	      }

	      #eliminate duplicate


	      if ($#{$array} >= 0 ) {
		$array=eliminateDuplicateRecords($array);

		for (my $indexArray=0;$indexArray<=$#{$array};$indexArray++ ) {
		  my $stringToPrint;
		  $stringToPrint=join("\t",($idInput,$chrInput,$coordInput));
		  $stringToPrint.="\t".$array->[$indexArray]->{'record'}->{'id'};

		  if (    $array->[$indexArray]->{'record'}->{'strand'} eq "+") {

		    if ( ($coordInput < $array->[$indexArray]->{'record'}->{'range'}->getStart() &&
			  $coordInput < $array->[$indexArray]->{'record'}->{'range'}->getEnd() ) ) {

		      $stringToPrint.=
			"\t"."upstream".
			  "\t".($array->[$indexArray]->{'record'}->{'range'}->getStart()-$coordInput).
			    "\t".($array->[$indexArray]->{'record'}->{'range'}->getEnd()-$coordInput);
		      if ($printFullFiles) {
			$allStats->addCustomUpstream();
		      }

		    } elsif ( ($coordInput > $array->[$indexArray]->{'record'}->{'range'}->getStart() &&
			       $coordInput > $array->[$indexArray]->{'record'}->{'range'}->getEnd() ) ) {

		      $stringToPrint.=
			"\t"."downstream".
			  "\t".($coordInput-$array->[$indexArray]->{'record'}->{'range'}->getStart()).
			    "\t".($coordInput-$array->[$indexArray]->{'record'}->{'range'}->getEnd());
		      if ($printFullFiles) {
			$allStats->addCustomDownstream();
		      }

		    } elsif ( $array->[$indexArray]->{'record'}->{'range'}->contains($coordInput) ) {

		      $stringToPrint.=
			"\t"."contained".
			  "\t".($coordInput-$array->[$indexArray]->{'record'}->{'range'}->getStart()).
			    "\t".($array->[$indexArray]->{'record'}->{'range'}->getEnd()-$coordInput);
		      if ($printFullFiles) {
			$allStats->addCustomContained();
		      }

		    } else {
		      print STDERR "Invalid state dump: ".Dumper($array->[$indexArray])."\n";
		      terminate(9);
		    }

		  } elsif ($array->[$indexArray]->{'record'}->{'strand'} eq "-") {

		    if ( ($coordInput < $array->[$indexArray]->{'record'}->{'range'}->getStart() &&
			  $coordInput < $array->[$indexArray]->{'record'}->{'range'}->getEnd() ) ) {

		      $stringToPrint.=
			"\t"."downstream".
			  "\t".($array->[$indexArray]->{'record'}->{'range'}->getEnd()-$coordInput).
			    "\t".($array->[$indexArray]->{'record'}->{'range'}->getStart()-$coordInput);
		      if ($printFullFiles) {
			$allStats->addCustomDownstream();
		      }

		    } elsif ( ($coordInput > $array->[$indexArray]->{'record'}->{'range'}->getStart() &&
			       $coordInput > $array->[$indexArray]->{'record'}->{'range'}->getEnd() ) ) {

		      $stringToPrint.=
			"\t"."upstream".
			  "\t".($coordInput - $array->[$indexArray]->{'record'}->{'range'}->getEnd()).
			    "\t".($coordInput - $array->[$indexArray]->{'record'}->{'range'}->getStart());
		      if ($printFullFiles) {
			$allStats->addCustomUpstream();
		      }

		    } elsif ( $array->[$indexArray]->{'record'}->{'range'}->contains($coordInput) ) {

		      $stringToPrint.=
			"\t"."contained".
			  "\t".($array->[$indexArray]->{'record'}->{'range'}->getEnd()-$coordInput).
			    "\t".($coordInput-$array->[$indexArray]->{'record'}->{'range'}->getStart());
		      if ($printFullFiles) {
			$allStats->addCustomContained();
		      }

		    } else {
		      print STDERR "Invalid state dump: ".Dumper($array->[$indexArray])."\n";
		      terminate(9);
		    }
		  } else {
		    print STDERR "Problem with strand in tree structure dump: ".Dumper($array->[$indexArray])."\n";
		    terminate(9);
		  }

		  if (!$printStatsOnly) {
		    print INPUTALLOUT $stringToPrint."\n";
		  }
		}
	      } else {
		#add site to stats
		if ($printFullFiles) {
		  $allStats->addToNoHits();
		}
		if (!$printStatsOnly) {
		  print INPUTNONEOUT join("\t",($idInput,$chrInput,$coordInput))."\n";
		}
	      }

	    }
	    #end if









	  } elsif ( $annotationMode == 2) { #range annotation mode

	    my $chrInput     ;
	    my $coordInput1  ;
	    my $coordInput2  ;
	    my $idInput      ;

	    if ($line =~ /^(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s*$/) {
	      $chrInput     = $1;
	      $coordInput1  = $2;
	      $coordInput2  = $3;
	      $idInput      = $4;

	      if ($coordInput1 > $coordInput2) {
		print STDERR "Line $line in $fileName the first coordinate is greater than the second\n";
		terminate(9);
	      }

	      if ( ($coordInput2 - $coordInput1) > $MAXBPFORINTERVALS ) {
		print STDERR "Line $line in $fileName the second coordinate cannot be more than $MAXBPFORINTERVALS from the first\n";
		terminate(9);
	      }

	    } else {		#end for each line
	      print STDERR "Line $line in $fileName did not parse\n";
	      terminate(9);
	    }


	    if (exists $hashOfIdInput{$idInput}) {
	      print STDERR "The input id ".$idInput." was found twice\n";
	      terminate(9);
	    } else {
	      $hashOfIdInput{$idInput}=0;
	    }

	    if ( !ChromosomeIndexing::chrEXISTS($speciesCode,$chrInput) ) {
	      print STDERR "Chromosome in file $fileName contains chromosome $chrInput which is not recognized for the current build\n";
	      terminate(9);
	    }

	    if (!$customDatabase) {
	      if ($printFullFiles) {
		$allStats->addSite();
		$singleStats->addSite();
		$bestStats->addSite();
	      }

	      #find coordinate in tree
	      my $array=[];


	      if ($arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$chrInput)]) {
		my $indexbinst=int(($coordInput1-$range)/$genomebinsize);
		my $indexbinen=int(($coordInput2+$range)/$genomebinsize);

		for (my $indexbin=$indexbinst;$indexbin<=$indexbinen;$indexbin++) {
		  if ($arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$chrInput)]->[$indexbin]) {
		    push( @{$array},  @{$arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$chrInput)]->[$indexbin]->seekRange($coordInput1-$range,$coordInput2+$range)} );
		  }
		}
	      }



	      if ($#{$array} >= 0 ) {
		$array=eliminateDuplicateTranscripts($array);

		my %hashOfGenes;

		for (my $indexArray=0;$indexArray<=$#{$array};$indexArray++   ) {
		  if ($printFullFiles) {
		    if (exists $hashOfGenes{$array->[$indexArray]->{'transcript'}->getGeneCounter()} ) {
		      push( @{$hashOfGenes{$array->[$indexArray]->{'transcript'}->getGeneCounter()}} , $array->[$indexArray]);
		    } else {
		      $hashOfGenes{$array->[$indexArray]->{'transcript'}->getGeneCounter()}=[$array->[$indexArray]];
		    }

		    my $hashForGeneFile ;

		    if ($singleTranscriptDatabase) {
		      $hashForGeneFile = {'idInput' => $idInput,
					  'trans'   => $array->[$indexArray]->{'transcript'}->getId()};
		    } else {
		      $hashForGeneFile = {'idInput' => $idInput,
					  'symbol'  => $array->[$indexArray]->{'transcript'}->getGeneSymbol(),
					  'trans'   => $array->[$indexArray]->{'transcript'}->getId()};
		    }

		    if (!$doNotPrintGenesOut) {
		      if (exists $hashOfGenesAll{$array->[$indexArray]->{'transcript'}->getGeneCounter()} ) {
			push( @{$hashOfGenesAll{$array->[$indexArray]->{'transcript'}->getGeneCounter()}} , $hashForGeneFile); #$array->[$indexArray] );
		      } else {
			$hashOfGenesAll{$array->[$indexArray]->{'transcript'}->getGeneCounter()}= [$hashForGeneFile]; #[$array->[$indexArray]];
		      }
		    }
		  }

		  $array->[$indexArray]->{'annoResult'} = $array->[$indexArray]->{'transcript'}->annotateRange($coordInput1,$coordInput2);

		  my $stringToPrint;
		  $stringToPrint=join("\t",($idInput,$chrInput,$coordInput1,$coordInput2))."\t".$array->[$indexArray]->{'transcript'}->getId()."\t";

		  if (!$singleTranscriptDatabase) {
		    $stringToPrint.=$array->[$indexArray]->{'transcript'}->getGeneSymbol()."\t";
		  }

		  $stringToPrint.=formatIntervalAnnotation($array->[$indexArray]->{'annoResult'} )."\n";


		  if (!$printStatsOnly) {
		    print INPUTALLOUT $stringToPrint;
		  }
		  if ($printFullFiles) {
		    $allStats->addAnnotationInterval($array->[$indexArray]->{'annoResult'});

		    if (      $array->[$indexArray]->{'annoResult'}->{'code'} == 1) {
		      if (!$printStatsOnly) {
			print INPUTEXONOUT    $stringToPrint;
		      }

		    } elsif ( $array->[$indexArray]->{'annoResult'}->{'code'} == 2) {
		      if (!$printStatsOnly) {
			print INPUTEXONOUT    $stringToPrint;
			print INPUTINTRONOUT  $stringToPrint; 
		      }
		    } elsif ( $array->[$indexArray]->{'annoResult'}->{'code'} == 3) {
		      if (!$printStatsOnly) {
			print INPUTINTRONOUT  $stringToPrint;
		      }

		    } elsif ( $array->[$indexArray]->{'annoResult'}->{'code'} == 4) {

		    } elsif ( $array->[$indexArray]->{'annoResult'}->{'code'} == 5) {
		      if (!$printStatsOnly) {
			print INPUTEXONOUT    $stringToPrint;
			print INPUTUPSTRMOUT  $stringToPrint;
		      }
		    } elsif ( $array->[$indexArray]->{'annoResult'}->{'code'} == 6) {
		      if (!$printStatsOnly) {
			print INPUTUPSTRMOUT  $stringToPrint;
			print INPUTEXONOUT    $stringToPrint;
			print INPUTINTRONOUT  $stringToPrint;
		      }
		    } elsif ( $array->[$indexArray]->{'annoResult'}->{'code'} == 7) {
		      if (!$printStatsOnly) {
			print INPUTDWNSTRMOUT $stringToPrint;
			print INPUTEXONOUT    $stringToPrint;
		      }
		    } elsif ( $array->[$indexArray]->{'annoResult'}->{'code'} == 8) {
		      if (!$printStatsOnly) {
			print INPUTDWNSTRMOUT $stringToPrint;
			print INPUTEXONOUT    $stringToPrint;
			print INPUTINTRONOUT  $stringToPrint;
		      }
		    } elsif ( $array->[$indexArray]->{'annoResult'}->{'code'} == 9) {
		      if (!$printStatsOnly) {
			print INPUTUPSTRMOUT  $stringToPrint;
		      }
		    } elsif ( $array->[$indexArray]->{'annoResult'}->{'code'} == 10) {
		      if (!$printStatsOnly) {
			print INPUTDWNSTRMOUT $stringToPrint;
		      }
		    } else {
		      print STDERR "Wrong code\n";
		      terminate(9);
		    }
		  }
		}


		if ($printFullFiles) {
		  my @allGenes;
		  my %singleGenes;
		  my $bestResult;

		  if (!$singleTranscriptDatabase) {
		    foreach my $geneCounter (keys %hashOfGenes) {
		      my $singleBestResult=selectSingleBestTranscript( @{$hashOfGenes{$geneCounter}} );
		      push(@allGenes,@{$hashOfGenes{$geneCounter}});


		      if (exists $singleGenes{$geneCounter} ) {
			print STDERR "Duplicate gene counters\n";
			terminate(9);
		      } else {
			$singleGenes{$geneCounter} = $singleBestResult->{'transcript'}->getId() ;
		      }

		      $singleStats->addAnnotationInterval($singleBestResult->{'annoResult'});


		      if (!$printStatsOnly) {
			print INPUTSINGLEOUT join("\t",($idInput,$chrInput,$coordInput1))."\t".$singleBestResult->{'transcript'}->getId()."\t".$singleBestResult->{'transcript'}->getGeneSymbol()."\t".formatIntervalAnnotation($singleBestResult->{'annoResult'})."\n";
		      }
		      my $hashForGeneFile = {'idInput' => $idInput,
					     'symbol'  => $singleBestResult->{'transcript'}->getGeneSymbol() ,
					     'trans'   => $singleBestResult->{'transcript'}->getId() };

		      if (!$doNotPrintGenesOut) {
			if (exists $hashOfGenesSingle{$singleBestResult->{'transcript'}->getGeneCounter()} ) {
			  push( @{$hashOfGenesSingle{$singleBestResult->{'transcript'}->getGeneCounter()}} , $hashForGeneFile); #$singleBestResult );
			} else {
			  $hashOfGenesSingle{$singleBestResult->{'transcript'}->getGeneCounter()}= [$hashForGeneFile]; #[$singleBestResult];
			}
		      }
		    }

		    $bestResult=selectSingleBestTranscript( @allGenes );
		    $bestStats->addAnnotationInterval($bestResult->{'annoResult'});

		    if (!$printStatsOnly) {
		      print INPUTBESTOUT join("\t",($idInput,$chrInput,$coordInput1))."\t".$bestResult->{'transcript'}->getId()."\t".$bestResult->{'transcript'}->getGeneSymbol()."\t".formatIntervalAnnotation( $bestResult->{'annoResult'} )."\n";
		    }
		    my $hashForGeneFile = {'idInput' => $idInput,
					   'symbol'  => $bestResult->{'transcript'}->getGeneSymbol(),
					   'trans'   => $bestResult->{'transcript'}->getId()  };

		    if (!$doNotPrintGenesOut) {
		      if (exists $hashOfGenesBest{$bestResult->{'transcript'}->getGeneCounter()} ) {
			push( @{$hashOfGenesBest{$bestResult->{'transcript'}->getGeneCounter()}} , $hashForGeneFile); #$bestResult );
		      } else {
			$hashOfGenesBest{$bestResult->{'transcript'}->getGeneCounter()}= [$hashForGeneFile]; #[$bestResult];
		      }
		    }

		  }else{ #for single transcripts

		    foreach my $geneCounter (keys %hashOfGenes) {
		      push(@allGenes,@{$hashOfGenes{$geneCounter}});
		    }

		    $bestResult=selectSingleBestTranscript( @allGenes );
		    $bestStats->addAnnotationInterval($bestResult->{'annoResult'});

		    if (!$printStatsOnly) {
		      print INPUTBESTOUT join("\t",($idInput,$chrInput,$coordInput1))."\t".$bestResult->{'transcript'}->getId()."\t".formatIntervalAnnotation( $bestResult->{'annoResult'} )."\n";
		    }
		    my $hashForGeneFile = {'idInput' => $idInput,
					   'trans'   => $bestResult->{'transcript'}->getId()  };

		    if (!$doNotPrintGenesOut) {
		      if (exists $hashOfGenesBest{$bestResult->{'transcript'}->getTranscriptCounter()} ) {
			push( @{$hashOfGenesBest{$bestResult->{'transcript'}->getTranscriptCounter()}} , $hashForGeneFile); #$bestResult );
		      } else {
			$hashOfGenesBest{$bestResult->{'transcript'}->getTranscriptCounter()}= [$hashForGeneFile]; #[$bestResult];
		      }
		    }

		  }



		  if (!$doNotPrintXML) {
		    if (!$printStatsOnly) {
		      printXML(0,$annotationMode,$INPUTALLXML,$coordInput1,$chrInput,$idInput,$bestResult,-1,\%singleGenes,\%hashOfGenes,\@arrayOfFieldHeader,$singleTranscriptDatabase);
		    }
		  }

		}
	      } else {
		if ($printFullFiles) {
		  $allStats->addToNoHits();
		  $bestStats->addToNoHits();
		  $singleStats->addToNoHits();
		}
		if (!$printStatsOnly) {
		  print INPUTNONEOUT join("\t",($idInput,$chrInput,$coordInput1))."\n";
		}
	      }
	    } else {

	      $allStats->addSite();

	      #find coordinate in tree
	      my $array=[];

	      if ($arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$chrInput)]) {
		my $indexbinst=int(($coordInput1-$range)/$genomebinsize);
		my $indexbinen=int(($coordInput2+$range)/$genomebinsize);

		for (my $indexbin=$indexbinst;$indexbin<=$indexbinen;$indexbin++) {
		  if ($arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$chrInput)]->[$indexbin]) {
		    push( @{$array},  @{$arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$chrInput)]->[$indexbin]->seekRange($coordInput1-$range,$coordInput2+$range)} );
		  }
		}
	      }

	      #eliminate duplicate


	      if ($#{$array} >= 0 ) {
		$array=eliminateDuplicateRecords($array);

		for (my $indexArray=0;$indexArray<=$#{$array};$indexArray++ ) {
		  my $stringToPrint;
		  $stringToPrint=join("\t",($idInput,$chrInput,$coordInput1,$coordInput2));
		  $stringToPrint.="\t".$array->[$indexArray]->{'record'}->{'id'};

		  if (    $array->[$indexArray]->{'record'}->{'strand'} eq "+") {
		    if ( ($coordInput1 < $array->[$indexArray]->{'record'}->{'range'}->getStart() &&
			  $coordInput2 < $array->[$indexArray]->{'record'}->{'range'}->getStart() ) ) {
		      $stringToPrint.="\t"."upstream";
		      if ($printFullFiles) {
			$allStats->addCustomUpstream();
		      }
		    } elsif ( ($coordInput1 > $array->[$indexArray]->{'record'}->{'range'}->getEnd() &&
			       $coordInput2 > $array->[$indexArray]->{'record'}->{'range'}->getEnd() ) ) {
		      $stringToPrint.="\t"."downstream";
		      if ($printFullFiles) {
			$allStats->addCustomDownstream();
		      }
		    } elsif ( $array->[$indexArray]->{'record'}->{'range'}->contains($coordInput1) &&
			      $array->[$indexArray]->{'record'}->{'range'}->contains($coordInput2) ) {
		      $stringToPrint.="\t"."contained";
		      if ($printFullFiles) {
			$allStats->addCustomContained();
		      }
		    } elsif ( $coordInput1 < $array->[$indexArray]->{'record'}->{'range'}->getStart() &&
			      $array->[$indexArray]->{'record'}->{'range'}->contains($coordInput2) ) {
		      $stringToPrint.="\t"."overlaps5P";
		      if ($printFullFiles) {
			$allStats->addCustomOverlaps5p();
		      }
		    } elsif ( $array->[$indexArray]->{'record'}->{'range'}->contains($coordInput1) &&
			      $coordInput2 > $array->[$indexArray]->{'record'}->{'range'}->getEnd()  ) {
		      $stringToPrint.="\t"."overlaps3P";
		      if ($printFullFiles) {
			$allStats->addCustomOverlaps3p();
		      }
		    } elsif ( ($coordInput1 < $array->[$indexArray]->{'record'}->{'range'}->getStart() &&
			       $coordInput2 > $array->[$indexArray]->{'record'}->{'range'}->getEnd() ) ) {
		      $stringToPrint.="\t"."spans";
		      if ($printFullFiles) {
			$allStats->addCustomSpans();
		      }
		    } else {
		      print STDERR "Invalid state dump: ".Dumper($array->[$indexArray])."\n";
		      terminate(9);
		    }

		  } elsif ($array->[$indexArray]->{'record'}->{'strand'} eq "-") {

		    if ( ($coordInput1 < $array->[$indexArray]->{'record'}->{'range'}->getStart() &&
			  $coordInput2 < $array->[$indexArray]->{'record'}->{'range'}->getStart() ) ) {
		      $stringToPrint.="\t"."downstream";
		      if ($printFullFiles) {
			$allStats->addCustomDownstream();
		      }
		    } elsif ( ($coordInput1 > $array->[$indexArray]->{'record'}->{'range'}->getEnd() &&
			       $coordInput2 > $array->[$indexArray]->{'record'}->{'range'}->getEnd() ) ) {
		      $stringToPrint.="\t"."upstream";
		      if ($printFullFiles) {
			$allStats->addCustomUpstream();
		      }
		    } elsif ( $array->[$indexArray]->{'record'}->{'range'}->contains($coordInput1) &&
			      $array->[$indexArray]->{'record'}->{'range'}->contains($coordInput2) ) {
		      $stringToPrint.="\t"."contained";
		      if ($printFullFiles) {
			$allStats->addCustomContained();
		      }
		    } elsif ( $coordInput1 < $array->[$indexArray]->{'record'}->{'range'}->getStart() &&
			      $array->[$indexArray]->{'record'}->{'range'}->contains($coordInput2) ) {
		      $stringToPrint.="\t"."overlaps3P";
		      if ($printFullFiles) {
			$allStats->addCustomOverlaps3p();
		      }
		    } elsif ( $array->[$indexArray]->{'record'}->{'range'}->contains($coordInput1) &&
			      $coordInput2 > $array->[$indexArray]->{'record'}->{'range'}->getEnd()  ) {
		      $stringToPrint.="\t"."overlaps5P";
		      if ($printFullFiles) {
			$allStats->addCustomOverlaps5p();
		      }
		    } elsif ( ($coordInput1 < $array->[$indexArray]->{'record'}->{'range'}->getStart() &&
			       $coordInput2 > $array->[$indexArray]->{'record'}->{'range'}->getEnd() ) ) {
		      $stringToPrint.="\t"."spans";
		      if ($printFullFiles) {
			$allStats->addCustomSpans();
		      }
		    } else {
		      print STDERR "Invalid state dump: ".Dumper($array->[$indexArray])."\n";
		      terminate(9);
		    }
		  } else {
		    print STDERR "Problem with strand in tree structure dump: ".Dumper($array->[$indexArray])."\n";
		    terminate(9);
		  }

		  if (!$printStatsOnly) {
		    print INPUTALLOUT $stringToPrint."\n";
		  }
		}
	      } else {
		#add site to stats
		if ($printFullFiles) {
		  $allStats->addToNoHits();
		}
		if (!$printStatsOnly) {
		  print INPUTNONEOUT join("\t",($idInput,$chrInput,$coordInput1,$coordInput2))."\n";
		}
	      }



	    }


	  } elsif ( $annotationMode == 3) { #coord annotation mode

	    # For SNP annotation
	    #chr15	67532214	A C myCoordinate1293
	    if ( (!$noRefBP && $line =~ /^(\S+)\s+(\d+)\s+([ACGT])\s+([ACGT])\s+(\S+)\s*$/) ||
		 ($noRefBP && $line =~ /^(\S+)\s+(\d+)\s+([ACGT])\s+(\S+)\s*$/)  ) {


	      my $chrInput     = $1;
	      my $coordInput   = $2;
	      my $refDNA       = $3;
	      my $readDNA      = $4;
	      my $idInput      = $5;

	      if (!$noRefBP) {
		$chrInput     = $1;
		$coordInput   = $2;
		$refDNA       = $3;
		$readDNA      = $4;
		$idInput      = $5;
	      } else {
		$chrInput     = $1;
		$coordInput   = $2;
		$refDNA       = "#";
		$readDNA      = $3;
		$idInput      = $4;
	      }
	      if (exists $hashOfIdInput{$idInput}) {
		print STDERR "The input id ".$idInput." was found twice\n";
		terminate(9);
	      } else {
		$hashOfIdInput{$idInput}=0;
	      }

	      if ( !ChromosomeIndexing::chrEXISTS($speciesCode,$chrInput) ) {
		print STDERR "Chromosome in file $fileName contains chromosome $chrInput which is not recognized for the current build\n";
		terminate(9);
	      }

	      if ($printFullFiles) {
		$allStats->addSite();
		$singleStats->addSite();
		$bestStats->addSite();
	      }
	      #find coordinate in tree
	      my $array=[];



	      if ($range == 0) {
		if ($arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$chrInput)]) {
		  my $indexbinst=int($coordInput/$genomebinsize);
		  my $indexbinen=int($coordInput/$genomebinsize);

		  for (my $indexbin=$indexbinst;$indexbin<=$indexbinen;$indexbin++) {
		    #print "indexbin $indexbin ".$coordInput."\n";

		    if ($arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$chrInput)]->[$indexbin]) {
		      push( @{$array} , @{$arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$chrInput)]->[$indexbin]->seekCoord($coordInput)} );
		    }
		  }

		}
	      } else {
		if ($arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$chrInput)]) {
		  my $indexbinst=int(($coordInput-$range)/$genomebinsize);
		  my $indexbinen=int(($coordInput+$range)/$genomebinsize);

		  for (my $indexbin=$indexbinst;$indexbin<=$indexbinen;$indexbin++) {
		    #print "indexbin $indexbin ".$coordInput."\n";
		    if ($arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$chrInput)]->[$indexbin]) {
		      push( @{$array},  @{$arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$chrInput)]->[$indexbin]->seekRange($coordInput-$range,$coordInput+$range)} );
		    }
		  }
		}
	      }






	      #find SNP in database
	      my $snpHash;
	      my $snpString="";

	      if ($compareWithDBSNP) {
		$snpHash   = $arrayOfSNPStructures->[ChromosomeIndexing::chr2index($speciesCode,$chrInput)]->annotateCoord($coordInput,$refDNA,$readDNA,$noRefBP);
		$snpString = formatDBSNP($snpHash);

		$allStats->addDBSNP($snpHash);
		$singleStats->addDBSNP($snpHash);
		$bestStats->addDBSNP($snpHash);

		print INPUTDBSNPOUT join("\t",($idInput,$chrInput,$coordInput))."\t".$snpString."\n";
	      }

	      #if any results were found
	      if ($#{$array} >= 0 ) {
		$array=eliminateDuplicateTranscripts($array);

		my %hashOfGenes;

		#foreach my $result ( @{$array} ) {
		for (my $indexArray=0;$indexArray<=$#{$array};$indexArray++   ) {
		  if ($printFullFiles) {
		    if (exists $hashOfGenes{$array->[$indexArray]->{'transcript'}->getGene()->getGeneCounter()} ) {
		      push( @{$hashOfGenes{$array->[$indexArray]->{'transcript'}->getGene()->getGeneCounter()}} , $array->[$indexArray]);
		    } else {
		      $hashOfGenes{$array->[$indexArray]->{'transcript'}->getGene()->getGeneCounter()}=[$array->[$indexArray]];
		    }

		    my $hashForGeneFile = {'idInput' => $idInput,
					   'symbol'  => $array->[$indexArray]->{'transcript'}->getGeneSymbol(),
					   'trans'   => $array->[$indexArray]->{'transcript'}->getId()};
		    if (!$doNotPrintGenesOut) {
		      if (exists $hashOfGenesAll{$array->[$indexArray]->{'transcript'}->getGene()->getGeneCounter()} ) {
			push( @{$hashOfGenesAll{$array->[$indexArray]->{'transcript'}->getGene()->getGeneCounter()}} , $hashForGeneFile); #$array->[$indexArray] );
		      } else {
			$hashOfGenesAll{$array->[$indexArray]->{'transcript'}->getGene()->getGeneCounter()}=[$hashForGeneFile]; #[$array->[$indexArray]];
		      }
		    }
		  }
		  $array->[$indexArray]->{'annoResult'} = $array->[$indexArray]->{'transcript'}->annotateCoordinate($coordInput);
		  $array->[$indexArray]->{'snpResult'}  = $array->[$indexArray]->{'transcript'}->annotateSNP($coordInput,$refDNA,$readDNA,$produceAASeq,$noRefBP);
		  my $stringToPrint=join("\t",($idInput,$chrInput,$coordInput))."\t".$array->[$indexArray]->{'transcript'}->getId()."\t".$array->[$indexArray]->{'transcript'}->getGeneSymbol()."\t".formatAnnotation(1,$array->[$indexArray]->{'annoResult'})."\t".formatSNPAnno($array->[$indexArray]->{'snpResult'})."\n";


		  if (!$printStatsOnly) {
		    print INPUTALLOUT $stringToPrint;
		  }
		  if ($printFullFiles) {
		    $allStats->addAnnotationSNP($array->[$indexArray]->{'annoResult'},$array->[$indexArray]->{'snpResult'});

		    if (     $array->[$indexArray]->{'annoResult'}->{'code'} == 1) {
		      if (!$printStatsOnly) {
			print INPUTEXONOUT    $stringToPrint;
		      }


		      if ($array->[$indexArray]->{'snpResult'}->{'coding'} ) {
			if (!$printStatsOnly) {
			  print INPUTCDEOUT $stringToPrint;
			}

			if ($array->[$indexArray]->{'snpResult'}->{'reliable'} ) {
			  if ($array->[$indexArray]->{'snpResult'}->{'synonymous'} == 1) {
			    if (!$printStatsOnly) {
			      print INPUTSYNOUT $stringToPrint;
			    }
			  } else {
			    if (!$printStatsOnly) {
			      print INPUTNSYNOUT  $stringToPrint;
			    }
			  }
			}
		      } else {
			if (!$printStatsOnly) {
			  print INPUTNCDEOUT $stringToPrint;
			}
		      }

		    } elsif ( $array->[$indexArray]->{'annoResult'}->{'code'} == 2) {
		      if (!$printStatsOnly) {
			print INPUTINTRONOUT  $stringToPrint;
		      }
		    } elsif ( $array->[$indexArray]->{'annoResult'}->{'code'} == 3) {
		      if (!$printStatsOnly) {
			print INPUTUPSTRMOUT  $stringToPrint;
		      }
		    } elsif ( $array->[$indexArray]->{'annoResult'}->{'code'} == 4) {
		      if (!$printStatsOnly) {
			print INPUTDWNSTRMOUT $stringToPrint;
		      }
		    } else {
		      print STDERR "Wrong code\n";
		      terminate(9);
		    }
		  }
		}		#end for each element in result array


		if ($printFullFiles) {

		  my @allGenes;
		  my %singleGenes;
		  foreach my $geneCounter (keys %hashOfGenes) {
		    my $singleBestResult=selectSingleBestTranscriptSNP( @{$hashOfGenes{$geneCounter}} );
		    push(@allGenes,@{$hashOfGenes{$geneCounter}});

		    if (exists $singleGenes{$geneCounter} ) {
		      print STDERR "Duplicate gene counters\n";
		      terminate(9);
		    } else {
		      $singleGenes{$geneCounter} = $singleBestResult->{'transcript'}->getId() ;
		    }

		    my $hashForGeneFile = {'idInput' => $idInput,
					   'symbol'  => $singleBestResult->{'transcript'}->getGeneSymbol(),
					   'trans'   => $singleBestResult->{'transcript'}->getId()};

		    $singleStats->addAnnotationSNP($singleBestResult->{'annoResult'},$singleBestResult->{'snpResult'});

		    if (!$printStatsOnly) {
		      print INPUTSINGLEOUT join("\t",($idInput,$chrInput,$coordInput))."\t".$singleBestResult->{'transcript'}->getId()."\t".$singleBestResult->{'transcript'}->getGeneSymbol()."\t".formatAnnotation(1,$singleBestResult->{'annoResult'})."\t".formatSNPAnno($singleBestResult->{'snpResult'})."\n";
		    }
		    if (!$doNotPrintGenesOut) {
		      if (exists $hashOfGenesSingle{$singleBestResult->{'transcript'}->getGene()->getGeneCounter()} ) {
			push( @{$hashOfGenesSingle{$singleBestResult->{'transcript'}->getGene()->getGeneCounter()}} , $hashForGeneFile); #$singleBestResult );
		      } else {
			$hashOfGenesSingle{$singleBestResult->{'transcript'}->getGene()->getGeneCounter()}=[$hashForGeneFile]; #[$singleBestResult];
		      }
		    }

		  }

		  my $bestResult=selectSingleBestTranscriptSNP( @allGenes );
		  $bestStats->addAnnotationSNP($bestResult->{'annoResult'},$bestResult->{'snpResult'});

		  if (!$printStatsOnly) {
		    print INPUTBESTOUT join("\t",($idInput,$chrInput,$coordInput))."\t".$bestResult->{'transcript'}->getId()."\t".$bestResult->{'transcript'}->getGeneSymbol()."\t".formatAnnotation(1,$bestResult->{'annoResult'})."\t".formatSNPAnno($bestResult->{'snpResult'})."\n";
		  }
		  my $hashForGeneFile = {'idInput' => $idInput,
					 'symbol'  => $bestResult->{'transcript'}->getGeneSymbol(),
					 'trans'   => $bestResult->{'transcript'}->getId() };

		  if (!$doNotPrintGenesOut) {
		    if (exists $hashOfGenesBest{$bestResult->{'transcript'}->getGene()->getGeneCounter()} ) {
		      push( @{$hashOfGenesBest{$bestResult->{'transcript'}->getGene()->getGeneCounter()}} , $hashForGeneFile); #$bestResult );
		    } else {
		      $hashOfGenesBest{$bestResult->{'transcript'}->getGene()->getGeneCounter()}=[$hashForGeneFile]; #[$bestResult];
		    }
		  }

		  if (!$printStatsOnly) {
		    if (!$doNotPrintXML) {
		      if ($compareWithDBSNP) {
			printXML(1,$annotationMode,$INPUTALLXML,$coordInput,$chrInput,$idInput,$bestResult,$snpString,\%singleGenes,\%hashOfGenes,\@arrayOfFieldHeader);
		      } else {
			printXML(1,$annotationMode,$INPUTALLXML,$coordInput,$chrInput,$idInput,$bestResult,-1,\%singleGenes,\%hashOfGenes,\@arrayOfFieldHeader);
		      }
		    }
		  }
		}
	      } else {
		if ($printFullFiles) {
		  $allStats->addToNoHits();
		  $bestStats->addToNoHits();
		  $singleStats->addToNoHits();
		}
		if (!$printStatsOnly) {
		  print INPUTNONEOUT join("\t",($idInput,$chrInput,$coordInput))."\n";
		}
	      }


	    } else {
	      print STDERR "Line $line in $fileName did not parse\n";
	      terminate(9);
	    }

	  } else {		#end mode 2
	    print STDERR "invalid mode\n";
	    terminate(9);
	  }

	}			# end for each line input

	close(INPUTFILE);




	if ($printFullFiles) {
	  if (!$doNotPrintGenesOut) {
	    #Print to gene files

	    foreach my $key (keys %hashOfGenesAll) {

	      if (!$singleTranscriptDatabase) {
		if (!$printStatsOnly) {
		  print GENESALLOUT $hashOfGenesAll{$key}->[0]->{'symbol'}."\t";
		}
	      }
	      my @tempArray;

	      for (my $index=0; $index<= $#{ $hashOfGenesAll{$key}} ; $index++   ) {
		push(@tempArray,$hashOfGenesAll{$key}->[$index]->{'trans'});
	      }
	      @tempArray =  keys %{{map {$_=>1} @tempArray  }};
	      if (!$printStatsOnly) {
		print GENESALLOUT join(",",@tempArray)."\t";
	      }

	      @tempArray=();
	      for (my $index=0; $index<= $#{$hashOfGenesAll{$key}} ; $index++   ) {
		push(@tempArray,$hashOfGenesAll{$key}->[$index]->{'idInput'});
	      }
	      @tempArray =  keys %{{map {$_=>1} @tempArray  }};
	      if (!$printStatsOnly) {
		print GENESALLOUT join(",",@tempArray)."\t";
		print GENESALLOUT "\n";
	      }
	    }



	    foreach my $key (keys %hashOfGenesBest) {
	      if (!$printStatsOnly) {
		print GENESBESTOUT $hashOfGenesBest{$key}->[0]->{'symbol'}."\t";
	      }
	      my @tempArray;
	      for (my $index=0; $index<= $#{ $hashOfGenesBest{$key} } ; $index++   ) {
		push(@tempArray,$hashOfGenesBest{$key}->[$index]->{'trans'});
	      }
	      @tempArray =  keys %{{map {$_=>1} @tempArray  }};
	      if (!$printStatsOnly) {
		print GENESBESTOUT join(",",@tempArray)."\t";
	      }

	      @tempArray=();
	      for (my $index=0; $index<= $#{ $hashOfGenesBest{$key} } ; $index++   ) {
		push(@tempArray,$hashOfGenesBest{$key}->[$index]->{'idInput'});
	      }
	      @tempArray =  keys %{{map {$_=>1} @tempArray  }};
	      if (!$printStatsOnly) {
		print GENESBESTOUT join(",",@tempArray)."\t";
		print GENESBESTOUT "\n";
	      }
	    }

	    if (!$singleTranscriptDatabase) {
	      foreach my $key (keys %hashOfGenesSingle) {
		if (!$printStatsOnly) {
		  print GENESSINGLEOUT $hashOfGenesSingle{$key}->[0]->{'symbol'}."\t";
		}
		my @tempArray;
		for (my $index=0; $index<= $#{ $hashOfGenesSingle{$key} } ; $index++   ) {
		  push(@tempArray,$hashOfGenesSingle{$key}->[$index]->{'trans'});
		}
		@tempArray =  keys %{{map {$_=>1} @tempArray  }};
		if (!$printStatsOnly) {
		  print GENESSINGLEOUT join(",",@tempArray)."\t";
		}

		@tempArray=();
		for (my $index=0; $index<= $#{ $hashOfGenesSingle{$key} } ; $index++   ) {
		  push(@tempArray,$hashOfGenesSingle{$key}->[$index]->{'idInput'});
		}
		@tempArray =  keys %{{map {$_=>1} @tempArray  }};
		if (!$printStatsOnly) {
		  print GENESSINGLEOUT join(",",@tempArray)."\t";
		  print GENESSINGLEOUT "\n";
		}

	      }

	    }#single transcript
	  }
	}

	if (!$customDatabase) {
	  if (!$printStatsOnly) {
	    close(INPUTNONEOUT);
	    close(INPUTALLOUT);
	  }
	  if ($printFullFiles) {
	    print INPUTALLSTAT    $allStats->printStats(" (warning: stats may not tally)");
	    if (!$doNotPrintXML) {
	      print INPUTALLSTATXML    $allStats->printStatsXML();
	    }

	    if (!$singleTranscriptDatabase) {
	      print INPUTSINGLESTAT $singleStats->printStats(" (warning: stats may not tally)");
	      print INPUTBESTSTAT   $bestStats->printStats("");
	      if (!$doNotPrintXML) {
		print INPUTSINGLESTATXML $singleStats->printStatsXML();
		print INPUTBESTSTATXML   $bestStats->printStatsXML();
	      }
	    }else{
	      print INPUTBESTSTAT   $bestStats->printStats("");
	      if (!$doNotPrintXML) {
		print INPUTBESTSTATXML   $bestStats->printStatsXML();
	      }
	    }

	    close(INPUTALLSTAT);

	    if (!$doNotPrintXML) {
	      if (!$printStatsOnly) {
		close($INPUTALLXML);
	      }
	      close(INPUTALLSTATXML);
	    }

	    if (!$singleTranscriptDatabase) {
	      if (!$printStatsOnly) {
		close(INPUTSINGLEOUT);
	      }
	      close(INPUTSINGLESTAT);


	      if (!$printStatsOnly) {
		close(INPUTBESTOUT);
	      }
	      close(INPUTBESTSTAT);

	      if (!$doNotPrintXML) {
		close(INPUTSINGLESTATXML);
		close(INPUTBESTSTATXML);
	      }
	    }else{


	      if (!$printStatsOnly) {
		close(INPUTBESTOUT);
	      }
	      close(INPUTBESTSTAT);

	      if (!$doNotPrintXML) {
		close(INPUTBESTSTATXML);
	      }

	    }

	    if (!$printStatsOnly) {
	      if (!$doNotPrintGenesOut) {
		close(GENESALLOUT);
		if (!$singleTranscriptDatabase) {
		  close(GENESSINGLEOUT);
		  close(GENESBESTOUT);
		}else{
		  close(GENESBESTOUT);
		}
	      }


	      close(INPUTEXONOUT);
	      close(INPUTINTRONOUT);
	      close(INPUTUPSTRMOUT);
	      close(INPUTDWNSTRMOUT);
	    }
	    if ($annotationMode == 3 ) {
	      if (!$printStatsOnly) {
		if ($compareWithDBSNP) {
		  close(INPUTDBSNPOUT);
		}

		close(INPUTCDEOUT);
		close(INPUTNCDEOUT);
		close(INPUTNSYNOUT);
		close(INPUTSYNOUT);
	      }
	    }
	  }
	} else {

	  print INPUTALLSTAT    $allStats->printStats("");
	  if (!$printStatsOnly) {
	    close(INPUTNONEOUT);
	    close(INPUTALLOUT);
	  }
	  if ($printFullFiles) {
	    close(INPUTALLSTAT);
	  }
	}

      }				#end for each range


      ############################################
      ##                                        ##
      ##      END ANNOTATION MODE 1,2 and 3     ##
      ##                                        ##
      ############################################








      ############################################
      ##                                        ##
      ##        BEGIN ANNOTATION MODE 4         ##
      ##                                        ##
      ############################################


    } elsif ( $annotationMode == 4) { #if not coord annotation mode 1 or 2 or 3

      my $resultOpen=open(INPUTFILE,$fileName) or die "Unable to open $fileName\n";

      if (!$printStatsOnly) {
	open(INPUTCLOSESTOUT, ">".$fileName.".input.closest.out");
	open(INPUTNONEOUT,    ">".$fileName.".input.none.out");
      }
      if ($printFullFiles) {
	open(INPUTCLOSESTSTAT, ">".$fileName.".input.closest.stat");
	if (!$doNotPrintXML) {
	  open(INPUTCLOSESTSTATXML, ">".$fileName.".input.closest.stat.xml");
	  if (!$printStatsOnly) {
	    open(INPUTCLOSESTXML, ">".$fileName.".input.closest.xml");
	  }
	}
      }

      #                        0           1     2       3            4        5        6              7           8
      my $arrayOfFieldHeader=["#INPUTID","CHR","COORD","TRANSCRIPT","GENE","POSITION","DISTANCE5p","DISTANCE3p","PARTIAL"];
      if (!$printStatsOnly) {
	print INPUTCLOSESTOUT join("\t",@{$arrayOfFieldHeader})."\n";
      }
      print INPUTNONEOUT    join("\t",("#INPUTID","CHR","COORD"))."\n";

      if (!$resultOpen) {
	print STDERR "Unable to open $fileName\n";
	terminate(9);
      }

      my %hashOfIdInput;


      my $stats;
      if ($printFullFiles) {
	$stats=StatiticsModule->new(genomeUsed      => $speciesCode,
				    databaseUsed    => $databaseName,
				    databaseFiles   => $databaseFiles,
				    range           => "###",
				    mode            => $annotationMode,
				    bpInBins        => $bpInBins,
				    numberBins      => $numberBins);
      }
      my $blankLines=0;

      while (my $line =<INPUTFILE>) {

	chomp($line);

	if (length($line) == 0) {
	  $blankLines++;
	  next;
	} else {
	  if ($blankLines>0) {
	    print STDERR "Blank lines found prior to line $line in $fileName \n";
	    terminate(9);
	  }
	}

	if ($line =~ /^(\S+)\s+(\d+)\s+(\S+)\s*$/) {
	  my $chrInput     = $1;
	  my $coordInput   = $2;
	  my $idInput      = $3;

	  if (exists $hashOfIdInput{$idInput}) {
	    print STDERR "The input id ".$idInput." was found twice\n";
	    terminate(9);
	  } else {
	    $hashOfIdInput{$idInput}=0;
	  }

	  if ( !ChromosomeIndexing::chrEXISTS($speciesCode,$chrInput) ) {
	    print STDERR "Chromosome in file $fileName contains chromosome $chrInput which is not recognized for the current build\n";
	    terminate(9);
	  }

	  my $closestTSS;
	  if ( $#{$arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$chrInput)]} > -1) {
	    $closestTSS=findClosestTSS($arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$chrInput)],$coordInput);

	    my $stringAnnoFormated=formatAnnotation(0,$closestTSS->{'transcript'}->annotateCoordinate($coordInput));
	    if (!$printStatsOnly) {
	      print INPUTCLOSESTOUT join("\t",($idInput,$chrInput,$coordInput))."\t".$closestTSS->{'transcript'}->getId()."\t".$closestTSS->{'transcript'}->getGeneSymbol()."\t". $stringAnnoFormated ."\n";
	    }

	    if ($printFullFiles) {

	      my $distance;
	      if ($closestTSS->{'transcript'}->getStrand() eq "+") {
		$distance = $coordInput-$closestTSS->{'transcript'}->get5Prime();
	      } elsif ($closestTSS->{'transcript'}->getStrand() eq "-") {
		$distance = $closestTSS->{'transcript'}->get5Prime()-$coordInput;
	      } else {
		print STDERR "Wrong strand for transcript ".Dumper($closestTSS->{'transcript'}->getStrand());
		terminate(9);
	      }

	      $stats->addClosestDistance($distance);
	      if (!$printStatsOnly) {
		if (!$doNotPrintXML) {
		  print INPUTCLOSESTXML "<INPUTID name=\"".$idInput." ".$chrInput." ".$coordInput ."\"  chr=\"".$chrInput."\"  coord=\"".$coordInput."\"   >"."\n";
		  print INPUTCLOSESTXML "<GENE name=\"".$closestTSS->{'transcript'}->getGeneSymbol()."\" link=\"NCBI\">"."\n";
		  print INPUTCLOSESTXML "<TRANS name=\"".$closestTSS->{'transcript'}->getId()."\"";
		}

		my @arraytemp=split("\t", $stringAnnoFormated );
		my $indexStart=5;
		if (!$doNotPrintXML) {
		  for (my $indexFIELD=$indexStart;$indexFIELD<=$#{$arrayOfFieldHeader};$indexFIELD++) {
		    print INPUTCLOSESTXML " ".lc($arrayOfFieldHeader->[$indexFIELD])."=\"". $arraytemp[$indexFIELD-$indexStart] ."\"";
		  }
		  print INPUTCLOSESTXML " color=\"".$colorNormal."\" ";
		  print INPUTCLOSESTXML " />"."\n";
		  print INPUTCLOSESTXML "</GENE>"."\n";
		  print INPUTCLOSESTXML "</INPUTID>"."\n";
		}
	      }
	    }
	  } else {
	    if (!$printStatsOnly) {
	      print INPUTNONEOUT join("\t",($idInput,$chrInput,$coordInput))."\n";
	    }
	  }


	} else {
	  print STDERR "Line $line in $fileName did not parse\n";
	  terminate(9);
	}

      }				# end for each line input

      if (!$printStatsOnly) {
	close(INPUTCLOSESTOUT);
	close(INPUTNONEOUT);
      }
      close(INPUTFILE);
      if ($printFullFiles) {
	print INPUTCLOSESTSTAT    $stats->printStats("");
	if (!$doNotPrintXML) {
	  print INPUTCLOSESTSTATXML $stats->printStatsXML();
	}

	close(INPUTCLOSESTSTAT);
	if (!$doNotPrintXML) {
	  close(INPUTCLOSESTSTATXML);
	  if (!$printStatsOnly) {
	    close(INPUTCLOSESTXML);
	  }
	}
      }

      ############################################
      ##                                        ##
      ##         END ANNOTATION MODE 4          ##
      ##                                        ##
      ############################################















      ############################################
      ##                                        ##
      ##        BEGIN ANNOTATION MODE 5         ##
      ##                                        ##
      ############################################

    } elsif ( $annotationMode == 5) { #INSERTION, DELETION and TRANSLOCATION

      foreach my $range (@arrayOfRanges) { #For each range to consider

	open(INPUTALLOUT,     ">".$fileName.".".$range.".".$databaseName.".input.all.out");
	if ($printFullFiles) {
	  if (!$doNotPrintXML) {
	    open($INPUTALLXML,    ">".$fileName.".".$range.".".$databaseName.".input.all.xml");
	  }

	  if (!$doNotPrintGenesOut) {
	    open(GENESALLOUT,     ">".$fileName.".".$range.".".$databaseName.".genes.all.out");
	  }
	}

	open(INPUTNONEOUT,    ">".$fileName.".".$range.".".$databaseName.".input.none.out");


	my $resultOpen=open(INPUTFILE,$fileName) or die "Unable to open $fileName\n";

	if (!$resultOpen) {
	  print STDERR "Unable to open $fileName\n";
	  terminate(9);
	}

	my %hashOfIdInput;
	my %hashOfGenesAll;
	my $linesXMLinsert    = "";
	my $linesXMLdeletion  = "";
	my $linesXMLtrans     = "";

	my $blankLines=0;


	while (my $line =<INPUTFILE>) {
	  #chr15	67532214	A C myCoordinate1293
	  my $hashInput;
	  #INS chr1 23493 GATGATGAT  insertion1239
	  #                      1       2        3          4

	  chomp($line);

	  if (length($line) == 0) {
	    $blankLines++;
	    next;
	  } else {
	    if ($blankLines>0) {
	      print STDERR "Blank lines found prior to line $line in $fileName \n";
	      terminate(9);
	    }
	  }

	  if ($line =~ /^INS\s+(\S+)\s+(\d+)\s+([ACGT]+)\s+(\S+)\s*$/) {


	    $hashInput={'type'     => 1,
			'chr'      => $1,
			'coord'    => $2,
			'seq'      => $3,
			'idInput'  => $4};

	    if (exists $hashOfIdInput{ $hashInput->{'idInput'} }) {
	      print STDERR "The input id ".$hashInput->{'idInput'}." was found twice\n";
	      terminate(9);
	    } else {
	      $hashOfIdInput{ $hashInput->{'idInput'} }=0;
	    }

	    if ( !ChromosomeIndexing::chrEXISTS($speciesCode,$hashInput->{'chr'}) ) {
	      print STDERR "Chromosome in file $fileName contains chromosome ".$hashInput->{'chr'}." which is not recognized for the current build\n";
	      terminate(9);
	    }

	    #DEL chr2 32432 32439  deletion543
	    #                         1       2       3       4
	  } elsif ($line =~ /^DEL\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s*$/) {

	    $hashInput={'type'     => 2,
			'chr'      => $1,
			'coord1'   => $2,
			'coord2'   => $3,
			'idInput'  => $4};
	    if (exists $hashOfIdInput{ $hashInput->{'idInput'} }) {
	      print STDERR "The input id ".$hashInput->{'idInput'}." was found twice\n";
	      terminate(9);
	    } else {
	      $hashOfIdInput{ $hashInput->{'idInput'} }=0;
	    }

	    if ( !ChromosomeIndexing::chrEXISTS($speciesCode,$hashInput->{'chr'}) ) {
	      print STDERR "Chromosome in file $fileName contains chromosome ".$hashInput->{'chr'}." which is not recognized for the current build\n";
	      terminate(9);
	    }

	    if ($hashInput->{'coord1'} > $hashInput->{'coord2'}) {
	      print STDERR "Line $line in $fileName did not parse, coordinate ".$hashInput->{'coord1'}." cannot be lesser than ".$hashInput->{'coord2'}."\n";
	      terminate(9);
	    }
	    #TRANS chr1 239499 + chr8 92344 - trans1031
	    #                           1       2        3         4       5        6         7
	  } elsif ($line =~ /^TRANS\s+(\S+)\s+(\d+)\s+([\+-]+)\s+(\S+)\s+(\d+)\s+([\+-]+)\s+(\S+)\s*$/) {
	    $hashInput={'type'     => 3,
			'chr1'     => $1,
			'coord1'   => $2,
			'strand1'  => $3,
			'chr2'     => $4,
			'coord2'   => $5,
			'strand2'  => $6,
			'idInput'  => $7};

	    if (exists $hashOfIdInput{ $hashInput->{'idInput'} }) {
	      print STDERR "The input id ".$hashInput->{'idInput'}." was found twice\n";
	      terminate(9);
	    } else {
	      $hashOfIdInput{ $hashInput->{'idInput'} }=0;
	    }

	    if ( !ChromosomeIndexing::chrEXISTS($speciesCode,$hashInput->{'chr1'}) ) {
	      print STDERR "Chromosome in file $fileName contains chromosome ".$hashInput->{'chr1'}." which is not recognized for the current build\n";
	      terminate(9);
	    }

	    if ( !ChromosomeIndexing::chrEXISTS($speciesCode,$hashInput->{'chr2'}) ) {
	      print STDERR "Chromosome in file $fileName contains chromosome ".$hashInput->{'chr2'}." which is not recognized for the current build\n";
	      terminate(9);
	    }

	  } else {
	    print STDERR "Line $line in $fileName did not parse\n";
	    terminate(9);
	  }












	  if (    $hashInput->{'type'} == 1) { #INSERTION

	    my $array=[];

	    #	    if ($range == 0) {
	    #	      if ( $arrayOfTrees->[ ChromosomeIndexing::chr2index($speciesCode,$hashInput->{'chr'}) ]) {
	    #		$array  =  $arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$hashInput->{'chr'})]->seekCoord($hashInput->{'coord'});
	    #	      }
	    #	    } else {
	    #	      if ($arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$hashInput->{'chr'})]) {
	    #		$array  =  $arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$hashInput->{'chr'})]->seekRange($hashInput->{'coord'}-$range,
	    #														       $hashInput->{'coord'}+$range);
	    #	      }
	    #	    }


	    if ($range == 0) {
	      if ($arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$hashInput->{'chr'})]) {
		my $indexbinst=int($hashInput->{'coord'}/$genomebinsize);
		my $indexbinen=int($hashInput->{'coord'}/$genomebinsize);

		for (my $indexbin=$indexbinst;$indexbin<=$indexbinen;$indexbin++) {
		  if ($arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$hashInput->{'chr'})]->[$indexbin]) {
		    push( @{$array} , @{$arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$hashInput->{'chr'})]->[$indexbin]->seekCoord($hashInput->{'coord'})} );
		  }
		}

	      }
	    } else {
	      if ($arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$hashInput->{'chr'})]) {
		my $indexbinst=int(($hashInput->{'coord'}-$range)/$genomebinsize);
		my $indexbinen=int(($hashInput->{'coord'}+$range)/$genomebinsize);

		for (my $indexbin=$indexbinst;$indexbin<=$indexbinen;$indexbin++) {
		  if ($arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$hashInput->{'chr'})]->[$indexbin]) {
		    push( @{$array},  @{$arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$hashInput->{'chr'})]->[$indexbin]->seekRange($hashInput->{'coord'}-$range,$hashInput->{'coord'}+$range)} );
		  }
		}
	      }
	    }



	    if ($#{$array} >= 0 ) {
	      $array=eliminateDuplicateTranscripts($array);
	      if ($printFullFiles) {
		$linesXMLinsert.="<INPUTID name=\"".$hashInput->{'idInput'}." ".$hashInput->{'chr'}." ".$hashInput->{'coord'}."\" chr=\"".$hashInput->{'chr'}."\" coord=\"".$hashInput->{'coord'}."\">"."\n";

		if (!$singleTranscriptDatabase) {
		  $linesXMLinsert.="<GENE name=\"".$array->[0]->{'transcript'}->getGeneSymbol()."\" link=\"NCBI\">"."\n";
		}
	      }
	      for (my $indexArray=0;$indexArray<=$#{$array};$indexArray++   ) {
		#print $INPUTALLXML "<TRANS name=\"".$array->[$indexArray]->{'transcript'}->getId()."\" ";
		if ($printFullFiles) {
		  $linesXMLinsert.="<TRANS name=\"".$array->[$indexArray]->{'transcript'}->getId()."\" ";

		  my $hashForGeneFile = {'idInput' => $hashInput->{idInput},
					 'symbol'  => $array->[$indexArray]->{'transcript'}->getGeneSymbol(),
					 'trans'   => $array->[$indexArray]->{'transcript'}->getId()};

		  if (!$doNotPrintGenesOut) {
		    if (exists $hashOfGenesAll{$array->[$indexArray]->{'transcript'}->getGene()->getGeneCounter()} ) {
		      push( @{$hashOfGenesAll{$array->[$indexArray]->{'transcript'}->getGene()->getGeneCounter()}} , $hashForGeneFile); #$array->[$indexArray] );
		    } else {
		      $hashOfGenesAll{$array->[$indexArray]->{'transcript'}->getGene()->getGeneCounter()}=[$hashForGeneFile]; #[$array->[$indexArray]];
		    }
		  }
		}

		$array->[$indexArray]->{'annoResult'}    = $array->[$indexArray]->{'transcript'}->annotateCoordinate($hashInput->{'coord'});
		$array->[$indexArray]->{'insertResult'}  = $array->[$indexArray]->{'transcript'}->annotateINSERT($hashInput->{'coord'},$hashInput->{'seq'});


		my $stringToPrint=join("\t",($hashInput->{'idInput'},$hashInput->{'chr'},$hashInput->{'coord'}))."\t".$array->[$indexArray]->{'transcript'}->getId()."\t".$array->[$indexArray]->{'transcript'}->getGeneSymbol()."\t".formatAnnotation(1,$array->[$indexArray]->{'annoResult'})."\t".formatInsertAnno($array->[$indexArray]->{'insertResult'},$array->[$indexArray]->{'annoResult'})."\n";

		print INPUTALLOUT $stringToPrint;

		if ($printFullFiles) {
		  my @arraytemp=split("\t",formatAnnotation(1,$array->[$indexArray]->{'annoResult'}));
		  my @arrayOfFieldHeader = ("POSITION","INDEX","DISTANCE5p","DISTANCE3p","PARTIAL");

		  for (my $indexFIELD=0;$indexFIELD<5;$indexFIELD++) {
		    #print $INPUTALLXML "".lc($arrayOfFieldHeader[$indexFIELD])."=\"". $arraytemp[$indexFIELD] ."\" ";
		    $linesXMLinsert.="".lc($arrayOfFieldHeader[$indexFIELD])."=\"". $arraytemp[$indexFIELD] ."\" ";
		  }

		  @arraytemp=split("\t", formatInsertAnno($array->[$indexArray]->{'insertResult'},$array->[$indexArray]->{'annoResult'}) );
		  if ($#arraytemp == 2) {
		    $linesXMLinsert.=" indexaa=\"".$arraytemp[0]."\"  seqref=\"".$arraytemp[1]."\"  seqread=\"".$arraytemp[2]."\"  comment=\"\" ";
		  } else {
		    $linesXMLinsert.=" indexaa=\"\"  seqread=\"\"  seqref=\"\"  comment=\"".$arraytemp[0]."\" ";
		  }

		  $linesXMLinsert.= " color=\"".$colorNormal."\" ";
		  $linesXMLinsert.="link=\"NCBI\"/>"."\n";
		}
	      }			#end for each element in result array

	      if ($printFullFiles) {
		if (!$singleTranscriptDatabase) {
		  #print $INPUTALLXML "</GENE>"."\n";
		  $linesXMLinsert.="</GENE>"."\n";
		}

		#print $INPUTALLXML "</INPUTID>"."\n";
		$linesXMLinsert.="</INPUTID>"."\n";
	      }
	    } else {
	      print INPUTNONEOUT join("\t",($hashInput->{'idInput'},$hashInput->{'chr'},$hashInput->{'coord'}))."\n";
	    }












	  } elsif ($hashInput->{'type'} == 2) { #DELETION
	    my $array=[];
	    #	    if ($arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$hashInput->{'chr'})]) {
	    #	      $array  =  $arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$hashInput->{'chr'})]->seekRange($hashInput->{'coord1'}-$range,
	    #														     $hashInput->{'coord2'}+$range);
	    #	    }

	    if ($arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$hashInput->{'chr'})]) {
	      my $indexbinst=int(($hashInput->{'coord1'}-$range)/$genomebinsize);
	      my $indexbinen=int(($hashInput->{'coord2'}+$range)/$genomebinsize);

	      for (my $indexbin=$indexbinst;$indexbin<=$indexbinen;$indexbin++) {
		if ($arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$hashInput->{'chr'})]->[$indexbin]) {
		  push( @{$array},  @{$arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$hashInput->{'chr'})]->[$indexbin]->seekRange($hashInput->{'coord1'}-$range,$hashInput->{'coord2'}+$range)} );
		}
	      }
	    }

	    if ($#{$array} >= 0 ) {
	      $array=eliminateDuplicateTranscripts($array);
	      if ($printFullFiles) {
		$linesXMLdeletion.="<INPUTID name=\"".$hashInput->{'idInput'}." ".$hashInput->{'chr'}." ".$hashInput->{'coord1'}."-".$hashInput->{'coord2'}."\" chr=\"".$hashInput->{'chr'}."\" coord=\"".$hashInput->{'coord1'}."-".$hashInput->{'coord2'}."\">"."\n";

		if (!$singleTranscriptDatabase) {
		  $linesXMLdeletion.="<GENE name=\"".$array->[0]->{'transcript'}->getGeneSymbol()."\" link=\"NCBI\">"."\n";
		}
	      }

	      for (my $indexArray=0;$indexArray<=$#{$array};$indexArray++   ) {
		if ($printFullFiles) {
		  $linesXMLdeletion.="<TRANS name=\"".$array->[$indexArray]->{'transcript'}->getId()."\" ";

		  my $hashForGeneFile = {'idInput' => $hashInput->{idInput},
					 'symbol'  => $array->[$indexArray]->{'transcript'}->getGeneSymbol(),
					 'trans'   => $array->[$indexArray]->{'transcript'}->getId()};

		  if (!$doNotPrintGenesOut) {
		    if (exists $hashOfGenesAll{$array->[$indexArray]->{'transcript'}->getGene()->getGeneCounter()} ) {
		      push( @{$hashOfGenesAll{$array->[$indexArray]->{'transcript'}->getGene()->getGeneCounter()}} , $hashForGeneFile); #$array->[$indexArray] );
		    } else {
		      $hashOfGenesAll{$array->[$indexArray]->{'transcript'}->getGene()->getGeneCounter()}=[$hashForGeneFile]; #[$array->[$indexArray]];
		    }
		  }

		}
		$array->[$indexArray]->{'deletionResult'}  = $array->[$indexArray]->{'transcript'}->annotateDELETION($hashInput->{'coord1'},
														     $hashInput->{'coord2'});
		$array->[$indexArray]->{'annoResult'}      = $array->[$indexArray]->{'transcript'}->annotateRange($hashInput->{'coord1'},
														  $hashInput->{'coord2'});

		my $stringToPrint=join("\t",($hashInput->{'idInput'},$hashInput->{'chr'},$hashInput->{'coord1'},$hashInput->{'coord2'}))."\t".$array->[$indexArray]->{'transcript'}->getId()."\t".$array->[$indexArray]->{'transcript'}->getGeneSymbol()."\t".formatIntervalAnnotation( $array->[$indexArray]->{'annoResult'} )."\t".formatDeletionAnno($array->[$indexArray]->{'deletionResult'},$array->[$indexArray]->{'annoResult'})."\n";
		print INPUTALLOUT $stringToPrint;
		if ($printFullFiles) {
		  my @arraytemp=split("\t",formatIntervalAnnotation( $array->[$indexArray]->{'annoResult'} ));
		  my @arrayOfFieldHeader = ("POSITION","INDEXexons","INDEXintrons","PARTIAL");

		  for (my $indexFIELD=0;$indexFIELD<4;$indexFIELD++) {
		    $linesXMLdeletion.="".lc($arrayOfFieldHeader[$indexFIELD])."=\"". $arraytemp[$indexFIELD] ."\" ";
		  }


		  @arraytemp=split("\t",formatDeletionAnno($array->[$indexArray]->{'deletionResult'},$array->[$indexArray]->{'annoResult'}));
		  if ($#arraytemp == 1) {
		    $linesXMLdeletion.=" seqref=\"".$arraytemp[0]."\" seqread=\"".$arraytemp[1]."\" comment=\"\"";
		  } else {
		    $linesXMLdeletion.=" seqref=\"\" seqread=\"\" comment=\"".$arraytemp[0]."\"";
		  }
		  $linesXMLdeletion.= " color=\"".$colorNormal."\" ";
		  $linesXMLdeletion.="link=\"NCBI\"/>"."\n";

		}

	      }			#end for each element in result array

	      if ($printFullFiles) {
		if (!$singleTranscriptDatabase) {
		  $linesXMLdeletion.="</GENE>"."\n";
		}

		$linesXMLdeletion.="</INPUTID>"."\n";
	      }
	    } else {
	      print INPUTNONEOUT join("\t",($hashInput->{'idInput'},$hashInput->{'chr'},$hashInput->{'coord1'},$hashInput->{'coord2'}))."\n";
	    }






	  } elsif ($hashInput->{'type'} == 3) { #TRANSLOCATION

	    my @resultsSameDir1;
	    my @resultsDiffDir1;
	    my @resultsSameDir2;
	    my @resultsDiffDir2;

	    if ( $arrayOfTrees->[ ChromosomeIndexing::chr2index($speciesCode,$hashInput->{'chr1'}) ]) {
	      my $array1=[];
	      #=  $arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$hashInput->{'chr1'})]->seekCoord($hashInput->{'coord1'});
	      my $indexbinst=int($hashInput->{'coord1'}/$genomebinsize);
	      my $indexbinen=int($hashInput->{'coord1'}/$genomebinsize);
	      for (my $indexbin=$indexbinst;$indexbin<=$indexbinen;$indexbin++) {
		if ($arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$hashInput->{'chr1'})]->[$indexbin]) {
		  push( @{$array1} , @{$arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$hashInput->{'chr1'})]->[$indexbin]->seekCoord($hashInput->{'coord1'})} );
		}
	      }

	      if (  ($#{$array1} >= 0 )  ) {
		$array1=eliminateDuplicateTranscripts($array1);
		for (my $indexArray=0;$indexArray<=$#{$array1};$indexArray++   ) {

		  my $hashForGeneFile = {'idInput' => $hashInput->{idInput},
					 'symbol'  => $array1->[$indexArray]->{'transcript'}->getGeneSymbol(),
					 'trans'   => $array1->[$indexArray]->{'transcript'}->getId()};

		  if (!$doNotPrintGenesOut) {
		    if (exists $hashOfGenesAll{$array1->[$indexArray]->{'transcript'}->getGene()->getGeneCounter()} ) {
		      push( @{$hashOfGenesAll{$array1->[$indexArray]->{'transcript'}->getGene()->getGeneCounter()}} , $hashForGeneFile); #$array->[$indexArray] );
		    } else {
		      $hashOfGenesAll{$array1->[$indexArray]->{'transcript'}->getGene()->getGeneCounter()}=[$hashForGeneFile]; #[$array->[$indexArray]];
		    }
		  }

		  if (      ($hashInput->{'strand1'} eq "+") && ($array1->[$indexArray]->{'transcript'}->getStrand() eq "+") ) {
		    $array1->[$indexArray]->{'transResult1'}  = $array1->[$indexArray]->{'transcript'}->annotateDELETION($hashInput->{'coord1'},
															 $array1->[$indexArray]->{'transcript'}->get3Prime());
		    my $hashTemp={'transID'   => $array1->[$indexArray]->{'transcript'}->getId(),
				  'geneID'    => $array1->[$indexArray]->{'transcript'}->getGeneSymbol(),
				  'transAnno' => formatTransAnno($array1->[$indexArray]->{'transResult1'}) };
		    push(@resultsSameDir1,$hashTemp);
		  } elsif (  ($hashInput->{'strand1'} eq "+") && ($array1->[$indexArray]->{'transcript'}->getStrand() eq "-") ) {
		    $array1->[$indexArray]->{'transResult1'}  = $array1->[$indexArray]->{'transcript'}->annotateDELETION($hashInput->{'coord1'},
															 $array1->[$indexArray]->{'transcript'}->get5Prime());
		    my $hashTemp={'transID'   => $array1->[$indexArray]->{'transcript'}->getId(),
				  'geneID'    => $array1->[$indexArray]->{'transcript'}->getGeneSymbol(),
				  'transAnno' => formatTransAnno($array1->[$indexArray]->{'transResult1'}) };

		    push(@resultsDiffDir1,$hashTemp);
		  } elsif (  ($hashInput->{'strand1'} eq "-") && ($array1->[$indexArray]->{'transcript'}->getStrand() eq "+") ) {
		    $array1->[$indexArray]->{'transResult1'}  = $array1->[$indexArray]->{'transcript'}->annotateDELETION($array1->[$indexArray]->{'transcript'}->get5Prime(),
															 $hashInput->{'coord1'} );
		    my $hashTemp={'transID'   => $array1->[$indexArray]->{'transcript'}->getId(),
				  'geneID'    => $array1->[$indexArray]->{'transcript'}->getGeneSymbol(),
				  'transAnno' => formatTransAnno($array1->[$indexArray]->{'transResult1'}) };

		    push(@resultsDiffDir1,$hashTemp);

		  } elsif (  ($hashInput->{'strand1'} eq "-") && ($array1->[$indexArray]->{'transcript'}->getStrand() eq "-") ) {

		    $array1->[$indexArray]->{'transResult1'}  = $array1->[$indexArray]->{'transcript'}->annotateDELETION($array1->[$indexArray]->{'transcript'}->get3Prime(),
															 $hashInput->{'coord1'} );
		    my $hashTemp={'transID'   => $array1->[$indexArray]->{'transcript'}->getId(),
				  'geneID'    => $array1->[$indexArray]->{'transcript'}->getGeneSymbol(),
				  'transAnno' => formatTransAnno($array1->[$indexArray]->{'transResult1'}) };

		    push(@resultsSameDir1,$hashTemp);

		  } else {
		    print STDERR "Invalid state for translocation 1\n";
		    terminate(9);
		  }
		}

	      }
	    }

	    if ( $arrayOfTrees->[ ChromosomeIndexing::chr2index($speciesCode,$hashInput->{'chr2'}) ]) {
	      my $array2=[];
	      #		=  $arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$hashInput->{'chr2'})]->seekCoord($hashInput->{'coord2'});
	      my $indexbinst=int($hashInput->{'coord2'}/$genomebinsize);
	      my $indexbinen=int($hashInput->{'coord2'}/$genomebinsize);

	      for (my $indexbin=$indexbinst;$indexbin<=$indexbinen;$indexbin++) {
		if ($arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$hashInput->{'chr2'})]->[$indexbin]) {
		  push( @{$array2} , @{$arrayOfTrees->[ChromosomeIndexing::chr2index($speciesCode,$hashInput->{'chr2'})]->[$indexbin]->seekCoord($hashInput->{'coord2'})} );
		}
	      }

	      if (  ($#{$array2} >= 0 )  ) {
		$array2=eliminateDuplicateTranscripts($array2);
		for (my $indexArray=0;$indexArray<=$#{$array2};$indexArray++   ) {


		  my $hashForGeneFile = {'idInput' => $hashInput->{idInput},
					 'symbol'  => $array2->[$indexArray]->{'transcript'}->getGeneSymbol(),
					 'trans'   => $array2->[$indexArray]->{'transcript'}->getId()};

		  if (!$doNotPrintGenesOut) {
		    if (exists $hashOfGenesAll{$array2->[$indexArray]->{'transcript'}->getGene()->getGeneCounter()} ) {
		      push( @{$hashOfGenesAll{$array2->[$indexArray]->{'transcript'}->getGene()->getGeneCounter()}} , $hashForGeneFile); #$array->[$indexArray] );
		    } else {
		      $hashOfGenesAll{$array2->[$indexArray]->{'transcript'}->getGene()->getGeneCounter()}=[$hashForGeneFile]; #[$array->[$indexArray]];
		    }
		  }

		  if (      ($hashInput->{'strand2'} eq "+") && ($array2->[$indexArray]->{'transcript'}->getStrand() eq "+") ) {
		    $array2->[$indexArray]->{'transResult2'}  = $array2->[$indexArray]->{'transcript'}->annotateDELETION($array2->[$indexArray]->{'transcript'}->get5Prime() , $hashInput->{'coord2'});

		    my $hashTemp={'transID'   => $array2->[$indexArray]->{'transcript'}->getId(),
				  'geneID'    => $array2->[$indexArray]->{'transcript'}->getGeneSymbol(),
				  'transAnno' => formatTransAnno($array2->[$indexArray]->{'transResult2'}) };

		    push(@resultsSameDir2,$hashTemp);
		  } elsif (  ($hashInput->{'strand2'} eq "+") && ($array2->[$indexArray]->{'transcript'}->getStrand() eq "-") ) {
		    $array2->[$indexArray]->{'transResult2'}  = $array2->[$indexArray]->{'transcript'}->annotateDELETION($array2->[$indexArray]->{'transcript'}->get3Prime(),$hashInput->{'coord2'} );
		    my $hashTemp={'transID'   => $array2->[$indexArray]->{'transcript'}->getId(),
				  'geneID'    => $array2->[$indexArray]->{'transcript'}->getGeneSymbol(),
				  'transAnno' => formatTransAnno($array2->[$indexArray]->{'transResult2'}) };

		    push(@resultsDiffDir2,$hashTemp);

		  } elsif (  ($hashInput->{'strand2'} eq "-") && ($array2->[$indexArray]->{'transcript'}->getStrand() eq "+") ) {
		    $array2->[$indexArray]->{'transResult2'}  = $array2->[$indexArray]->{'transcript'}->annotateDELETION($hashInput->{'coord2'},
															 $array2->[$indexArray]->{'transcript'}->get3Prime() );
		    my $hashTemp={'transID'   => $array2->[$indexArray]->{'transcript'}->getId(),
				  'geneID'    => $array2->[$indexArray]->{'transcript'}->getGeneSymbol(),
				  'transAnno' => formatTransAnno($array2->[$indexArray]->{'transResult2'}) };
		    push(@resultsDiffDir2,$hashTemp);

		  } elsif (  ($hashInput->{'strand2'} eq "-") && ($array2->[$indexArray]->{'transcript'}->getStrand() eq "-") ) {
		    $array2->[$indexArray]->{'transResult2'}  = $array2->[$indexArray]->{'transcript'}->annotateDELETION($hashInput->{'coord2'},
															 $array2->[$indexArray]->{'transcript'}->get5Prime());
		    my $hashTemp={'transID'   => $array2->[$indexArray]->{'transcript'}->getId(),
				  'geneID'    => $array2->[$indexArray]->{'transcript'}->getGeneSymbol(),
				  'transAnno' => formatTransAnno($array2->[$indexArray]->{'transResult2'}) };
		    push(@resultsSameDir2,$hashTemp);
		  } else {
		    print STDERR "Invalid state for translocation 2\n";
		    terminate(9);
		  }

		}
	      }
	    }

	    my $stringToPrint="";
	    my $stringToPrintXML="";

	    for (my $index=0;$index<=max($#resultsSameDir1,$#resultsSameDir2);$index++) {
	      $stringToPrint.=join("\t",($hashInput->{'idInput'},$hashInput->{'chr1'},$hashInput->{'coord1'},$hashInput->{'chr2'},$hashInput->{'coord2'}))."\t";
	      if ($index <=$#resultsSameDir1) {
		$stringToPrint.="SD1=".$resultsSameDir1[$index]->{'transID'}."\t".$resultsSameDir1[$index]->{'geneID'}."\t".$resultsSameDir1[$index]->{'transAnno'}."\t";

		if ($printFullFiles) {
		  $stringToPrintXML.="<GENE name=\"".$resultsSameDir1[$index]->{'geneID'} ."\" link=\"NCBI\">"."\n";
		  $stringToPrintXML.= "<TRANS name=\"".$resultsSameDir1[$index]->{'transID'}."\"";
		  $stringToPrintXML.= " position=\"same strand as coord 1\" comment=\"".$resultsSameDir1[$index]->{'transAnno'}."\" ";
		  $stringToPrintXML.= " color=\"".$colorNormal."\" ";
		  $stringToPrintXML.= " />"."\n";
		  $stringToPrintXML.="</GENE>"."\n";
		}
	      }
	      if ($index <=$#resultsSameDir2) {
		$stringToPrint.="SD2=".$resultsSameDir2[$index]->{'transID'}."\t".$resultsSameDir2[$index]->{'geneID'}."\t".$resultsSameDir2[$index]->{'transAnno'}."\t";
		if ($printFullFiles) {
		  $stringToPrintXML.="<GENE name=\"".$resultsSameDir2[$index]->{'geneID'} ."\" link=\"NCBI\">"."\n";
		  $stringToPrintXML.= "<TRANS name=\"".$resultsSameDir2[$index]->{'transID'}."\"";
		  $stringToPrintXML.= " position=\"same strand as coord 2\" comment=\"".$resultsSameDir2[$index]->{'transAnno'}."\" ";
		  $stringToPrintXML.= " color=\"".$colorNormal."\" ";
		  $stringToPrintXML.= " />"."\n";
		  $stringToPrintXML.="</GENE>"."\n";
		}

	      }
	      $stringToPrint.="\n";
	    }

	    for (my $index=0;$index<=max($#resultsDiffDir1,$#resultsDiffDir2);$index++) {
	      $stringToPrint.=join("\t",($hashInput->{'idInput'},$hashInput->{'chr1'},$hashInput->{'coord1'},$hashInput->{'chr2'},$hashInput->{'coord2'}))."\t";
	      if ($index <=$#resultsDiffDir1) {
		$stringToPrint.="DD1=".$resultsDiffDir1[$index]->{'transID'}."\t".$resultsDiffDir1[$index]->{'geneID'}."\t".$resultsDiffDir1[$index]->{'transAnno'}."\t";
		if ($printFullFiles) {
		  $stringToPrintXML.="<GENE name=\"".$resultsDiffDir1[$index]->{'geneID'} ."\" link=\"NCBI\">"."\n";
		  $stringToPrintXML.= "<TRANS name=\"".$resultsDiffDir1[$index]->{'transID'}."\"";
		  $stringToPrintXML.= " position=\"different strand as coord 1\" comment=\"".$resultsDiffDir1[$index]->{'transAnno'}."\" ";
		  $stringToPrintXML.= " color=\"".$colorNormal."\" ";
		  $stringToPrintXML.= " />"."\n";
		  $stringToPrintXML.="</GENE>"."\n";
		}
	      }
	      if ($index <=$#resultsDiffDir2) {
		$stringToPrint.="DD2=".$resultsDiffDir2[$index]->{'transID'}."\t".$resultsDiffDir2[$index]->{'geneID'}."\t".$resultsDiffDir2[$index]->{'transAnno'}."\t";
		if ($printFullFiles) {
		  $stringToPrintXML.="<GENE name=\"".$resultsDiffDir2[$index]->{'geneID'} ."\" link=\"NCBI\">"."\n";
		  $stringToPrintXML.= "<TRANS name=\"".$resultsDiffDir2[$index]->{'transID'}."\"";
		  $stringToPrintXML.= " position=\"different strand as coord 2\"  comment=\"".$resultsDiffDir2[$index]->{'transAnno'}."\" ";
		  $stringToPrintXML.= " color=\"".$colorNormal."\" ";
		  $stringToPrintXML.= " />"."\n";
		  $stringToPrintXML.="</GENE>"."\n";
		}
	      }
	      $stringToPrint.="\n";
	    }

	    if ( (max($#resultsSameDir1,$#resultsSameDir2) != -1) || (max($#resultsDiffDir1,$#resultsDiffDir2) != -1)  ) {
	      print INPUTALLOUT $stringToPrint;
	      if ($printFullFiles) {
		$linesXMLtrans.="<INPUTID name=\"".$hashInput->{'idInput'}." ".$hashInput->{'chr1'}." ".$hashInput->{'coord1'}." ". $hashInput->{'chr2'} .$hashInput->{'coord2'}."\" chr=\"".$hashInput->{'chr1'}."-".$hashInput->{'chr2'}."\" coord=\"".$hashInput->{'coord1'}."-".$hashInput->{'coord2'}."\">"."\n";
		$linesXMLtrans.=$stringToPrintXML."\n";
		$linesXMLtrans.="</INPUTID>"."\n";
	      }
	    } else {
	      print INPUTNONEOUT join("\t",($hashInput->{'idInput'},$hashInput->{'chr1'},$hashInput->{'coord1'},$hashInput->{'chr2'},$hashInput->{'coord2'}))."\n";
	    }

	  }

	}			# end for each line input

	if ($printFullFiles) {
	  if (!$doNotPrintGenesOut) {
	    foreach my $key (keys %hashOfGenesAll) {
	      print GENESALLOUT $hashOfGenesAll{$key}->[0]->{'symbol'}."\t";
	      my @tempArray;
	      for (my $index=0; $index<= $#{ $hashOfGenesAll{$key} } ; $index++   ) {
		push(@tempArray,$hashOfGenesAll{$key}->[$index]->{'trans'});
	      }
	      @tempArray =  keys %{{map {$_=>1} @tempArray  }};
	      print GENESALLOUT join(",",@tempArray)."\t";

	      @tempArray=();
	      for (my $index=0; $index<= $#{ $hashOfGenesAll{$key} } ; $index++   ) {
		push(@tempArray,$hashOfGenesAll{$key}->[$index]->{'idInput'});
	      }
	      @tempArray =  keys %{{map {$_=>1} @tempArray  }};
	      print GENESALLOUT join(",",@tempArray)."\t";

	      print GENESALLOUT "\n";
	    }
	  }

	  if (!$doNotPrintXML) {
	    print $INPUTALLXML "<INSERTIONS>"."\n";
	    if ($linesXMLinsert eq "") {
	      print $INPUTALLXML "<fake></fake>\n";
	    } else {
	      print $INPUTALLXML $linesXMLinsert."\n";
	    }
	    print $INPUTALLXML "</INSERTIONS>"."\n";

	    print $INPUTALLXML "<DELETIONS>"."\n";
	    if ($linesXMLdeletion eq "") {
	      print $INPUTALLXML "<fake></fake>\n";
	    } else {
	      print $INPUTALLXML $linesXMLdeletion."\n";
	    }
	    print $INPUTALLXML "</DELETIONS>"."\n";

	    print $INPUTALLXML "<TRANSLOCATIONS>"."\n";
	    if ($linesXMLtrans eq "") {
	      print $INPUTALLXML "<fake></fake>\n";
	    } else {
	      print $INPUTALLXML $linesXMLtrans."\n";
	    }
	    print $INPUTALLXML "</TRANSLOCATIONS>"."\n";
	  }
	}
	close(INPUTALLOUT);
	if (!$doNotPrintXML) {
	  close($INPUTALLXML);
	}
	if (!$doNotPrintGenesOut) {
	  close(GENESALLOUT);
	}
	close(INPUTNONEOUT);

      }				#for each range



    } else {

      ############################################
      ##                                        ##
      ##         END ANNOTATION MODE 5          ##
      ##                                        ##
      ############################################

      print STDERR "invalid mode\n";
      terminate(9);
    }


  }				#end for each file to annotate
  print "done\nProgram ran successfully\n";

  ############################################
  ##                                        ##
  ##                                        ##
  ##             END ANNOTATION             ##
  ##                                        ##
  ##                                        ##
  ############################################

  #print "Time to annotate: ".((time - $startTIME)*1000)." \n";


  terminate(9);

  ################################################
  ################################################
  ##                                            ##
  ##                                            ##
  ##            END MAIN PROGRAM                ##
  ##                                            ##
  ##                                            ##
  ################################################
  ################################################






  #######################################################
  #######################################################
  ####                                               ####
  ####                                               ####
  ####             END ANNOTATION MODE               ####
  ####                                               ####
  ####                                               ####
  #######################################################
  #######################################################

}





















