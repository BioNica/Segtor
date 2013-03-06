=head1 NAME

   ChromosomeIndexing::ChromosomeIndexing.pm


=head1 Disclaimer

This software and all the data therein were developed by the Bioinformatics and Computational Biology Lab (Laboratorio de Bioinformatica e Biologia Computational (LBBC)) of the Brazilian National Cancer Institute (Instituto Nacional de Cancer (INCA)) in Rio de Janeiro, RJ, Brazil. The INCA is part of the Health Ministry (Ministerio da Saude) of the Federal Government of Brazil. This software was developed by the staff of the LBBC as part of their official duties. According to federal regulations, this software is therefore public domain.

This software is provided as is and the INCA does make any guarantee regarding the accuracy/quality/reliability of the results. The INCA cannot be held accountable for any hardware or software damage that might be caused through the use or misuse of this software.

The LBBC encourages researchers using our software to acknowledge us. The software can be distributed freely and modified given that any derivative bears an acknowledgment that it was derived from our original software.

=head1 Author

Gabriel Renaud grenaud@inca.gov.br



=cut


package ChromosomeIndexing;


use strict;
use warnings;


use Data::Dumper;


my $hashSpecies2Indexing={};
my $baseDir="/storage/data/ucsc/";






=item getBaseDir

   Subroutine to return the base directory
   where the UCSC data is located

=cut
sub getBaseDir{
  return $baseDir;
}



=item setBaseDir

   Subroutine to set the base directory
   where the UCSC data is located

=cut
sub setBaseDir{
  my ($baseDirToUse)=@_;
  $baseDir=$baseDirToUse;
}






=item indexExistsForSpecies

   Subroutine to determine if the index file required by this program exists or not
   Returns 1 if it exists, 0 otherwise

=cut
sub indexExistsForSpecies{
  my ($species)=@_;
  if( (-e $baseDir.$species."/chromosomes/index.txt" ) ){
    return 1;
  }else{
    return 0;
  }
}




=item speciesWasIndexed

   Subroutine to determine if a species has been already indexed in memory by the program.
   Returns 1 if it was, 0 otherwise

=cut
sub speciesWasIndexed{
  my ($species)=@_;
  if(exists $hashSpecies2Indexing->{$species}){
    return 1;
  }else{
    return 0;
  }
}




=item indexSpecies

   Subroutine to create the hash structure hashSpecies2Indexing
   based on the species code given by : $speciesToIndex

   The hash hashSpecies2Indexing contains the main data
   structure that works by having the species symbol as key
   which points to 6 different hashes:

   $hashSpecies2Indexing={'human-hg18'}={'file2chr'    =  'scaffolds.fa'  => 'scafford3294'
                                         'chr2file'    =  'scafford3294'  => 'scaffolds.fa'
                                         'file2index'  =  'scaffolds.fa'  => '1'
                                         'index2file'  =  '1'             => 'scaffolds.fa'
                                         'chr2index'   =  'scafford3294'  => '12'
                                         'index2chr'   =  '12'            => 'scafford3294'
                                         };


=cut

sub indexSpecies{
  my ($speciesToIndex)=@_;

  if(exists $hashSpecies2Indexing->{$speciesToIndex}){
    die "ChromosomeIndexing.pm Species $speciesToIndex exists\n";
  }

  $hashSpecies2Indexing->{$speciesToIndex}={};

  $hashSpecies2Indexing->{$speciesToIndex}->{'file2chr'}={};
  $hashSpecies2Indexing->{$speciesToIndex}->{'chr2file'}={};

  $hashSpecies2Indexing->{$speciesToIndex}->{'file2index'}={};
  $hashSpecies2Indexing->{$speciesToIndex}->{'index2file'}={};

  $hashSpecies2Indexing->{$speciesToIndex}->{'chr2index'}={};
  $hashSpecies2Indexing->{$speciesToIndex}->{'index2chr'}={};


  if(! (-e $baseDir.$speciesToIndex."/chromosomes/index.txt" ) ){
    die "Error, unable to open index file ".$baseDir.$speciesToIndex."/chromosomes/index.txt\n";
  }

  my $indexFile=$baseDir.$speciesToIndex."/chromosomes/index.txt";
  open(INFO, $indexFile) or die ("Can't open $indexFile");

  my $currentFile="#";
  my $currentFasta=[];

  my $fileIndex=0;
  my $fastaRecordIndex=0;

  while(my $line = <INFO>){





    if($line =~ /^(\w\S+)\.fa$/){
      my $file=$1;

      if ($currentFile ne "#") {
	#Key = file value = anonymous array of chr
	if (exists $hashSpecies2Indexing->{$speciesToIndex}->{'file2chr'}->{$currentFile}) {
	  die "Problem #7: duplicate keys in file2chr\n";
	} else {
	  $hashSpecies2Indexing->{$speciesToIndex}->{'file2chr'}->{$currentFile}=$currentFasta;
	}
      }

      $currentFasta=[];
      $currentFile=$file;
      $fileIndex++;
      #print "$file $fileIndex\n";

      #Key = file value = index of file
      if(exists $hashSpecies2Indexing->{$speciesToIndex}->{'file2index'}->{$file}){
	die "Problem #5: duplicate keys in file2index\n";
      }else{
	$hashSpecies2Indexing->{$speciesToIndex}->{'file2index'}->{$file}=$fileIndex;
      }

      #Key = index of file value = file
      if(exists $hashSpecies2Indexing->{$speciesToIndex}->{'index2file'}->{$fileIndex}){
	die "Problem #6: duplicate keys in index2file\n";
      }else{
	$hashSpecies2Indexing->{$speciesToIndex}->{'index2file'}->{$fileIndex}=$file;
      }









    }elsif($line =~ /^>(\w\S+)$/){
      my $chr=$1;
      push(@{$currentFasta},$chr);
      $fastaRecordIndex++;
      #print "\t$1 $fastaRecordIndex\n";

      #Key = fasta record value = index of chromosome
      if(exists $hashSpecies2Indexing->{$speciesToIndex}->{'chr2index'}->{$chr}){
	die "Problem #1: duplicate $chr keys in chr2index\n";
      }else{
	$hashSpecies2Indexing->{$speciesToIndex}->{'chr2index'}->{$chr}=$fastaRecordIndex;
      }

      #Key =  index of chromosome value = fasta record
      if(exists $hashSpecies2Indexing->{$speciesToIndex}->{'index2chr'}->{$fastaRecordIndex}){
	die "Problem #2: duplicate keys in index2chr\n";
      }else{
	$hashSpecies2Indexing->{$speciesToIndex}->{'index2chr'}->{$fastaRecordIndex}=$chr;
      }

      if ($currentFasta eq "#") {
	die "Problem #3: The currentFasta variable was not set\n";
      } else {
	#Key = fasta record value = file containing the fasta record
	if (exists $hashSpecies2Indexing->{$speciesToIndex}->{'chr2file'}->{$chr}) {
	  die "Problem #4: duplicate keys in chr2file\n";
	} else {
	  $hashSpecies2Indexing->{$speciesToIndex}->{'chr2file'}->{$chr}=$currentFile;
	}
      }


    }else{
      die "Wrong line $line in index file $indexFile\n";
    }
  }

  close(INFO);


  if ($currentFile ne "#") {
    #Key = file value = anonymous array of chr
    if (exists $hashSpecies2Indexing->{$speciesToIndex}->{'file2chr'}->{$currentFile}) {
      die "Problem #7: duplicate keys in file2chr\n";
    } else {
      $hashSpecies2Indexing->{$speciesToIndex}->{'file2chr'}->{$currentFile}=$currentFasta;
    }
  }else{
    die "File $indexFile is empty\n";
  }

  #die Dumper($hashSpecies2Indexing);
}




=item file2chr

  Subroutine returns an anonymous array
  of chromosomes for the given file

=cut
sub file2chr{
  my ($species,$file)=@_;

  if(!speciesWasIndexed($species)){
    indexSpecies($species);
  }

  if (exists $hashSpecies2Indexing->{$species}->{'file2chr'}->{$file}) {
    return $hashSpecies2Indexing->{$species}->{'file2chr'}->{$file};
  }else{
    die "ChromosomeIndexing.pm: Key $file was not found in file2chr\n";
  }
}


=item file2chrArray

  Subroutine returns an array
  of chromosomes for the given file

=cut
sub file2chrArray{
  my ($species,$file)=@_;
  my @arrayToReturn;

  if(!speciesWasIndexed($species)){
    indexSpecies($species);
  }

  if (exists $hashSpecies2Indexing->{$species}->{'file2chr'}->{$file}) {
    @arrayToReturn = $hashSpecies2Indexing->{$species}->{'file2chr'}->{$file};
  }else{
    die "ChromosomeIndexing.pm: Key $file was not found in file2chr\n";
  }

  return @arrayToReturn;
}



=item chr2file

  Subroutine to return the file
  which contains the given chromosome

=cut
sub chr2file{
  my ($species,$chr)=@_;

  if(!speciesWasIndexed($species)){
    indexSpecies($species);
  }

  if (exists $hashSpecies2Indexing->{$species}->{'chr2file'}->{$chr}) {
    return $hashSpecies2Indexing->{$species}->{'chr2file'}->{$chr};
  }else{
    die "ChromosomeIndexing.pm: Key $chr was not found in chr2file\n";
  }
}

=item file2index

  Subroutine to return the index
  for the given file

=cut
sub file2index{
  my ($species,$file)=@_;

  if(!speciesWasIndexed($species)){
    indexSpecies($species);
  }

  if (exists $hashSpecies2Indexing->{$species}->{'file2index'}->{$file}) {
    return $hashSpecies2Indexing->{$species}->{'file2index'}->{$file};
  }else{
    die "ChromosomeIndexing.pm: Key $file was not found in file2index\n";
  }
}


=item index2file

  Subroutine to return the file
  for the given index

=cut
sub index2file{
  my ($species,$index)=@_;

  if(!speciesWasIndexed($species)){
    indexSpecies($species);
  }

  if (exists $hashSpecies2Indexing->{$species}->{'index2file'}->{$index}) {
    return $hashSpecies2Indexing->{$species}->{'index2file'}->{$index};
  }else{
    die "ChromosomeIndexing.pm: Key $index was not found in index2file\n";
  }

}



=item chr2index

  Subroutine to return the index
  for the given chromosome

=cut
sub chr2index{
  my ($species,$chr)=@_;

  if(!speciesWasIndexed($species)){
    indexSpecies($species);
  }

  if (exists $hashSpecies2Indexing->{$species}->{'chr2index'}->{$chr}) {
    return $hashSpecies2Indexing->{$species}->{'chr2index'}->{$chr};
  }else{
    die "ChromosomeIndexing.pm: Key $chr was not found in chr2index\n";
  }
}



=item index2chr

  Subroutine to return the chromosome
  for the given index

=cut
sub index2chr{
  my ($species,$index)=@_;

  if(!speciesWasIndexed($species)){
    indexSpecies($species);
  }

  if (exists $hashSpecies2Indexing->{$species}->{'index2chr'}->{$index}) {
    return $hashSpecies2Indexing->{$species}->{'index2chr'}->{$index};
  }else{
    die "ChromosomeIndexing.pm: Key $index was not found in index2chr\n";
  }

}






















=item fileEXISTS

  Subroutine to verify if a file is part of a given
  build $species. Returns 1 if so, 0 otherwise

=cut
sub fileEXISTS{
  my ($species,$file)=@_;

  if(!speciesWasIndexed($species)){
    indexSpecies($species);
  }

  if (exists $hashSpecies2Indexing->{$species}->{'file2index'}->{$file}) {
    return 1;
  }else{
    return 0;
  }

}



=item chrEXISTS

  Subroutine to verify if a chr index is part of a given
  build $species. Returns 1 if so, 0 otherwise

=cut
sub chrEXISTS{
  my ($species,$chr)=@_;

  if(!speciesWasIndexed($species)){
    indexSpecies($species);
  }

  if (exists $hashSpecies2Indexing->{$species}->{'chr2index'}->{$chr}) {
    return 1;
  }else{
    return 0;
  }
}




=item fileIndexEXISTS

  Subroutine to verify if a file index is part of a given
  build $species. Returns 1 if so, 0 otherwise

=cut
sub fileIndexEXISTS{
  my ($species,$index)=@_;

  if(!speciesWasIndexed($species)){
    indexSpecies($species);
  }

  if (exists $hashSpecies2Indexing->{$species}->{'index2file'}->{$index}) {
    return 1;
  }else{
    return 0;
  }

}



=item chrIndexEXISTS

  Subroutine to verify if a chr is part of a given
  build $species. Returns 1 if so, 0 otherwise

=cut
sub chrIndexEXISTS{
  my ($species,$index)=@_;

  if(!speciesWasIndexed($species)){
    indexSpecies($species);
  }

  if (exists $hashSpecies2Indexing->{$species}->{'index2chr'}->{$index}) {
    return 1;
  }else{
    return 0;
  }
}




=item listIndexChr

  Returns the list of all the chromosome
  indices for the given species

=cut
sub listIndexChr{
  my ($species)=@_;

  if(!speciesWasIndexed($species)){
    indexSpecies($species);
  }

  return  sort {$a <=> $b} (keys %{$hashSpecies2Indexing->{$species}->{'index2chr'}});
}




=item listChr

  Returns the list of all the chromosome
  names for the give species

=cut
sub listChr{
  my ($species)=@_;

  if(!speciesWasIndexed($species)){
    indexSpecies($species);
  }
  my @arrayToReturn;
  my @arrayTemp=sort {$a <=> $b} (keys %{$hashSpecies2Indexing->{$species}->{'index2chr'}});
  foreach my $index (@arrayTemp){
    push(@arrayToReturn,$hashSpecies2Indexing->{$species}->{'index2chr'}->{$index});
  }
  return @arrayToReturn;
}




=item listIndexFile

  Returns the list of all the file
  indices for the given species

=cut
sub listIndexFile{
  my ($species)=@_;

  if(!speciesWasIndexed($species)){
    indexSpecies($species);
  }

  return  sort {$a <=> $b} (keys %{$hashSpecies2Indexing->{$species}->{'index2file'}});
}




=item listFile

  Returns the list of all the file
  names for the given species

=cut
sub listFile{
  my ($species)=@_;

  if(!speciesWasIndexed($species)){
    indexSpecies($species);
  }
  my @arrayToReturn;
  my @arrayTemp=sort {$a <=> $b} (keys %{$hashSpecies2Indexing->{$species}->{'index2file'}});
  foreach my $index (@arrayTemp){
    push(@arrayToReturn,$hashSpecies2Indexing->{$species}->{'index2file'}->{$index});
  }
  return @arrayToReturn;
}











