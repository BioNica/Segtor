=head1 NAME

   SegmentTree::SegmentTree.pm


=head1 Disclaimer

This software and all the data therein were developed by the Bioinformatics and Computational Biology Lab (Laboratorio de Bioinformatica e Biologia Computational (LBBC)) of the Brazilian National Cancer Institute (Instituto Nacional de Cancer (INCA)) in Rio de Janeiro, RJ, Brazil. The INCA is part of the Health Ministry (Ministerio da Saude) of the Federal Government of Brazil. This software was developed by the staff of the LBBC as part of their official duties. According to federal regulations, this software is therefore public domain.

This software is provided as is and the INCA does make any guarantee regarding the accuracy/quality/reliability of the results. The INCA cannot be held accountable for any hardware or software damage that might be caused through the use or misuse of this software.

The LBBC encourages researchers using our software to acknowledge us. The software can be distributed freely and modified given that any derivative bears an acknowledgment that it was derived from our original software.

=head1 Author

Gabriel Renaud grenaud@inca.gov.br



=cut

package SegmentTree;

use strict;
use warnings;

use SegmentTreeRange;
use SegmentTreeNode;

use Data::Dumper;

use POSIX qw(ceil floor);

=head1 NAME



=head1 SYNOPSIS


=head1 DESCRIPTION



=head1 AUTHOR

Gabriel Renaud grenaud@inca.gov.br

=cut


my $MAXINT=2147483647;
my $MAXFORARRAY=50;


=item createSubTreeWithArray

  Subroutine used to build the basic tree using the sorted endpoints
  in $array. The subroutine recursively divides the $array in two and stops
  when a single element is found. The single elements becomes a leaf while
  the subroutine calls that generate recursive call generate internal nodes.

  Input :                      1 2 3 4 5

  createSubTreeWithArray :   [ 1 2 3 4 5 ]
  createSubTreeWithArray :   [ 1 2 3 ] [ 4 5 ]
  createSubTreeWithArray :   [ [ 1 2 ] [ 3 ] ]      [ [ 4 ] [ 5 ] ]
  createSubTreeWithArray :   [ [ [ 1 ] [ 2 ] ] [ 3 ] ]      [ [ 4 ] [ 5 ] ]

                                    []                   []
                                   / \  [3]             / \
                                   1 2                  4 5

                                    []
                                       \               []
                                   / \  3             / \
                                   1 2                4 5


                                             []
                                            /  \
                                           /    \
                                          []     []
                                            \    / \
                                        / \  3   4 5
                                        1 2



@param  $array Array containing the endpoints
@param  $index1 left  index on which to launch the subroutine
@param  $index2 right index on which to launch the subroutine


=cut
sub createSubTreeWithArray{
  my ($array,$index1,$index2)=@_;

  if($index1 == $index2){
    die "Big problem createSubTreeWithArray\n";
  }elsif($index1 == ($index2-1)){
    my $node = new SegmentTreeNode();
    #    $node->setKey($array->[$index1]);
    $node->setLeaf(1);
    $node->setInterval(new SegmentTreeRange(start => $array->[$index1],end => $array->[$index2]) );
    return $node;
  }else{
    my $node = new SegmentTreeNode();
    my $middleIndex=floor( ($index2+$index1)/2 );

    $node->setLeaf(0);
    my $leftSubtree=createSubTreeWithArray($array,$index1,$middleIndex);
    my $rightSubtree=createSubTreeWithArray($array,$middleIndex,$index2);

    $node->setLeftChild(  $leftSubtree  );
    $node->setRightChild( $rightSubtree );
    $leftSubtree->setParent($node);
    $rightSubtree->setParent($node);

    $node->setInterval(new SegmentTreeRange(start => $leftSubtree->getInterval()->getStart(),end => $rightSubtree->getInterval()->getEnd() ) );
    return $node;
  }

}


=item insertInTreeRecursive

  Subroutine used to recursively add the record to insert $recordToInsert
  into the node once the createSubTreeWithArray() established the nodes
  given the endpoints. If the range in $recordToInsert does not perfectly
  contains the interval of the node, we will launch recursively on the
  child nodes the range in $recordToInsert overlaps the interval of 
  the child nodes

@param  $recordToInsert The record to insert
@param  $node  The current node


=cut
sub insertInTreeRecursive{
  my ($recordToInsert,$node)=@_;

  #If the range perfectly contains the interval of the $node
  if($recordToInsert->{'range'}->contains( $node->getInterval() )){
    $node->addRecord($recordToInsert);
  }else{ #if not, launch recursively on the child nodes
    if($node->getLeftChild()){
      if( $recordToInsert->{'range'}->overlaps( $node->getLeftChild()->getInterval() ) ){
	insertInTreeRecursive($recordToInsert,$node->getLeftChild());
      }
    }

    if($node->getRightChild()){
      if( $recordToInsert->{'range'}->overlaps( $node->getRightChild()->getInterval() ) ){
	insertInTreeRecursive($recordToInsert,$node->getRightChild());
      }
    }
  }

}


=item new ()

  Constructor to build the tree
@param arrayOfInputs Array of inputs where each contain the 'id' and the 'range'

=cut
sub new{
  my ($class,%arg)=(@_);
  my ($self) =bless{_arrayOfInputs  => $arg{arrayOfInputs}}, $class;

  if ( $#{$self->{_arrayOfInputs}} <= $MAXFORARRAY ) {
    $self->{_useOfArray}=1;
  } else {
    $self->{_useOfArray}=0;

    my $arrayOfStEndCoords;
    foreach my $record (@{$self->{_arrayOfInputs}}) {
      push(@{$arrayOfStEndCoords},$record->{'range'}->getStart());
      push(@{$arrayOfStEndCoords},$record->{'range'}->getEnd());
    }

    my %hashTemp  =  map { $_ => 1 } @${arrayOfStEndCoords};
    my @arrayTemp =  sort {$a <=> $b} (keys %hashTemp);
    unshift(@arrayTemp,-$MAXINT);
    push(@arrayTemp,$MAXINT);

    $self->{_root}=createSubTreeWithArray(\@arrayTemp,0,$#arrayTemp);
    $self->{_root}->setAsRoot();
    foreach my $record (@{$self->{_arrayOfInputs}}) {
      insertInTreeRecursive($record,$self->{_root});
    }
  }
  return $self;
}



=item seekCoordRecursive ()

  Subroutine used to recursively find the ranges that intersect the coord
  by putting the records in the current node in $arrayRecords and launching
  a recursion on the child nodes if their intervals intersect the $coord

@param  $node  The current node
@param  $coord The coordinate to use
@param  $arrayRecords The array of records that was found thus far.

@return The array of records overlapping the $coord

=cut
sub seekCoordRecursive{
  my ($node,$coord,$arrayRecords)=@_;

  push(@{$arrayRecords},@{$node->getRecords()});
  if(! $node->getLeaf()){
    if($node->getLeftChild()->getInterval()->containsCoord($coord)){
      seekCoordRecursive($node->getLeftChild(),$coord,$arrayRecords);
    }

    if($node->getRightChild()->getInterval()->containsCoord($coord)){
      seekCoordRecursive($node->getRightChild(),$coord,$arrayRecords);
    }

  }
}


=item seekCoord ()

  Subroutine used to find the ranges that intersect the coord
  by calling seekCoordRecursive() on the root
@param  $coord Coordinate to look for
@return The array of records overlapping the $coord

=cut
sub seekCoord{
  my ($self,$coord)=@_;

  if($coord < (-$MAXINT) ){
    die "Error SegmentTree.pm seekCoord, the coordinate $coord cannot be lesser than -$MAXINT\n";
  }

  if($coord > $MAXINT ){
    die "Error SegmentTree.pm seekCoord, the coordinate $coord cannot be greater than $MAXINT\n";
  }

  my $arrayToReturn=[];
  if( $self->{_useOfArray} ){
    foreach my $record (@{$self->{_arrayOfInputs}}) {
      if($record->{'range'}->containsCoord($coord)){
	push(@{$arrayToReturn},$record);
      }
    }
  }else{
    if( $self->{_root}->getInterval()->containsCoord($coord) ){
      seekCoordRecursive($self->{_root},$coord,$arrayToReturn);
    }
  }
  return $arrayToReturn;
}


=item seekRange ()

  Subroutine used to find the ranges that intersect the range
  provided as argument using the start and end coordinates
@param  $startCoord Start coordinate to look for
@param  $endCoord   End coordinate to look for
@return The array of records overlapping the $coord

=cut
sub seekRange{
  my ($self,$startCoord,$endCoord)=@_;

  if($startCoord < (-$MAXINT) ){
    die "Error SegmentTree.pm seekRange, the $startCoord coordinate cannot be lesser than -$MAXINT\n";
  }
  if($endCoord > $MAXINT){
    die "Error SegmentTree.pm seekRange, the $endCoord coordinate cannot be greater than $MAXINT\n";
  }

  my $arrayToReturn=[];
  my $range=new SegmentTreeRange(start => $startCoord, end   => $endCoord);
  if( $self->{_useOfArray} ){

    foreach my $record (@{$self->{_arrayOfInputs}}) {
      if($record->{'range'}->overlaps($range)){
	push(@{$arrayToReturn},$record);
      }
    }

  } else {

    if ( $self->{_root}->getInterval()->contains($range) ) { #the tree contains the entire range
      seekRangeRecursive($self->{_root},$range,$arrayToReturn);

    } elsif ( $self->{_root}->getInterval()->containsCoord($startCoord) &&
	      !($self->{_root}->getInterval()->containsCoord($endCoord)) ) { #the tree contains the start coordinate only
      seekRangeRecursiveLeftTree($self->{_root},$startCoord,$arrayToReturn);

    } elsif ( !($self->{_root}->getInterval()->containsCoord($startCoord)) &&
	      ($self->{_root}->getInterval()->containsCoord($endCoord)) ) { #the tree contains the end coordinate only
      seekRangeRecursiveRightTree($self->{_root},$endCoord,$arrayToReturn);

    } elsif ( !($self->{_root}->getInterval()->containsCoord($startCoord)) &&
	      !($self->{_root}->getInterval()->containsCoord($endCoord)) ) { #the tree does not overlap the range provided
      #nothing to do
    } else {
      die "Invalid state in seekRange()\n";
    }
  }

  return $arrayToReturn;
}


=item seekRangeRecursive ()

  Subroutine used to recursively find the ranges
  that intersect a given range. Called by seekRange()
@param  $node The current node
@param  $rangeCoord   The SegmentTreeRange.pm object 
@param  $arrayToReturn   Reference to the array used to store the results

=cut
sub seekRangeRecursive{
  my ($node,$rangeCoord,$arrayRecords)=@_;


  push(@{$arrayRecords},@{$node->getRecords()});

  if (! $node->getLeaf()) {
    if ($node->getLeftChild()->getInterval()->contains($rangeCoord)) {
      seekRangeRecursive($node->getLeftChild(),$rangeCoord,$arrayRecords);
    } elsif ($node->getRightChild()->getInterval()->contains($rangeCoord)) {
      seekRangeRecursive($node->getRightChild(),$rangeCoord,$arrayRecords);
    } elsif ( $node->getLeftChild()->getInterval()->containsCoord( $rangeCoord->getStart() ) &&
	      $node->getRightChild()->getInterval()->containsCoord( $rangeCoord->getEnd() ) ) {
      seekRangeRecursiveLeftTree(  $node->getLeftChild(),  $rangeCoord->getStart(), $arrayRecords);
      seekRangeRecursiveRightTree( $node->getRightChild(), $rangeCoord->getEnd(), $arrayRecords);
    } else {
      die "Invalid state in seekRangeRecursive()\n";
    }
  }

}


=item seekRangeRecursiveLeftTree ()

  Subroutine used to recursively find the ranges
  that intersect a given range with the lowest coordinate.
  Called by seekRangeRecursive()
@param  $node The current node
@param  $coord   The coordinate to find
@param  $arrayToReturn   Reference to the array used to store the results

=cut
sub seekRangeRecursiveLeftTree{
  my ($node,$coord,$arrayRecords)=@_;

  push(@{$arrayRecords},@{$node->getRecords()});

  if(! $node->getLeaf()){
    if($node->getLeftChild()->getInterval()->containsCoord($coord)){
      seekRangeRecursiveLeftTree($node->getLeftChild(),$coord,$arrayRecords);
      returnAllInTree($node->getRightChild(),$arrayRecords);
    }elsif($node->getRightChild()->getInterval()->containsCoord($coord)){
      seekRangeRecursiveLeftTree($node->getRightChild(),$coord,$arrayRecords);
    }else{
      die "Invalid state in seekRangeRecursiveLeftTree()\n";
    }
  }
}


=item seekRangeRecursiveLeftTree ()

  Subroutine used to recursively find the ranges
  that intersect a given range with the highest coordinate.
  Called by seekRangeRecursive()
@param  $node The current node
@param  $coord   The coordinate to find
@param  $arrayToReturn   Reference to the array used to store the results

=cut
sub seekRangeRecursiveRightTree{
  my ($node,$coord,$arrayRecords)=@_;

  push(@{$arrayRecords},@{$node->getRecords()});

  if(! $node->getLeaf()){
    if($node->getLeftChild()->getInterval()->containsCoord($coord)){
      seekRangeRecursiveRightTree($node->getLeftChild(),$coord,$arrayRecords);
    }elsif($node->getRightChild()->getInterval()->containsCoord($coord)){
      seekRangeRecursiveRightTree($node->getRightChild(),$coord,$arrayRecords);
      returnAllInTree($node->getLeftChild(),$arrayRecords);
    }else{
      die "Invalid state in seekRangeRecursiveLeftTree()\n";
    }
  }
}


=item returnAllInTree ()

  Subroutine used to recursively return the records
  contained in a subtree.
  Called by seekRangeRecursiveLeftTree() or seekRangeRecursiveRightTree()
@param  $node            The current node
@param  $arrayToReturn   Reference to the array used to store the results

=cut
sub returnAllInTree{
  my ($node,$arrayRecords)=@_;

  push(@{$arrayRecords},@{$node->getRecords()});

  if(! $node->getLeaf()){
    returnAllInTree($node->getLeftChild(),$arrayRecords);
    returnAllInTree($node->getRightChild(),$arrayRecords);
  }
}


1;
