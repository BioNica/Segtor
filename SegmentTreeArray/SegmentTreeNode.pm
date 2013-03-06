=head1 NAME

   SegmentTree::SegmentTreeNode.pm


=head1 Disclaimer

This software and all the data therein were developed by the Bioinformatics and Computational Biology Lab (Laboratorio de Bioinformatica e Biologia Computational (LBBC)) of the Brazilian National Cancer Institute (Instituto Nacional de Cancer (INCA)) in Rio de Janeiro, RJ, Brazil. The INCA is part of the Health Ministry (Ministerio da Saude) of the Federal Government of Brazil. This software was developed by the staff of the LBBC as part of their official duties. According to federal regulations, this software is therefore public domain.

This software is provided as is and the INCA does make any guarantee regarding the accuracy/quality/reliability of the results. The INCA cannot be held accountable for any hardware or software damage that might be caused through the use or misuse of this software.

The LBBC encourages researchers using our software to acknowledge us. The software can be distributed freely and modified given that any derivative bears an acknowledgment that it was derived from our original software.

=head1 Author

Gabriel Renaud grenaud@inca.gov.br



=cut



package SegmentTreeNode;


use strict;
use warnings;

use SegmentTreeRange;
use Data::Dumper;





=head1 NAME



=head1 SYNOPSIS



=head1 DESCRIPTION



=head1 AUTHOR

Gabriel Renaud grenaud@inca.gov.br

=cut





=item new ()

  Constructor to build the node. The constructor 
  does not take any parameters

  It contains the following variables:

    _key         Value of the key for the node
    _interval    SegmentTreeRange object corresponding to the interval
    _leftChild   Pointer to a SegmentTreeNode representing the left child
    _rightChild  Pointer to a SegmentTreeNode representing the right child
    _leaf        binary flag to determine if the node is a leaf
    _records     Array of records contained in the node


=cut
sub new{
  my ($class,%arg)=(@_);
  my ($self) =bless{#_key         => undef,
		    _interval    => undef,
		    _leftChild   => undef,
		    _rightChild  => undef,
		    _parent      => undef,
		    _leaf        => undef,
		    _records     => []   ,
		    _root        => 0   }, $class;
  return $self;
}



=item getInterval ()

  Subroutine to retrieve the _interval parameter

=cut
sub getInterval{
  my($self)=@_;
  return $self->{_interval};
}


=item getLeftChild ()

  Subroutine to retrieve the _leftChild parameter

=cut
sub getLeftChild{
  my($self)=@_;
  return $self->{_leftChild};
}

=item getRightChild ()

  Subroutine to retrieve the _rightChild parameter

=cut
sub getRightChild{
  my($self)=@_;
  return $self->{_rightChild};
}

=item getLeaf ()

  Subroutine to retrieve the _leaf parameter

=cut
sub getLeaf{
  my($self)=@_;
  return $self->{_leaf};
}

=item getRecords ()

  Subroutine to retrieve the _records parameter

=cut
sub getRecords{
  my($self)=@_;
  return $self->{_records};
}


=item getParent ()

  Subroutine to retrieve the parent of the node

=cut
sub getParent{
  my($self)=@_;
  return $self->{_parent};
}


=item isRoot ()

  Subroutine to determine if the node is the root or not

=cut
sub isRoot{
  my($self)=@_;
  return $self->{_root};
}


=item setAsRoot ()

  Subroutine to set the root flag

=cut
sub setAsRoot{
  my($self)=@_;
  $self->{_root}=1;
}










=item setInterval ()

  Subroutine to set the _interval parameter

@param $interval  SegmentTreeRange object corresponding to the interval

=cut
sub setInterval{
  my($self,$interval)=@_;
  $self->{_interval}= $interval  ;
}

=item setLeftChild ()

  Subroutine to set the _leftChild parameter

@param $leftChild Pointer to a SegmentTreeNode representing the left child

=cut
sub setLeftChild{
  my($self,$leftChild)=@_;
  $self->{_leftChild}= $leftChild  ;
}

=item setRightChild ()

  Subroutine to set the _rightChild parameter

@param $rightChild Pointer to a SegmentTreeNode representing the right child

=cut
sub setRightChild{
  my($self,$rightChild)=@_;
  $self->{_rightChild}= $rightChild   ;
}




=item setLeaf ()

  Subroutine to set the _leaf parameter

@param $leaf  binary flag to determine if the node is a leaf

=cut
sub setLeaf{
  my($self,$leaf)=@_;
  $self->{_leaf}= $leaf  ;
}



=item addRecord ()

  Subroutine to set the _records parameter

@param $record Array of records contained in the node

=cut
sub addRecord {
  my($self,$record)=@_;
  push(@{$self->{_records}},$record);
}


=item print ()

  Subroutine to return a string representation of the node

=cut
sub print{
  my($self)=@_;
  return "key ".$self->{_interval}->print();
}


=item setParent ()

  Subroutine to set the parent for the current node

@param $parent Parent node for the given node

=cut
sub setParent {
  my($self,$parent)=@_;
  $self->{_parent}= $parent ;

}




1;
