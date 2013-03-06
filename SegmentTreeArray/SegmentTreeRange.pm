=head1 NAME

   SegmentTree::SegmentTreeRange.pm


=head1 Disclaimer

This software and all the data therein were developed by the Bioinformatics and Computational Biology Lab (Laboratorio de Bioinformatica e Biologia Computational (LBBC)) of the Brazilian National Cancer Institute (Instituto Nacional de Cancer (INCA)) in Rio de Janeiro, RJ, Brazil. The INCA is part of the Health Ministry (Ministerio da Saude) of the Federal Government of Brazil. This software was developed by the staff of the LBBC as part of their official duties. According to federal regulations, this software is therefore public domain.

This software is provided as is and the INCA does make any guarantee regarding the accuracy/quality/reliability of the results. The INCA cannot be held accountable for any hardware or software damage that might be caused through the use or misuse of this software.

The LBBC encourages researchers using our software to acknowledge us. The software can be distributed freely and modified given that any derivative bears an acknowledgment that it was derived from our original software.

=head1 Author

Gabriel Renaud grenaud@inca.gov.br



=cut



package SegmentTreeRange;


use strict;
use warnings;

use Data::Dumper;

sub new{
  my ($class,%arg)=(@_);
  my ($self) =bless{_start  => $arg{start},
                    _end    => $arg{end}}, $class;

  return $self;
}




sub getStart  {
  my($self)=@_;
  return $self->{_start};
}

sub getEnd  {
  my($self)=@_;
  return $self->{_end};
}

sub print  {
  my($self)=@_;
  return "SegmentTreeRange.pm ".$self->{_start}."-".$self->{_end};
}





#returns 1 if otherRange is contained in self
# otherRange    |-----------|
# self       |-----------------|
sub contains  {
  my($self,$otherRange)=@_;
  #  print Dumper($self);
  #  print Dumper($otherRange);
  #  die;
  if( ($self->{_start} <= $otherRange->{_start} ) &&
      ($self->{_end}   >= $otherRange->{_end} ) ){
    return 1;
  }else{
    return 0;
  }
}


#returns 1 if otherRange is overlaps in self
# otherRange                    |-----------|
# self         |-------------|
#
# otherRange    |-----------|
# self                         |--------------|

sub overlaps  {
  my($self,$otherRange)=@_;
  if( ($self->{_end}    < $otherRange->{_start} ) ||
      ($self->{_start}  > $otherRange->{_end} ) ){
    return 0;
  }else{
    return 1;
  }
}


# self                         |--------------|
# coord                               |
sub containsCoord  {
  my($self,$coord)=@_;
  if( ($self->{_start}    <= $coord ) &&
      ($self->{_end}      >= $coord ) ){
    return 1;
  }else{
    return 0;
  }
}


1;
