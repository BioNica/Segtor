=head1 NAME

   GeneDataStructure::Range.pm


=head1 Disclaimer

This software and all the data therein were developed by the Bioinformatics and Computational Biology Lab (Laboratorio de Bioinformatica e Biologia Computational (LBBC)) of the Brazilian National Cancer Institute (Instituto Nacional de Cancer (INCA)) in Rio de Janeiro, RJ, Brazil. The INCA is part of the Health Ministry (Ministerio da Saude) of the Federal Government of Brazil. This software was developed by the staff of the LBBC as part of their official duties. According to federal regulations, this software is therefore public domain.

This software is provided as is and the INCA does make any guarantee regarding the accuracy/quality/reliability of the results. The INCA cannot be held accountable for any hardware or software damage that might be caused through the use or misuse of this software.

The LBBC encourages researchers using our software to acknowledge us. The software can be distributed freely and modified given that any derivative bears an acknowledgment that it was derived from our original software.

=head1 Author

Gabriel Renaud grenaud@inca.gov.br



=cut

package Range;

use strict;
use warnings;

use Data::Dumper;

=item new

   Constructor that requires the following parameters:
   * start
   * end

   The start coordinate must be lesser or equal to the end coordinate

   Example of usage:
   my $range = Range->new(start => 5,
                          end   => 15);

=cut
sub new{
  my ($class,%arg)=(@_);

  if($arg{start} !~ /^\d+$/  &&
     $arg{start} !~ /^\d+\.\d+$/){
    die "Range.pm The start coordinate ".$arg{start}." must be numerical\n";
  }

  if($arg{end} !~ /^\d+$/  &&
     $arg{end} !~ /^\d+\.\d+$/){
    die "Range.pm The end coordinate ".$arg{end}." must be numerical\n";
  }


  my ($self) =bless{_start  => $arg{start},
		    _end    => $arg{end}}, $class;
  if($self->{_end} < $self->{_start} ){
    die "Range.pm The end ".$self->{_end}." of the range cannot be lesser than the start ".$self->{_start}."\n";
  }
  return $self;
}


=item getStart

   Retrieve the lowest coordinate (start)

=cut
sub getStart{
  my($self)=@_;
  return $self->{_start};
}

=item getEnd

   Retrieve the highest coordinate (start)

=cut
sub getEnd{
  my($self)=@_;
  return $self->{_end};
}



=item setStart

   Sets the lowest coordinate (start)

=cut
sub setStart{
  my($self,$start)=@_;
  if($self->{_end} < $start ){
    die "Range.pm The end ".$self->{_end}." of the range cannot be lesser than the start ".$start."\n";
  }
  $self->{_start}=$start;
}


=item setEnd

   Sets the highest coordinate (start)

=cut
sub setEnd{
  my($self,$end)=@_;
  if($end < $self->{_start} ){
    die "Range.pm The end ".$end." of the range cannot be lesser than the start ".$self->{_start}."\n";
  }

  $self->{_end}=$end;
}





=item getLength

   Returns the length of the range.

=cut
sub getLength{
  my($self)=@_;
  return ($self->{_end}-($self->{_start}-1));
}

=item print

   Returns a string in the format [start]-[end]

=cut
sub print{
  my($self)=@_;
  return $self->{_start}."-".$self->{_end};
}

=item overlaps

   Determines if the range overlaps
   another range object. Returns 1 if it does, 0 otherwise

=cut

sub overlaps{
  my($self,$otherRange)=@_;
  if( ($otherRange->getEnd()   < $self->{_start} ) ||
      ($otherRange->getStart() > $self->{_end} ) ){
    return 0;
  }else{
    return 1;
  }
}


=item contains

   Determines if the range overlaps
   a given coordinate. Returns 1 if it does, 0 otherwise

=cut
sub contains{
  my($self,$coordinate)=@_;
  if( ($self->{_start} <= $coordinate) &&
      ($coordinate     <= $self->{_end}) ){
    return 1;
  }else{
    return 0;
  }
}








1;
