=head1 NAME

   FastaRecord.pm


=head1 Disclaimer

This software and all the data therein were developed by the Bioinformatics and Computational Biology Lab (Laboratorio de Bioinformatica e Biologia Computational (LBBC)) of the Brazilian National Cancer Institute (Instituto Nacional de Cancer (INCA)) in Rio de Janeiro, RJ, Brazil. The INCA is part of the Health Ministry (Ministerio da Saude) of the Federal Government of Brazil. This software was developed by the staff of the LBBC as part of their official duties. According to federal regulations, this software is therefore public domain.

This software is provided as is and the INCA does make any guarantee regarding the accuracy/quality/reliability of the results. The INCA cannot be held accountable for any hardware or software damage that might be caused through the use or misuse of this software.

The LBBC encourages researchers using our software to acknowledge us. The software can be distributed freely and modified given that any derivative bears an acknowledgment that it was derived from our original software.

=head1 Author

Gabriel Renaud grenaud@inca.gov.br



=cut

package FastaRecord;

use strict;
use warnings;

use Data::Dumper;

sub new{
  my ($class,%arg)=(@_);
  my ($self) =bless{_defline  => $arg{defline},
		    _sequence  => $arg{sequence}}, $class;


  return $self;
}

sub getDefline{
  my($self)=@_;
  return $self->{_defline};
}

sub setDefline{
  my($self,$defline)=@_;
  $self->{_defline}=$defline;
}



sub getSequence{
  my($self)=@_;
  return $self->{_sequence};
}

sub print{
  my($self)=@_;
  return formatFasta($self->{_defline},$self->{_sequence});
}

sub printWithLength{
  my($self)=@_;
  return formatFasta($self->{_defline}."::".length($self->{_sequence}),$self->{_sequence});
}



sub formatFasta{
  my($defline,$sequence)=@_;
  my $stringToReturn;
  $stringToReturn=$defline."\n";
  my $numberOfChar=0;
  foreach my $char (split(//,$sequence)){
    $stringToReturn.=$char;
    $numberOfChar++;
    if($numberOfChar == 60){
      $stringToReturn.="\n";
      $numberOfChar=0;
    }
  }
  $stringToReturn.="\n";

  return $stringToReturn;
}



1;
