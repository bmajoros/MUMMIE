package Numbers;
use strict;

######################################################################
#
# Numbers.pm bmajoros
#
# 
# 
#
# Attributes:
#
# Methods:
#   $numbers=new Numbers();
#   $exponent=Numbers::getExponent(1e-4 or 0.006);
#   $sci=Numbers::scientificNotation($value) -> 1e-8
#   $printable=Numbers::addCommas(1827317263) -> 1,827,317,263
#   
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
sub new
{
  my ($class)=@_;
  
  my $self={};
  bless $self,$class;

  return $self;
}
#---------------------------------------------------------------------
#   $sci=Numbers::scientificNotation($value) -> 1e-8
sub scientificNotation
  {
    my ($value)=@_;
    my $exponent=getExponent($value);
    return "1e$exponent";
  }
#---------------------------------------------------------------------
#   $exponent=Numbers::getExponent(1e-4 or 0.006);
sub getExponent
  {
    my ($P)=@_;
    my $logP=int(log($P)/log(10)-0.5);
    return $logP;
  }
#---------------------------------------------------------------------
#   $printable=Numbers::addCommas(1827317263) -> 1,827,317,263
sub addCommas
  {
    my ($x)=@_;
    my $r="";
    while($x=~/(.+)(...)\s*$/)
      {
	my ($left,$right)=($1,$2);
	$r=",$right$r";
	$x=$left;
      }
    $r="$x$r";
    return $r;
  }
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------






#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------

1;

