#!/usr/bin/perl
use strict;
use TempFilename;

my $VERBOSE=1;

my $WANT_BULGE=0;
my $WANT_CONSERVATION=1;
my $TARGETSCAN_VARIANCE=3.0;
my $TARGETSCAN_BG_MEAN=0.0;
my $TARGETSCAN_FG_MEAN=3.0;


die "microMUMMIE.pl <mature-miRNAs.txt> <genome.2bit> <paralyzer-output-dir> <library-name> <out.gff> <posterior-decoding:0/1> <UTRs.txt>\n"
  unless @ARGV==7;
my ($mature,$twoBitFile,$dataDir,$libraryName,$outfile,$wantPost,$UTR_FILE)=@ARGV;
#my $DASH_G= "-g";
my $DASH_G=$wantPost ? "-I" : "-g";

System("date");

print "
# -------------------------------------------
# EXTRACTING SEEDS FROM MATURE MICRO-RNA LIST
# -------------------------------------------
";

System("get-seeds.pl $mature seeds.fasta");

print "
# -------------------------------------------
# BUILDING THE SITE SUBMODEL
# -------------------------------------------
";

System("make-42state-site.pl seeds.fasta .33 .33 .33 4 cons.schema 7 bg.hmm site.hmm '' N N");

print "
# -------------------------------------------
# COMBINING SUBMODELS INTO FULL MODEL
# -------------------------------------------
";

System("model-combiner metamodel.txt submodels.txt PARCLIP.hmm");
if($WANT_CONSERVATION)
  { 
  System("hmm-edit PARCLIP.hmm MEAN all -- 1 $TARGETSCAN_BG_MEAN");
  System("hmm-edit PARCLIP.hmm VAR  all -- 1 $TARGETSCAN_VARIANCE");
  System("hmm-edit PARCLIP.hmm MEAN 3   -- 1 $TARGETSCAN_FG_MEAN");
  }
else {
    System("hmm-edit PARCLIP.hmm DTRK targetscan");
}

print "
# -------------------------------------------
# PREPARING INPUT FILES
# -------------------------------------------
";


print "
# -------------------------------------------
# RUN THE MODEL
# -------------------------------------------
";

System("rm predictions-var*.gff");
my @vars=(0.5,   0.25,  0.2,  0.15,  0.1,  0.01);
my @sens=(0.12,  0.17,  0.2,  0.27,  0.42, 0.62);
my @SNR= (15.7, 12.04,  9.95, 7.07,  5.09, 2.24);
my $N=@vars;
for(my $i=0 ; $i<$N ; ++$i) {
  my $var=$vars[$i];  my $sens=$sens[$i];  my $snr=$SNR[$i];
  print "Running at variance $var\n";
  System("hmm-edit PARCLIP.hmm VAR all 0 $var");
  System("parse $DASH_G 5-45 -p -d PARCLIP.hmm chunk-targetscan > chunk-preds.gff");
  System("identify-miRNAs.pl chunk-preds.gff seeds.fasta > identify.tmp");
  System("combine-miRNA-predictions.pl identify.tmp > chunk-preds.gff");
  System("get-chrom-coords_try.pl chunk-preds.gff > predictions-var$var.gff");
  addScores("predictions-var$var.gff",$sens,$snr);

  if($WANT_BULGE) {
    print "Running bulge models\n";
    my $type="";
    for(my $bulgeType=1 ; $bulgeType<=3 ; ++$bulgeType) {
      $type.="I";
      System("bulge.pl chunk-targetscan seeds.fasta $bulgeType $var");
      System("get-chrom-coords_try.pl bulge-predictions-var$var.gff > chunk-preds.gff");
      System("mv chunk-preds.gff bulge$type-predictions-var$var.gff");
      addScores("bulge$type-predictions-var$var.gff",$sens,$snr);
    }
  }
}

System("cat predictions-var*.gff > $outfile");

print "
# -------------------------------------------
# DONE.  OUTPUT IS IN: $outfile
# -------------------------------------------
";

System("date");


sub addScores
  {
    my ($filename,$sens,$snr)=@_;
    my $tempName=TempFilename::generate();
    open(IN,$filename) || die "can't open file: $filename\n";
    open(OUT,">$tempName") || die "can't write to file: $tempName\n";
    while(<IN>) {
      chomp;
      print OUT "${_}sens=$sens;SNR=$snr;\n";
    }
    close(OUT);
    close(IN);
    System("mv $tempName $filename");
  }


sub System
{
    my ($cmd)=@_;
    if($VERBOSE) { print "$cmd\n" }
    system($cmd);
}

