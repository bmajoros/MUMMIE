#!/usr/bin/perl
use strict;
use SummaryStats;
use ProgramName;

my $name=ProgramName::get();
die "$name <dir>\n" unless @ARGV==1;
my ($dir)=@ARGV;

#process("MUMMIE","mummie");
process("NAIVE","naive");

sub process {
  my ($label,$stem)=@_;
  my (@snr,@ppv,@fdr,@sn,@sp,@f);
  print "$label: ";
  my @files=`ls $dir`;
  foreach my $file (@files) {
    chomp $file;
    if($file=~/$stem.eval$/) {
      open(IN,"$dir/$file") || die $file;
      while(<IN>) {
	if(/SNR=(\S+)\s+PPV=(\S+)\s+FDR=(\S+)/) {
	  push @snr,$1;
	  push @ppv,$2;
	  push @fdr,$3;
	}
	elsif(/Sn=(\S+)\s+Sp=(\S+)\s+F=(\S+)/) {
	  push @sn,$1;
	  push @sp,$2;
	  push @f,$3;
	}
      }
      close(IN);
    }
  }
  report("SNR",\@snr);
  report("PPV",\@ppv);
  report("FDR",\@fdr);
  report("Sn",\@sn);
  report("Sp",\@sp);
  report("F1",\@f);
}



sub report {
  my ($label,$array)=@_;
  my ($mean,$stddev,$min,$max)=SummaryStats::roundedSummaryStats(\@$array);
  print "\t$label: $mean +/- $stddev\n";
}



