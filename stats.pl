#!/usr/bin/perl -w
#Modified 
use strict;
use Getopt::Long;

use vars qw($VERSION);
$VERSION = '0.01';

eval {
    require Pod::Usage;
    Pod::Usage->import();
    1;
} or do {
    *pod2usage = sub {
        die "Error in command line.\n";
    };
};


=head1 OPTIONS

=item  - How to use the program

This program is used to generate the statistics from the output.gff file. 

perl stats.pl output.gff --column 11 > stats.txt

Generates, Statistics as follows:

Stringency  Predictions Genes miRNA

First column represents stringency value (Senstivity or SNR as chossen when running the program)
Second column represents total number of predictions for the chosen stringency
Third column represents list of Genes for chosen stringecny
Fourth column represents total number of miRNA predicted for the specific gene at chosen stringency  

=item B<--column> - set the column to part on

This is the zero-based number of the column.
Multiple columns may be given.

=item B<--separator> - set the column separator

This is the separator for the columns. It defaults
to a tab character ("\t").


=item B<--version> - output version information


=cut

GetOptions(
    'usage'             => \my $usage,
    'column=i'          => \my @col,
    'separator=s'       => \my $sep,
    'verbose'           => \my $verbose,
    'filename-sep=s'    => \my $filename_sep,
    'header-line:s'     => \my $header,
    'help'              => \my $help,
    'version'           => \my $version,
) or pod2usage(2);
pod2usage(1) if $help;
if (defined $version) {
    print "$VERSION\n";
    exit 0;
};
pod2usage("$0: No files given.")  if ((@ARGV == 0) && (-t STDIN));


if (! defined $sep) {
    $sep = "\t";
};

$filename_sep ||= "012";
if ($filename_sep =~ /^\d{3}$/) {
    $filename_sep = chr oct $filename_sep
};

my %lines;
my %hash;
my $number;
if (defined $header) {
    $header ||= 1;
    if ($header =~ /^\d+$/) {
        my $count = $header;
        $header = "";
        $header .= <>
            while $count--;
    };
};

while (<>) {
    s/\r?\n$//;
    s/;/\t/g;
    s/sens=//g;
    s/SNR=//g;
    my @c = split /$sep/o;
    my $key = join $sep, @c[ @col ];
    if (not defined $lines{ $key }) {
        $lines{ $key } ||= [];
    };
    push @{ $lines{$key}}, $_
}

foreach my $ke (sort {$a<=>$b} keys %lines) {

my @len = @{$lines{$ke}};
my $count = scalar(@len);
print "Stringency\tPredictions\tGenes\tmiRNA","\n";
foreach my $line (@len){
my @a=split("\t",$line);
my $i=$a[0];
my $mi="NA"; my $seq="NA"; 
  if ($a[8] =~ m/seq\=(.+)\;miRNA\=(.+)\;/) { $seq=$1; $mi=$2;}
my @mirna = split(";",$mi);
push @{ $hash{$i} }, $mirna[0];
}

my $gene;
my $identifier;

foreach my $k (sort keys %hash) {
   $number = @{$hash{$k}};
    $gene = $k;
    $identifier = $number;
print "$ke\t$count\t$gene\t$identifier","\n";
}
my( $j, $v);
while(( $j, $v) = each %hash ) {
  delete $hash{$j};
}

}

