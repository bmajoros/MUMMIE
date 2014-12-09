#!/usr/bin/perl
use strict;

die "$0 \"<page title>\"\n" unless @ARGV==1;
my ($title)=@ARGV;


my $html;
$html.="<html><title>$title</title>\n";

my @enums;
my @allClassNames;
my $ls=`ls *.H`;
my @files=split/\s+/,$ls;
foreach my $file (@files)
  {
    my $classes=parseFile($file);
    foreach my $class(@$classes)
      {
	my $name=$class->{name};
	my $base=$class->{base};
	my $members=$class->{members};
	my $template=$class->{template};
	push @allClassNames,$name;

	# NAMED ANCHOR FOR CLASS===============================================
	$html.="<p><br><br>";
	my $anchor="class_$name";
	$html.="<a name=\"$anchor\">";

	# BEGINNING OF CLASS TABLE===========================================
	$html.="<table align=\"center\" width=\"80\%\" border=\"1\" cellspacing=\"0\" cellpadding=\"0\">\n";

	my $header="$name";
	if(defined($base)) {$header.="$base\n"}

	# TEMPLATE<CLASS T>====================================================
	if(defined($template)) 
	  {
	    $template=~s/&/&amp;/g;
	    $template=~s/</&lt;/g;
	    $html.=" <tr><td bgcolor=\"\#000000\" colspan=\"2\" rowspan=\"1\" valign=\"top\"><font color=\"\#FFFFFF\"><center><b>$template</b></center></font></td></tr>\n";
	  }

	# CLASS NAME:===========================================================
	$header=~s/&/&amp;/g;
	$header=~s/</&lt;/g;
	$html.=" <tr><td bgcolor=\"\#cc0000\" colspan=\"2\" rowspan=\"1\" valign=\"top\"><font color=\"\#ffffff\"><center><b>class $header</b></center></font></td></tr>\n";

	# CLASS DESCRIPTION:=================================================
	$html.=" <tr><td bgcolor=\"\#3366ff\" colspan=\"2\" rowspan=\"1\" valign=\"top\"><font color=\"\#FFFFFF\">.</font></td></tr>\n";

	# MEMBERS:===========================================================
	$html.="\n";
	foreach my $member (@$members)
	  {
	    $member=~s/&/&amp;/g;
	    $member=~s/</&lt;/g;

	    $html.=" <tr>\n";

	    # MEMBER NAME:====================================================
	    if($member=~/(.*)(\b\S+)(\(.*)/) 
	      {$member="$1<font color=\"\#3366ff\"><u>$2</u></font>$3"}
	    $member=~s/operator/<font color=\"\#3366ff\"><u>operator<\/u><\/font>/g;
	    $html.="   <td bgcolor=\"\#ccccff\" colspan=\"1\" rowspan=\"1\" valign=\"top\"><font color=\"\#000000\"><b>$member</b></font></td>\n";
	       # face=\"Courier New, Courier, monospace\" 

	    # MEMBER DESCRIPTION:=============================================
	    my $description=".";
	    if($member=~/~.*$name</) {$description="destructor"}
	    elsif($member=~/$name</) {$description="constructor"}
	    $html.="   <td bgcolor=\"\#ffcccc\" width=\"50\%\" colspan=\"1\" rowspan=\"1\" valign=\"top\"><font color=\"\#000000\">$description</font></td>\n";
	    $html.=" </tr>\n";
	  }
	$html.="</table>\n";
      }
  }

# ALPHABETICAL INDEX=====================================================
my $alphaIndex="<p><br><br><table align=\"center\" width=\"80\%\" border=\"1\" cellspacing=\"0\" cellpadding=\"0\">\n";
$alphaIndex.=" <tr><td bgcolor=\"\#3366ff\" colspan=\"2\" rowspan=\"1\" valign=\"top\"><font color=\"\#ffffff\"><center><b>Class Index (Alphabetical) </b></center></font></td></tr>\n";
$alphaIndex.="<tr>\n";
my $i=1;
foreach my $name (@allClassNames)
  {
    $alphaIndex.="<td bgcolor=\"\#ccccff\" colspan=\"1\" rowspan=\"1\" valign=\"top\"><font color=\"\#000000\"><b>$i. <a href=\"#class_$name\">$name</a></b></font></td>\n";
    ++$i;
    if($i%2==1) {$alphaIndex.="</tr><tr>\n"}
  }
$alphaIndex.="</tr>\n";
$alphaIndex.="</table>\n";
$alphaIndex.="<p><br><br></html>\n";

# SUBJECT INDEX========================================================
my $subjIndex="<p><br><br><table align=\"center\" width=\"80\%\" border=\"1\" cellspacing=\"0\" cellpadding=\"0\">\n";
$subjIndex.=" <tr><td bgcolor=\"\#3366ff\" colspan=\"2\" rowspan=\"1\" valign=\"top\"><font color=\"\#ffffff\"><center><b>Classes By Category </b></center></font></td></tr>\n";
$subjIndex.="<tr>\n";
my $n=@allClassNames;
for(my $i=0 ; $i<$n ; ++$i)
  {
    my $name=$allClassNames[$i];
    $subjIndex.="<td bgcolor=\"\#ccccff\" colspan=\"1\" rowspan=\"1\" valign=\"top\"><font color=\"\#000000\"><b><a href=\"#class_$name\">$name</a></b></font></td>\n";
    if(($i+1)%2==0) {$subjIndex.="</tr><tr>\n"}
  }
$subjIndex.="</tr>\n";
$subjIndex.="</table>\n";
$subjIndex.="<p><br><br></html>\n";

# ENUMS=================================================================
my $enums="<p><br><br><table align=\"center\" width=\"80\%\" border=\"1\" cellspacing=\"0\" cellpadding=\"0\">\n";
$enums.="<tr><td bgcolor=\"\#3366ff\" colspan=\"2\" rowspan=\"1\" valign=\"top\"><font color=\"\#ffffff\"><center><b>Enums</b></center></font></td></tr>\n";
foreach my $enum (@enums)
  {
    my $name=$enum->{name};
    $enums.="<tr>\n";
    $enums.="<td bgcolor=\"\#ccccff\" colspan=\"1\" rowspan=\"1\" valign=\"top\"><font color=\"\#000000\"><b>$name</b></font></td>\n";
    my $description;
    my $members=$enum->{members};
    my $numMembers=@$members;
    for(my $i=0 ; $i<$numMembers ; ++$i)
      {
	my $member=$members->[$i];
	$description.="$member";
	if($i+1<$numMembers) {$description.="&nbsp;&nbsp;&nbsp; "}
      }
    $enums.="<td bgcolor=\"\#ffcccc\" colspan=\"1\" rowspan=\"1\" valign=\"top\"><font face=\"Courier New, Courier, monospace\" color=\"\#000000\">$description</font></td>\n";
    $enums.="</tr>\n";
  }
$enums.="</table>\n";

# TOP LEVEL MENU======================================================
my $menu="<center>\n";
$menu.="<br><br><a href=\"#overview\"><big><big><big>Overview</big></big></big></a>\n";
$menu.="<br><br><big><big><big>* <a href=\"#bycat\">Clases by Category</a> *</big></big></big>\n";
$menu.="<br><br><a href=\"#byalpha\"><big><big><big>Alphabetical Index</big></big></big></a>\n";
$menu.="</center>\n";
$menu.="<big><br><br><br><br><br><br><br><br><br><br><br><br><br><br></big><hr>\n";

# OVERVIEW=============================================================
my $overview="<h1>Overview</h1>\n";
$overview.="<center><table align=\"center\" width=\"80\%\" border=\"0\" cellspacing=\"0\" cellpadding=\"0\">\n";
$overview.="<tr><td>No overview available yet.  Sorry.</td></tr>";
$overview.="</table></center>\n";
$overview.="<big><br><br><br><br><br><br><br><br><br><br><br><br><br><br><hr></big>\n";

# ASSEMBLE THE WHOLE PAGE=============================================
my $page="<hr><h1><center>$title</center></h1><hr>";
$page.="$menu";
$page.="<a name=\"overview\"></a>$overview";
$page.="<a name=\"bycat\"></a>$subjIndex";
$page.="<a name=\"byalpha\"></a>$alphaIndex";
$page.="$enums";
$page.="$html";
$page.="<br><hr><center>Contact: <a href=\"mailto:bmajoros\@tigr.org\">bmajoros\@tigr.org</a></center><br>";
print "$page";

#-------------------------------------------------------
sub parseFile
  {
    my ($file)=@_;
    open(IN,$file) || die "can't open $file\n";
    my $templateLine;
    my $classes=[];
    while(<IN>)
      {
	chomp;
	if(/^template\s*</) {$templateLine=$_}
	elsif(/^class\s+(\S+)(\s*:.*)?/)
	  {
	    next if(/;/);
	    my ($className,$base)=($1,$2);
	    my $class={name=>$className,members=>[]};
	    if(defined($base)) {$class->{base}=$base}
	    if(defined($templateLine)) {$class->{template}=$templateLine}
	    parseClass(\*IN,$class);
	    push @$classes,$class;
	    undef $templateLine;
	  }
	elsif(/^enum\s*(\S+)/)
	  {
	    my $name=$1;
	    my $enum={name=>$name,members=>[]};
	    push @enums,$enum;
	    while(<IN>)
	      {
		$_=~s/\/\/.*//g;
		if(/([a-zA-Z_]+)\s*,?/)
		  {
		    my $member=$1;
		    push @{$enum->{members}},$member;
		  }
		elsif(/;/) {last}
	      }
	  }
      }
    close(IN);
    return $classes;
  }
#-------------------------------------------------------
sub parseClass
  {
    my ($file,$class)=@_;
    my $public=0;
    my $members=$class->{members};
    while(<$file>)
      {
	chomp;
	if(/^};/) {last}
	if(/public:/) {$public=1;next}
	if(/private:/) {$public=0;next}
	if(/protected:/) {$public=0;next}
	next unless $public;
	next if(/^#/);# ignore preprocessor directives
	$_=~s/\s*\/\/.*//g;# get rid of comments
	next unless(/\S/);# get rid of blank lines

	# handle multi-line function prototypes:
	while($_!~/[;}]/) {my $next=<$file>;chomp $next;$_.=$next}

	$_=~s/\s*{.*}\s*/;/g;# get rid of inline code

	push @$members,$_;
      }
  }
#-------------------------------------------------------
#-------------------------------------------------------
#-------------------------------------------------------
#-------------------------------------------------------
#-------------------------------------------------------
#-------------------------------------------------------
#-------------------------------------------------------
#-------------------------------------------------------
#-------------------------------------------------------


