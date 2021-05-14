#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
#use Math::CDF qw(:all);
#use Statistics::Multtest qw(:all);
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $Time_Start = &sub_format_datetime(localtime($BEGIN_TIME));
print "Program Starts Time:$Time_Start\n";
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($gff,$fOut,$tsv);
GetOptions(
				"help|?" =>\&USAGE,
				"o:s"=>\$fOut,

				"gff:s"=>\$gff,
				"tsv:s"=>\$tsv,
				) or &USAGE;
&USAGE unless ($gff and $fOut and $tsv);

sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
Version:	$version
Contact:	zhangyw <yz46606\@uga.edu> 
Program Date:   2016/09/11
Description:	this program is used to find the methylation level of each gene
  Options:
  -gff <file>  input file,gff,forced 
  -tsv <file>  input dir,tsv,forced 
  -o <dir>  output file,forced  
  -h         Help

USAGE
	print $usage;
	exit;
}
mkdir $fOut if (! -d $fOut);
my %context=(
"CGA"=>"CG","CGC"=>"CG","CGG"=>"CG","CGT"=>"CG","CG"=>"CG","CGN"=>"CG","CHG"=>"CHG","CHH"=>"CHH",
"CAG"=>"CHG","CCG"=>"CHG","CTG"=>"CHG","CAA"=>"CHH","CAC"=>"CHH","CAT"=>"CHH",
"CCA"=>"CHH","CCC"=>"CHH","CCT"=>"CHH","CTA"=>"CHH","CTC"=>"CHH","CTT"=>"CHH",
);
my @code=("CG","CHG","CHH");

######################################1. reading gff and protein_primaryTranscriptOnly.fa files to get the cds location of primary transcripts

my %class;###rearrange the cds order based on strand state
my %geneloc;###recording cds in original order of gff files
my %ordergene;
$/="\n";

open (IN, $gff) or die $!;
my $chr;
while (<IN>) {
	chomp;
	next if (/^$/||/^\#/);
	my @lines=split/\s+/,$_;
	if ($lines[0]=~/Chr(\d+)/) {
		$chr=$1;
		if ($lines[2]=~/gene/||$lines[2]=~/transposable_element_geneq/) {
			my $gene;
			if ($lines[8]=~/ID\=(.*?)\;/) {
				$gene=$1;
				$geneloc{$chr}{$gene}="$lines[3] $lines[4]";
				$class{$gene}=$lines[2];
				push @{$ordergene{$chr}},$gene;
			}
		}
	}
}
close IN;


my $name=basename($tsv);
my %record=();
my %pos=();
open (IN, $tsv) or die $!;
while (<IN>) {
	chomp;
	next if (/^$/||/^\#/||$.==1||/total/);
	my @lines=split/\s+/,$_;
	if ($lines[5]>=3) {
		$record{$lines[0]}{$lines[1]}=$_;
		push @{$pos{$lines[0]}},$lines[1];
	}
}
close IN;

my %Mstat=();my %Tstat=();
my %MAstat=();my %TAstat=();
my %Mreads=();my %Treads=();
my %MAreads=();my %TAreads=();
foreach my $chr (sort {$a<=> $b}keys %geneloc) {
	my $shift=0;
	MM:foreach my $gene (@{$ordergene{$chr}}) {
		my ($start,$end)=split/\s+/,$geneloc{$chr}{$gene};
		if ($pos{$chr}&&$shift<@{$pos{$chr}}) {
			for (my $i=$shift;$i<@{$pos{$chr}};$i++) {
				if ($pos{$chr}[$i]>=$start&&$pos{$chr}[$i]<=$end) {
					my @info=split/\s+/,$record{$chr}{$pos{$chr}[$i]};
					if ($context{$info[3]}) {
						$Mreads{$gene}{$context{$info[3]}}+=$info[4];
						$Treads{$gene}{$context{$info[3]}}+=$info[5];
						$shift=$i+1;
					}
				}
				elsif ($pos{$chr}[$i]>$end) {
					$shift=$i+1;
					next MM;
				}
				elsif ($pos{$chr}[$i]<$start) {
						$shift=$i+1;
				}
			}
		}
	}
}
open (OUT, ">$fOut/$name.methyl.out") or die $!;

foreach my $chr (sort {$a<=> $b}keys %geneloc) {
	foreach my $gene (@{$ordergene{$chr}}) {
		if ($Treads{$gene}) {
			print OUT "$gene\t$class{$gene}\t$chr\t$geneloc{$chr}{$gene}";
			my $ratio;
			foreach my $cont (@code) {
				if (!$Mreads{$gene}{$cont}) {
					$Mreads{$gene}{$cont}="0";
					$ratio="0";
				}
				if (!$Treads{$gene}{$cont}) {
					$Treads{$gene}{$cont}="0";
					$ratio="0";
				}
				else{
					$ratio=$Mreads{$gene}{$cont}/$Treads{$gene}{$cont};
				}
				print OUT "\t$cont\t$Mreads{$gene}{$cont}\t$Treads{$gene}{$cont}\t$ratio";
			}
			print OUT "\n";
		}
	}
}

#######################################################################################
my $Time_End   = sub_format_datetime(localtime(time()));
print STDOUT "Program Ends Time:$Time_End\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
#######################################################################################

sub ABSOLUTE_DIR{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir\n";
		exit;
	}
	chdir $cur_dir;
	return $return;
}

#######################################################################################

sub max{#&max(lists or arry);
	#求列表中的最大值
	my $max=shift;
	my $temp;
	while (@_) {
		$temp=shift;
		$max=$max>$temp?$max:$temp;
	}
	return $max;
}

#######################################################################################

sub min{#&min(lists or arry);
	#求列表中的最小值
	my $min=shift;
	my $temp;
	while (@_) {
		$temp=shift;
		$min=$min<$temp?$min:$temp;
	}
	return $min;
}

#######################################################################################

sub revcom(){#&revcom($ref_seq);
	#获取字符串序列的反向互补序列，以字符串形式返回。ATTCCC->GGGAAT
	my $seq=shift;
	$seq=~tr/ATCGatcg/TAGCtagc/;
	$seq=reverse $seq;
	return uc $seq;			  
}

#######################################################################################

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

#######################################################################################

sub sub_format_datetime {#Time calculation subroutine
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


