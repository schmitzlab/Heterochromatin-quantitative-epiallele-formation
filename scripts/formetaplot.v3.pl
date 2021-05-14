#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
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
my ($gff,$allc,$fOut,$list);
GetOptions(
				"help|?" =>\&USAGE,
				"o:s"=>\$fOut,
				"gff:s"=>\$gff,
				"allc:s"=>\$allc,
                "l:s"=>\$list,
				) or &USAGE;
&USAGE unless ($gff and $allc and $fOut and $list);

sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
Version:	$version
Contact:	zhangyw <yz46606\@uga.edu> 
Program Date:   2016/09/11
Description:	this program is used to generate table for metaplot
Usage:
  Options:
  -gff <file>  gff input file,forced  
  -l <file>  gene gbM,teM,UM classification list,forced  
  -allc <file>  tsv input file,forced  
  -o <dir>  output file,forced  
  -h         Help

USAGE
	print $usage;
	exit;
}
mkdir $fOut if (! -d $fOut);
my %type;
open (IN, $list) or die $!;
while (<IN>) {
    chomp;
    next if ($.==1);
    my @lines=split/\s+/,$_;
    $type{$lines[3]}=$lines[7];
}
close IN;

my $winsize=50;
my $binnum=20;
mkdir $fOut if (! -d $fOut);
my %context=(
"CGA"=>"CG","CGC"=>"CG","CGG"=>"CG","CGT"=>"CG","CG"=>"CG","CGN"=>"CG","CHG"=>"CHG","CHH"=>"CHH",
"CAG"=>"CHG","CCG"=>"CHG","CTG"=>"CHG","CAA"=>"CHH","CAC"=>"CHH","CAT"=>"CHH",
"CCA"=>"CHH","CCC"=>"CHH","CCT"=>"CHH","CTA"=>"CHH","CTC"=>"CHH","CTT"=>"CHH",
);
my @class=("gbM","teM","UM");

&gene_bin_site_methylation();

sub gene_bin_site_methylation{
	open (IN, $gff) or die $!;
	my $up_max     = $winsize * $binnum; 
	my $up_stop    = $binnum;
	my $gb_start   = $binnum+1;
	my $gb_stop    = $binnum*2;
	my $gb_bin     = $binnum;
	my $down_start = $gb_stop+1;
	my $down_stop  = $down_start+$binnum-1;
	my %chr_site;
	my $chr;
	while (<IN>) {
		chomp;
		next if (/\#/||$_=~/^$/);
		my $count1 =0;	## Reset parameter for gene body
		my $count2 =0;	## Reset papameter for downstream
		my @lines=split/\t/,$_;## chr	TAIR	type	start	stop	.	strand	.	genename
		if($lines[0]=~/Chr(\d+)/){
			$chr=$1;
			if ($lines[2]=~/gene/) {
				(my $id)=($lines[-1])=~/ID\=(.*?)\;/;
                if($type{$id}){
                    $lines[2]=$type{$id};
                    for my $i (1..$up_stop) {
                        my $binnum;
                        if ($lines[6] eq "+") {
                            $binnum=$i;
                        }
                        elsif ($lines[6] eq "-") {
                            $binnum=$down_stop + 1 - $i;
                        }
                        my $start=$lines[3]-$up_max+($i-1)*$winsize;
                        my $stop=$start+$winsize-1;
                        for my $j ($start..$stop) {			## print the one bin information
                            $chr_site{"$chr\t$j"}="$lines[2]\t$binnum";	## Gene_group/Bin/
                            #print $chr_site{"$lines[0]]\t$j"};die;
                        }
                    }
				#####################body region
                    my $laststop;
                    for my $i ($gb_start..($gb_stop - 1)){
                        my $binnum;
                        if ($lines[6] eq "+") {
                            $binnum=$i;
                        }
                        elsif ($lines[6] eq "-") {
                            $binnum=$down_stop + 1 - $i;
                        }
                        $count1++;
                        my $range=int(($lines[4]-$lines[3])/$gb_bin);
                        my $start = $lines[3] + $range * ($count1 - 1);
                        my $stop  = $start + $range - 1;
                        $laststop=$stop;
                        for my $j ($start..$stop) {
                            $chr_site{"$chr\t$j"}="$lines[2]\t$binnum";## Gene_group/Bin/
                        }
                    }
                    my $binlast; #####last bin of genebody
                    if ($lines[6] eq "+") {
                            $binlast=$gb_stop;
                        }
                        elsif ($lines[6] eq "-") {
                            $binlast=$down_stop + 1 - $gb_stop;
                    }

                    $laststop++;
                    for my $j ($laststop..$lines[4]) {
                        $chr_site{"$chr\t$j"}="$lines[2]\t$binlast";## Gene_group/Bin/
                    }
				    ##down stream
                    for my $i($down_start..$down_stop) {
                        my $binnum;
                        if ($lines[6] eq "+") {
                            $binnum=$i;
                        }
                        elsif ($lines[6] eq "-") {
                            $binnum=$down_stop + 1 - $i;
                        }
                        $count2++;
                        my $start = $lines[4] + $winsize * ($count2 - 1) + 1;	
                        my $stop  = $start + $winsize - 1; 
                        for my $j ($start..$stop){
                            $chr_site{"$chr\t$j"}="$lines[2]\t$binnum";## Gene_group/Bin/
                        }
                    }
                }
			}
		}
	}
	close IN;
	#print Dumper %chr_site;die;
	my %mCbase;
	my %allbase;
	my $tsvname=basename($allc);
	$tsvname=~s/\.tsv//;
	open (IN, $allc) or die $!;
	while (<IN>) {
		chomp;
		next if ($.==1);
		my @lines=split/\s+/,$_;
		$lines[0]=~s/Chr//;
		#print $_;die;
		if ($chr_site{"$lines[0]\t$lines[1]"}&&$context{$lines[3]}) {
			#print $_;die;
			my ($group,$bin)=split/\s+/,$chr_site{"$lines[0]\t$lines[1]"};
			#print "$group\t$bin";die;
			$mCbase{$group}{$bin}{$context{$lines[3]}}+=$lines[4];
			$allbase{$group}{$bin}{$context{$lines[3]}}+=$lines[5];
		}
	}
	close IN;
	my @text=("CG","CHG","CHH");
	open (OUT, ">$fOut/$tsvname\.bin.table") or die $!;
	print OUT "Group\tbinnum\tCG\tCHG\tCHH\n";
	foreach my $group (@class) {
		foreach my $bin (sort{$a<=>$b}keys %{$allbase{$group}}) {
			print OUT "$group\t$bin";
			my $methylation_level;
			foreach my $con (@text) {
				if ($mCbase{$group}{$bin}{$con}&&$allbase{$group}{$bin}{$con}) {
					$methylation_level = $mCbase{$group}{$bin}{$con}/$allbase{$group}{$bin}{$con};
					print OUT "\t$methylation_level";
				}
				elsif ($allbase{$group}{$bin}{$con}) {
					print OUT "\t0";
				}
				else{
					print OUT "\t0";
				}
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
	#���б��е����ֵ
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
	#���б��е���Сֵ
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
	#��ȡ�ַ������еķ��򻥲����У����ַ�����ʽ���ء�ATTCCC->GGGAAT
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


