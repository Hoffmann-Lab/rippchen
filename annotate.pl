#! /usr/bin/env perl
# (c) Konstantin Riege
use strict;
use warnings;
use feature ":5.10";

if ($#ARGV < 2){
	say "usage: annotate.pl [gtf.descr|0] [gtf|0] [deseq.csv|heatmap.ps ..]";
	say "0 or gtf.descr => 4 tab seperated columns: geneID geneName biotype description";
	say "0 or gtf => parsed for gene_id gene_name gene_biotype";
	exit 1;	
}

my %m;
my %mps;
if ($ARGV[0]){
	open I,"<$ARGV[0]" or die $!;
	while(<I>){
		chomp;
		my @l=split/\t/;
		push @l,"" while $#l<3;
		$l[3]=~s/\s+\[Source.+//;
		$m{$l[0]}='"'.$l[0].'","'.$l[1].'","'.$l[2].'","'.$l[3].'"';
		$mps{$l[0]}=$l[1];
	}
	close I;
}

if ($ARGV[1]){
	open I,"<$ARGV[1]" or die $!;
		while(<I>){
			chomp;
			my @l=split/\t/;
			next unless $l[2] eq 'gene';
			my ($g,$n,$b) = ('""','""','""');
			if ($l[-1]=~/gene_id\s+\"([^\"]+)/){
				$g = $1;
				next if exists $mps{$g};
				if ($l[-1]=~/gene_name\s+\"([^\"]+)/){
					$n = '"'.$1.'"';
					$mps{$g} = $n;
				}
				if ($l[-1]=~/gene_biotype\s+\"([^\"]+)/){
					$b = '"'.$1.'"';
				}
				$m{$g} = join",",('"'.$g.'"',$n,$b,'""');
			}
		}
	close I;
}

for (2..$#ARGV){
	my $f=$ARGV[$_];
	my @o=split/\./,$f;
	my $e=$o[-1];
	$o[-1]="annotated";
	push @o,$e;
	my $o=join".",@o;
	open I,"<$f" or die $!;
	open O,">$o" or die $!;
	my $ps;
	while(<I>){
		chomp;
		my @l=split/,/;
		if ($.==1){
			if (/PS-Adobe/){
				$ps=1;
				say O;
			} else {
				$l[0]='"geneID","geneName","biotype","description"';
				say O join",",@l;
			}
			next;
		}
		if ($ps) {
			if (/\((\S+)(\s\+|\s-)*\)\s+(0|0\.25)\s+0\s+t$/) {
				my $n = $mps{$1};
				s/$1/$n/ if $n;
			}
			say O;
		} else {
			$l[0]=~/([^\"]+)/;
			$l[0] = $m{$1};
			$l[0] = '"'.$1.'","","",""' unless $l[0];
			say O join",",@l;
		}

	}
	close I;
	close O;
}
