use strict;
use warnings;

# post-processing for bedtools intersect applied to
#   -a BED feature coverage (variable column number)
#   -b BED6-formatted base coverage runs

# variables
my ($nMapUniq, $nBedCol) = @ARGV; # nBedCol allows generic BED call
my $ftMaxI  = $nBedCol - 1;
my $runCovI = $ftMaxI + 5; # 13 for transcriptome
my $nBasesI = $ftMaxI + 7; # 15 for transcriptome
my $mrc = $nMapUniq / 1e6;
my $prec = 1e5;
my ($pfc, $pft, $pln) = (0);

# thread rows; sequential rows encompass all coverage runs in a feature
while(my $line = <STDIN>){
    
    # strip out just the gene feature bed
    chomp $line;
    my @f  = split("\t", $line);
    my $ft = join("\t", @f[0..$ftMaxI]);

    # commit a run of bases in a given feature
    if ($pft){
        if($pft ne $ft){
            printFeature();
            $pfc = 0;
            $pft = $ft;
            $pln = $f[2] - $f[1];    
        }
    } else { # initialize the first feature
        $pft = $ft;
        $pln = $f[2] - $f[1];        
    }
    
    # sum the fractional coverages for all covered bases within the feature    
    $pfc += $f[$runCovI] * $f[$nBasesI]; # multiple per-base coverage times number of overlapping bases
}

# finish the last feature
printFeature();

# add length, count, density and RPKM to gene feature bed
sub printFeature {
    $pln and print join("\t",
        $pft,
        $pln,
        int($pfc                 * $prec + 0.5) / $prec,        
        int($pfc/$pln            * $prec + 0.5) / $prec,        
        int($pfc/($pln/1e3)/$mrc * $prec + 0.5) / $prec
    ), "\n";
}

#=========================================================================
# FEATURE FILE FORMAT
# extension of BED, called xxx.bed.bgz
# compressed by bgzip, indexed by tabix
#-------------------------------------------------------------------------
#  1  chrom
#  2  start (0-based)
#  3  end (1-based, i.e. half-open)
#  4  name (gene identifier or na)
#  5  score (per BED, 0 to 1000, only used for display purposes)
#  6  strand (+ or -)
#-------------------------------------------------------------------------
#  7  gene length (span of gene, including all exons) OR number of genes crossed
#-------------------------------------------------------------------------
#  8  feature type (e.g. exon, intron_antisense, gene)
#  9  splice index in gene (left-to-right order, 1-based, 0=na)
#-------------------------------------------------------------------------
#  10 feature length ((summed) length of input feature(s))
#  11 count (fractional number of reads in (merged) feature(s))
#  12 density = count/length
#  13 RPKM = count/(length/1E3)/(mapped read count/1E6)
#=========================================================================
