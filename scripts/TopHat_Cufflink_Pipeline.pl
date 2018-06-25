use lib $ARGV[0];

BEGIN {
	require $ARGV[1];
		$ARGV[1]->import( qw(TGFAM_CONFIG) );
	require $ARGV[2];
		$ARGV[2]->import( qw(TGFAM_PATH) );
	}
use strict;
use Getopt::Long;

my ($stBowtieIndex,$stRead1,$stRead2,$stOut,$nMismatches,$nMinIntron,$nMaxIntron,$nMaxDeletion,$nMaxInsertion,$stCPU,$stGff,$stNoNovelJunction,$stNoNovelIndels,$stMicroExonSearch,$stRefGenome,$nCPU,$stCovSearch);
my $stOptions="";
my $stCuffOptions="";
my $stMergeOptions="";
my @stTopHatOptions=();
my @stCufflinksOptions=();
my @stCuffMergeOptions=();
my @stAssemblyList=();
my $stMergeOut="";
my $stLibrary="";
my $stUsageFlag;
&GetOptions(
	'h|help' => \$stUsageFlag,
	'Bowtie2_Index=s' => \$stBowtieIndex,
	'For=s' => \$stRead1,
	'Rev=s' => \$stRead2,
	'Out=s' => \$stOut,
	'Splice-mismatches=i' => \$nMismatches,
	'min-intron-length=i' => \$nMinIntron,
	'max-intron-length=i' => \$nMaxIntron,
	'max-insertion-length=i' => \$nMaxInsertion,
	'max-deletion-length=i' => \$nMaxDeletion,
	'num-threads=i' => \$nCPU,
	'Gff-Gtf=s' => \$stGff,
	'no-novel-juncs=s' => \$stNoNovelJunction,
	'no-novel-indels=s' => \$stNoNovelIndels,
	'coverage-search=s' => \$stCovSearch,
	'microexon-search=s' => \$stMicroExonSearch,
	'ref-genome=s' => \$stRefGenome,
	'CuffMergeOut=s' => \$stMergeOut,
	'Library=s' => \$stLibrary
);
my $help = qq^
#####################################################################################################################################################################################################
#  Usage: 
#  perl TopHat_Cufflink_Pipeline.pl [Options] -Bowtie2_Index Reference_Index -For Read1,[Read2,...] -Rev Read1,[Read2,...]
#
#  If you want to perform multiple execution in each tissue, try this:
#  
#  Usage of multiple tissues: 
#  perl TopHat_Cufflink_Pipeline.pl [Options] -Bowtie2_Index Reference_Index -For Tis1_Read1,[Tis1_Read2,...]|Tis2_Read1,[Tis2_Read2,...] -Rev Tis1_Read1,[Tis1_Read2,...]|Tis2_Read1,[Tis2_Read2,...]
#
#  Options:
#  -h or -help: help message
#  -Out:			<String>		[ default: tophat_0.out(tophat), Cufflinks_0.out(cufflinks)]
#  	Note: If you perform multiple execution in each tissue, try this:
#  	-0ut OutName1,OutName2,...
#  -CuffMergeOut:		<String>		[ default: ./merge_asm				]
#  -Gff-Gtf:			<Ref.(gff3/gtf)>	[ default: FALSE				]
#  -num-threads:		<int>			[ default: 1					]
#  -min-intron-length:		<int>			[ default: 50(tophat), 50(cufflinks)		]
#  -max-intron-length:		<int>			[ default: 500000(tophat), 300000(cufflinks)	]
#  -max-insertion-length:	<int>			[ default: 3					]
#  -max-deletion-length:	<int>			[ default: 3					]
#  -Splice-mismatches:		<0-2>			[ default: 0(tophat)				]
#  -no-novel-juncs			
#  -no-novel-indels
#  -coverage-search					[ default: FALSE				]
#  -microexon-search					[ default: FALSE				]
#  -ref-genome:			<RefGenome.fa>		[ default: FALSE				]
#####################################################################################################################################################################################################
    ^;
if($stUsageFlag)
{
	die "$help\n";
}
if($stRefGenome ne "")
{
	push(@stCuffMergeOptions,"-s $stRefGenome");
}
if($nMismatches ne "")
{
	push(@stTopHatOptions,"-m $nMismatches");
}
if($nMinIntron ne "")
{
	push(@stTopHatOptions,"-i $nMinIntron");
	push(@stCufflinksOptions,"--min-intron-length $nMinIntron");
}
if($nMaxIntron ne "")
{
	push(@stTopHatOptions,"-I $nMaxIntron");
	push(@stCufflinksOptions,"-I $nMaxIntron");
}
if($nMaxInsertion ne "")
{
	push(@stTopHatOptions,"--max-insertion-length $nMaxInsertion");
}
if($nMaxDeletion ne "")
{
	push(@stTopHatOptions,"--max-deletion-length $nMaxDeletion");
}
if($nCPU ne "")
{
	push(@stTopHatOptions,"-p $nCPU");
	push(@stCufflinksOptions,"-p $nCPU");
	push(@stCuffMergeOptions, "-p $nCPU");
}
if($stGff ne "")
{
	push(@stTopHatOptions,"-G $stGff");
	#push(@stCufflinksOptions,"-g $stGff");
	#push(@stCuffMergeOptions,"-g $stGff");
}
if($stNoNovelJunction ne "")
{
	push(@stTopHatOptions,"--no-novel-juncs");
}
if($stNoNovelIndels ne "")
{
	push(@stTopHatOptions,"--no-novel-indels");
}
if($stCovSearch ne "")
{
	push(@stTopHatOptions,"--coverage-search");
}
if($stMicroExonSearch ne "")
{
	push(@stTopHatOptions,"--microexon-search");
}
if($stLibrary ne "")
{
	push(@stTopHatOptions, "--library-type $stLibrary");
	push(@stCufflinksOptions, "--library-type $stLibrary");
}
for(my $i=0; $i<@stTopHatOptions; $i++)
{
	$stOptions = $stOptions."$stTopHatOptions[$i] ";
}
for(my $k=0; $k<@stCufflinksOptions; $k++)
{
	$stCuffOptions = $stCuffOptions."$stCufflinksOptions[$k] ";
}
for(my $b=0; $b<@stCuffMergeOptions; $b++)
{
	$stMergeOptions = $stMergeOptions."$stCuffMergeOptions[$b] ";
}

my @stRead1 = split(/\|/,$stRead1);
my @stRead2 = split(/\|/,$stRead2);
my @stOut = split(/,/,$stOut);
for(my $j=0; $j<@stRead1; $j++)
{
	if($stRead2[$j] ne "")
	{
		if($stOptions ne "")
		{
			if($stOut[$j] eq "")
			{
				system("$TOPHAT_PATH/tophat -o tophat_$j.out $stOptions"."$stBowtieIndex $stRead1[$j] $stRead2[$j]");
			}
			elsif($stOut[$j] ne "")
			{
				system("$TOPHAT_PATH/tophat -o TopHat_$stOut[$j] $stOptions"."$stBowtieIndex $stRead1[$j] $stRead2[$j]");
			}
		}
		elsif($stOptions eq "")
		{
			if($stOut[$j] eq "")
			{
				system("$TOPHAT_PATH/tophat -o tophat_$j.out $stBowtieIndex $stRead1[$j] $stRead2[$j]");
			}
			elsif($stOut[$j] ne "")
			{
				system("$TOPHAT_PATH/tophat -o TopHat_$stOut[$j] $stBowtieIndex $stRead1[$j] $stRead2[$j]");
			}
		}
	}
	elsif($stRead2[$j] eq "")
	{
		if($stOptions ne "")
		{
			if($stOut[$j] eq "")
			{
				system("$TOPHAT_PATH/tophat -o tophat_$j.out $stOptions"."$stBowtieIndex $stRead1[$j]");
			}
			elsif($stOut[$j] ne "")
			{
				system("$TOPHAT_PATH/tophat -o TopHat_$stOut[$j] $stOptions"."$stBowtieIndex $stRead1[$j]");
			}
		}
		elsif($stOptions eq "")
		{
			if($stOut[$j] eq "")
			{
				system("$TOPHAT_PATH/tophat -o tophat_$j.out $stBowtieIndex $stRead1[$j]");	
			}
			elsif($stOut[$j] ne "")
			{
				system("$TOPHAT_PATH/tophat -o TopHat_$stOut[$j] $stBowtieIndex $stRead1[$j]");	
			}
		}
	}
	if($stCuffOptions ne "")
	{
		if($stOut[$j] eq "")
		{
			system("$CUFFLINKS_PATH/cufflinks -o Cufflinks_$j.out $stCuffOptions"."tophat_$j.out/accepted_hits.bam");
			push(@stAssemblyList, "Cufflinks_$j.out");
		}
		elsif($stOut[$j] ne "")
		{
			system("$CUFFLINKS_PATH/cufflinks $stCuffOptions"."-o Cufflinks_$stOut[$j] TopHat_$stOut[$j]/accepted_hits.bam");
			push(@stAssemblyList, "Cufflinks_$stOut[$j]");
		}
	}
	elsif($stCuffOptions eq "")
	{
		if($stOut[$j] eq "")
		{
			system("$CUFFLINKS_PATH/cufflinks -o Cufflinks_$j.out tophat_$j.out/accepted_hits.bam");
			push(@stAssemblyList, "Cufflinks_$j.out");
		}
		elsif($stOut[$j] ne "")
		{
			system("$CUFFLINKS_PATH/cufflinks -o Cufflinks_$stOut[$j] TopHat_$stOut[$j]/accepted_hits.bam");
			push(@stAssemblyList, "Cufflinks_$stOut[$j]");
		}
	}
}
