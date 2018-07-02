use lib $ARGV[0];

BEGIN {
	require $ARGV[1];
		$ARGV[1]->import( qw(TGFAM_CONFIG) );
	require $ARGV[2];
		$ARGV[2]->import( qw(TGFAM_PATH) );
	}

use strict;

my $stData = $TARGET_GENOME; ### genome fasta file (tempID)
my $stGff3 = $GFF3_OF_TARGET_GENOME;
my $stPrefix = $OUTPUT_PREFIX;
my $stPrefixData = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME";

my $stOut = $stData;
$stOut =~ s/\.fna//;
$stOut =~ s/\.fasta//;
$stOut =~ s/\.fa//;
my $stMerge = "";

my $stData = $TARGET_GENOME; ### genome
my $stID="";
my $newID="";
my $seq="";
my $num=0;

open FH1, ">", "$stPrefixData.tempID.fasta" || die ("Can't open myfile");
open(DATA, "$stData");
open FH, ">", "$stPrefixData.changeID.Info";
while(my $stLine = <DATA>)
{
        chomp($stLine);
        if($stLine =~ /^>([^\s]+)/)
        {
                if($seq ne "")
		{
			$seq=uc($seq);
			print FH1 ">$newID\n$seq\n";
			$seq="";
		}

		$stID=$1;
                $newID="tempCh$num";
		print FH "$stID	$newID\n";
                $num++
        }
	else
	{
		$seq = $seq."$stLine";
	}
}
close DATA;
close FH;

if($seq ne "")
{
	$seq=uc($seq);
	print FH1 ">$newID\n$seq\n";
}

close FH1;

my $stGff = $GFF3_OF_TARGET_GENOME;

my @stGene;
my @stCDS;
my $nIdx1=0;

open FH5, ">", "$RUNNING_PATH/$stPrefixData.sort" || die ("Can't open myfile");

open(DATA, "$stGff");
while(my $stLine = <DATA>)
{
	chomp($stLine);
	$nIdx1++;
	my @stInfo = split /[\s\t]+/, $stLine;
	if ($stInfo[2] eq "gene")
	{
		push(@stGene,$stLine);
	}
	elsif ($stInfo[2] eq "mRNA")
	{
		push(@stGene,$stLine);
	}
	elsif(($stInfo[2] eq "cds")||($stInfo[2] eq "CDS"))
	{
		push(@stCDS, $stLine);
	}
	if($stInfo[1] eq "")
	{
		next if($nIdx1 == 1);
		@stCDS = sort
		{
			($a =~ /^[^\s]+	[^\s]+	[^\s]+	([0-9]+)/)[0]<=>($b =~ /^[^\s]+	[^\s]+	[^\s]+	([0-9]+)/)[0]
		}@stCDS;
		for(my $i=0; $i<@stGene; $i++)
		{
			print FH5 "$stGene[$i]\n";
		}
		for(my $i=0; $i<@stCDS; $i++)
		{
			my @stTemp = split /[\s\t]+/, $stCDS[$i];
			$stTemp[2] = "exon";	
			my $stExon = join ("	", @stTemp);
			print FH5 "$stExon\n";
			print FH5 "$stCDS[$i]\n";
		}
		print FH5 "$stLine\n";
		@stCDS = ();
		@stGene = ();
	}
}
close(DATA);
if($#stCDS != -1)
{
		@stCDS = sort
		{
			($a =~ /^[^\s]+	[^\s]+	[^\s]+	([0-9]+)/)[0]<=>($b =~ /^[^\s]+	[^\s]+	[^\s]+	([0-9]+)/)[0]
		}@stCDS;
		for(my $i=0; $i<@stGene; $i++)
		{
			print FH5 "$stGene[$i]\n";
		}
		for(my $i=0; $i<@stCDS; $i++)
		{
			my @stTemp = split /[\s\t]+/, $stCDS[$i];
			$stTemp[2] = "exon";	
			my $stExon = join ("	", @stTemp);
			print FH5 "$stExon\n";
			print FH5 "$stCDS[$i]\n";
		}
		print FH5 "\n";
		@stCDS = ();
		@stGene = ();
}

close FH5;

my $stData = "$RUNNING_PATH/$stPrefixData.sort"; ### gff
my $stInfo = "$stPrefixData.changeID.Info"; ### changeID Info
my %ID;

open FH2, ">", "$stPrefixData.tempID.gff3" || die ("Can't open myfile");

open(DATA, "$stInfo");
while(my $stLine = <DATA>)
{
        chomp($stLine);
        my @stList = split /[\t]/, $stLine;
        $ID{$stList[0]}=$stList[1];
}
close DATA;

open(DATA, "$stData");
while(my $stLine = <DATA>)
{
        chomp($stLine);
        if($stLine =~ /^[^\s]+/)
        {
		my @stList = split /[\t]/, $stLine;
                my $newID = $ID{$stList[0]};
                
                print FH2 "$newID	$stList[1]	$stList[2]	$stList[3]	$stList[4]	$stList[5]	$stList[6]	$stList[7]	$stList[8]\n";
        }
        else
        {
                print FH2 "$stLine\n";
        }
}
close DATA;
close FH2;


###     Six frame translation   ###

my $stData= "$TGFAM_SCRIPTS_PATH/CodonUsage"; #codon usage
my $stGenome= "$stPrefixData.tempID.fasta"; #genome fasta file, sequences should be connected 
my ($count, $count1, $countX, $num)=(0,0,0,0);
my $length="";
my $last="";
my $str="";
my $end="";
my $chr="";
my $amino="";
my %codon;
my $totalseq="";


open FH3, ">", "$OUTPUT_PREFIX.SixFrameTranslation" || die ("Can't open myfile");

open (DATA, "$stData");
while(my $stLine = <DATA>)
{
	chomp($stLine);
	my @stList = split /[\t]/, $stLine;
	$codon{$stList[0]}="$stList[1]";	
}
close DATA;

open (DATA, "$stGenome");
while(my $stLine = <DATA>)
{
	chomp($stLine);
	
	if($stLine =~ /^>/)
	{
		$stLine =~ />([^\s]+)/;
		$chr = $1;
	}
	else
	{
		next if (length($stLine)<100); # Seungill added to prevent abnormal quit
		for(my $i=0; $i<3; $i++)
		{
			my $printCheck=0;
			$length = length($stLine) - $i;
                        if($length % 3 == 0)
                        {
                                $last = length($stLine);
                        }
                        elsif($length % 3 == 1)
                        {
                                $last = length($stLine) - 1;
                        }
                        elsif($length % 3 == 2)
                        {
                                $last = length($stLine) - 2;
                        }

			for(my $j=$i; $j<$last; $j=$j+3)
			{
				$printCheck=0;
				$amino = $codon{substr($stLine, $j, 3)};
				
				if(($amino ne "*")&&($amino ne ""))
				{
					$count1++;
					$countX=0;
					if($count1 == 1)
					{
						$num= $num+1;
						$str= $j+1;
						print FH3 "$chr	$num	$amino"
					}
					else
					{
						print FH3 "$amino";
					}
				}
				elsif($amino eq "")
				{
					$countX++;
					$count1++;
					if($count1 == 1)
					{
						$str= $j+1;
						$num= $num+1;
						print FH3 "$chr	$num	X";
					}
					elsif($countX == 100)
					{	
						$end= $j+3; #####
						print FH3 "X	$str	$end	+\n";
						$printCheck++;
						$countX=0;
						$count1=0;
					}
					else
					{
						print FH3 "X";
					}
				}
				else
				{
					$end= $j+3; ##### 
					$countX=0;
					if($count1 > 0)
					{
						print FH3 "*	$str	$end	+\n";
						$printCheck++;
					}
					$count1=0;
				}
			}

			if(($amino ne "*")&&($printCheck == 0))
			{
#				my $last= length($stLine);
				print FH3 "	$str	$last	+\n";
			}
			$count1=0;
			$countX=0;
		}
	
		my $revcom= reverse($stLine); 
		$revcom =~ tr/ATGCatgc/TACGtacg/;
		$totalseq= length($stLine);

		for(my $k=0; $k<3; $k++)
		{
			my $printCheck=0;

			$length = length($stLine) - $k;
                        if($length % 3 == 0)
                        {
                                $last = length($stLine);
                        }
                        elsif($length % 3 == 1)
                        {
                                $last = length($stLine) - 1;
                        }
                        elsif($length % 3 == 2)
                        {
                                $last = length($stLine) - 2;
                        }

			for(my $h=$k; $h<$last; $h=$h+3)
			{
				$printCheck=0;
				$amino = $codon{substr($revcom, $h, 3)};
				if(($amino ne "*")&&($amino ne ""))
				{
					$count1++;
					$countX=0;
					if($count1 == 1)
					{
						$num= $num+1;
						$str= $h+1;
						print FH3 "$chr	$num	$amino";
					}
					else
					{	
						print FH3 "$amino";
					}
				}
				elsif($amino eq "")
				{
					$countX++;
					$count1++;
					if($count1 == 1)
					{
						$str= $h+1;
						$num= $num+1;
						print FH3 "$chr	$num	X";
					}
					elsif($countX == 100)
					{		
						$end= $h+3; #####
						my $mlength= $end-$str;
						$str= $totalseq+1-$end;
						$end= $str+$mlength;
						print FH3 "X	$str	$end	-\n";
						$printCheck++;
						$countX=0;
						$count1=0;
					}
					else
					{
						print FH3 "X";
					}
				} 
				else
				{
					$end= $h+3; #####
					$countX=0;
					if($count1 > 0)
					{
						my $mlength= $end-$str;
						$str= $totalseq+1-$end;
                                                $end= $str+$mlength;
						print FH3 "*	$str	$end	-\n";
						$printCheck++;
					}
					$count1=0;
				}
			}
			
			if(($amino ne "*")&&($printCheck == 0))
                        {
#              			my $last= length($stLine);
				my $mlength= $last-$str;
				$str= $totalseq+1-$last;
                                $end= $str+$mlength;
                                print FH3 "	$str	$end	-\n";
			}
			$count1=0;
			$countX=0;
		}
	}
}
close DATA;
close FH3;


my $stData = "$OUTPUT_PREFIX.SixFrameTranslation";

open FH4, ">", "$OUTPUT_PREFIX.SixFrameTranslation.PEP.fasta" || die ("Can't open myfile");

open(DATA, "$stData");
while(my $stLine = <DATA>)
{
	chomp($stLine);
	my @stList = split /[\t]/, $stLine;

	print FH4 ">$stList[0]_$stList[1]_$stList[3]_$stList[4]_$stList[5]\n$stList[2]\n";
}
close DATA;
close FH4;

system("$BLAST_BIN_PATH/makeblastdb -in $RESOURCE_PROTEIN -dbtype prot -out $BLAST_DB_NAME");
