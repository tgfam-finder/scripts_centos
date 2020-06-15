use lib $ARGV[0];

BEGIN {
	require $ARGV[1];
		$ARGV[1]->import( qw(TGFAM_CONFIG) );
	require $ARGV[2];
		$ARGV[2]->import( qw(TGFAM_PATH) );
	}

use strict;

print "########## 5.Generating_FinalGeneModel.pl is started ###########\n";

my $stAssembly = "$CDS_OF_TARGET_GENOME";
my $stCodon = "$TGFAM_SCRIPTS_PATH/CodonUsage";
my $stMergeOut = "$RUNNING_PATH/$OUTPUT_PREFIX.PEP.Consensus";

my $stName = "";
my $stSeq = "";
my %stCodon = {};
my $nCut = 0;
my $stTemp = "tempaaaaa";

open(OUT, ">$stTemp");
print OUT ">AA";
close(OUT);

system ("cat $stAssembly $stTemp > NewAssembly.fa");

$stAssembly = "NewAssembly.fa";

open(DATA, "$stCodon");
while(my $stLine = <DATA>)
{
	chomp($stLine);
	my @stInfo = split /[\s\t]+/, $stLine;
	$stCodon{$stInfo[0]} = $stInfo[1];
}
close(DATA);

open(DATA, "$stAssembly");
open(OUT, ">$stMergeOut");

while(my $stLine = <DATA>)
{
	chomp($stLine);
	if($stLine =~ /^>/)
	{
		my @stTrans;
		my @stSubSeq;
		my $stTrans1;
		if($stName ne "")
		{
			for(my $i=0; $i<length($stSeq); $i= $i+3)
			{
				if($i+3<=length($stSeq))
				{
					$stSubSeq[0] = substr($stSeq,$i,3);
					$stSubSeq[1] = substr($stSeq,$i+1,3);
					$stSubSeq[2] = substr($stSeq,$i+2,3);
					if($stCodon{$stSubSeq[0]} eq "")
					{
						$stCodon{$stSubSeq[0]} = "X";
					}
					if($stCodon{$stSubSeq[1]} eq "")
					{
						$stCodon{$stSubSeq[1]} = "X";
					}
					if($stCodon{$stSubSeq[2]} eq "")
					{
						$stCodon{$stSubSeq[2]} = "X";
					}

					$stTrans[0] = $stTrans[0].$stCodon{$stSubSeq[0]};
					$stTrans[1] = $stTrans[1].$stCodon{$stSubSeq[1]};
					$stTrans[2] = $stTrans[2].$stCodon{$stSubSeq[2]};
				}
			}
			for(my $i=0; $i<3; $i++)
			{
				if ($stTrans[$i] !~ /\*./)
				{
					$stTrans1 = $stTrans[$i];
					last;
				}
			}
			if((length($stSeq) % 3 == 0 )&&($stTrans[0] !~ /\*./))
			{
				$stTrans1 = $stTrans[0];
			}

			if($stTrans1 !~ /\*./)
			{
				if(($stTrans1 =~ /^M/) && ($stTrans1 =~ /\*$/))
				{
					$stName = $stName."	FULL";
				}
				else
				{
					$stName = $stName."	Partial";
				}
				print OUT ">$stName\n$stTrans1\n";

			}
			else
			{
				print ">$stName\n$stTrans1\n";
			}
		}
		$stName = $stLine;
		$stName =~ s/>//g;
		$stSeq = "";
	}
	else
	{
		$stLine =~ tr/acgt/ACGT/;
		$stSeq = $stSeq.$stLine;
	}
}
close(DATA);
close(OUT);

my $stPfamID=$TARGET_DOMAIN_ID; #
my $stPubGff="$RUNNING_PATH/$OUTPUT_PREFIX.tempID.gff3"; #
my $stPubPEP="$RUNNING_PATH/$OUTPUT_PREFIX.PEP.Consensus";
my $stPubTSV="$RUNNING_PATH/$OUTPUT_PREFIX.tsv.addition";
my $stExceptID="$EXCLUDED_DOMAIN_ID"; #

my $stISGAP_Gff=glob("*.Final.GeneModel.Whole.Repre.sort.gff3.filter");
my $stISGAP_PEP=glob("*.Final.GeneModel.Whole.Repre.sort.PEP.fa.filter");
my $stISGAP_TSV=glob("*.Final.GeneModel.Whole.Repre.sort.PEP.fa.noStop.tsv");
my $stIS_AUG_Gff=glob("*.ISGAP.Augustus.Whole.Repre.sort.gff3");
my $stIS_AUG_PEP=glob("*.ISGAP.Augustus.Whole.Repre.sort.PEP.fa");
my $stIS_AUG_TSV=glob("*.ISGAP.Augustus.Whole.Repre.sort.PEP.fa.noStop.tsv");
my $stPM_Gff=glob("*.RefPEP.Consensus.Repre.gff3");
my $stPM_PEP=glob("*.RefPEP.Consensus.Repre.PEP.fa");
my $stPM_TSV=glob("*.RefPEP.Consensus.Repre.PEP.fa.noStop.tsv");
my $stAug_Gff=glob("*.augustus.gff3.filter");
my $stAug_PEP=glob("*.augustus.PEP.fa.filter");
my $stAug_TSV=glob("*.augustus.PEP.fa.filter.noStop.tsv");
my $stPM_Aug_Gff=glob("*.gff3.PM.Augustus");
my $stPM_Aug_PEP=glob("*.PEP.fa.PM.Augustus");
my $stIS_Aug_Par_Gff=glob("*.gff3.ISGAP.Aug.Partial_Augustus");
my $stIS_Aug_Par_PEP=glob("*.PEP.fa.ISGAP.Aug.Partial_Augustus");

my %stGff;

my %stMerge1 = &Merge($stPubGff, $stPubPEP, $stPubTSV, $stPfamID, "Full", \%stGff, $stExceptID, "PUBLIC");
my %stMerge2 = &Merge($stISGAP_Gff, $stISGAP_PEP, $stISGAP_TSV, $stPfamID, "Full", \%stMerge1, $stExceptID, "TGFam.ISGAP");
my %stMerge3 = &Merge($stIS_AUG_Gff, $stIS_AUG_PEP, $stIS_AUG_TSV, $stPfamID, "Full", \%stMerge2, $stExceptID, "TGFam.ISGAP");
my %stMerge4 = &Merge($stIS_Aug_Par_Gff, $stIS_Aug_Par_PEP, $stAug_TSV, $stPfamID, "Full", \%stMerge3, $stExceptID, "TGFam.ISGAP");
my %stMerge5 = &Merge($stPM_Aug_Gff, $stPM_Aug_PEP, $stAug_TSV, $stPfamID, "Full", \%stMerge4, $stExceptID, "TGFam.PM");
my %stMerge6 = &Merge($stPM_Gff, $stPM_PEP, $stPM_TSV, $stPfamID, "Full", \%stMerge5, $stExceptID, "TGFam.PM");
my %stMerge7 = &Merge($stPubGff, $stPubPEP, $stPubTSV, $stPfamID, "Partial", \%stMerge6, $stExceptID, "PUBLIC");
my %stMerge8 = &Merge($stISGAP_Gff, $stISGAP_PEP, $stISGAP_TSV, $stPfamID, "Partial", \%stMerge7, $stExceptID, "TGFam.ISGAP");
my %stMerge9 = &Merge($stIS_AUG_Gff, $stIS_AUG_PEP, $stIS_AUG_TSV, $stPfamID, "Partial", \%stMerge8, $stExceptID, "TGFam.ISGAP");
my %stMerge10 = &Merge($stIS_Aug_Par_Gff, $stIS_Aug_Par_PEP, $stAug_TSV, $stPfamID, "Partial", \%stMerge9, $stExceptID, "TGFam.ISGAP");
my %stMerge11 = &Merge($stAug_Gff, $stAug_PEP, $stAug_TSV, $stPfamID, "Full", \%stMerge10, $stExceptID, "TGFam.AUG");
my %stMerge12 = &Merge($stPM_Aug_Gff, $stPM_Aug_PEP, $stAug_TSV, $stPfamID, "Partial", \%stMerge11, $stExceptID, "TGFam.PM");
my %stMerge13 = &Merge($stPM_Gff, $stPM_PEP, $stPM_TSV, $stPfamID, "Partial", \%stMerge12, $stExceptID, "TGFam.PM");
my %stMerge14 = &Merge($stAug_Gff, $stAug_PEP, $stAug_TSV, $stPfamID, "Partial", \%stMerge13, $stExceptID, "TGFam.AUG");
#my %stMerge10 = &Merge($stAug_Gff, $stAug_PEP, $stAug_TSV, $stPfamID, "Partial", \%stMerge9, $stExceptID);


open(FH, ">$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.Merge.gff3");
foreach my $stKey (sort{($a =~ /tempCh([0-9]+)/)[0] <=> ($b =~ /tempCh([0-9]+)/)[0] || ($a =~ /[^\s]+_([^\s]+)_[^\s]+/)[0] <=> ($b =~ /[^\s]+_([^\s]+)_[^\s]+/)[0]}keys %stMerge14)
{
	my @stInfo=split(/\n/,$stMerge13{$stKey});
	print FH "$stInfo[0]\n" if ($stInfo[0] =~ /[a-zA-Z]/);
	my $stmRNA=$stInfo[0];
	$stmRNA =~ s/	gene	/	mRNA	/;
	print FH "$stmRNA\n" if ($stInfo[0] =~ /[a-zA-Z]/);
	my @stInfo1=split(/\t/,$stInfo[0]);
	my $stAttri=$stInfo1[8];
	for(my $i=1; $i<@stInfo; $i++)
	{
		my @stInfo2=split(/\t/,$stInfo[$i]);
		if($stInfo[$i] =~ /	CDS	[0-9]+/)
		{
			print FH "$stInfo2[0]	$stInfo2[1]	exon	$stInfo2[3]	$stInfo2[4]	$stInfo2[5]	$stInfo2[6]	$stInfo2[7]	$stAttri\n";
			print FH "$stInfo2[0]	$stInfo2[1]	CDS	$stInfo2[3]	$stInfo2[4]	$stInfo2[5]	$stInfo2[6]	$stInfo2[7]	$stAttri\n";
		}
		elsif($stInfo[$i] =~ /	cds	[0-9]+/)
		{
			print FH "$stInfo2[0]	$stInfo2[1]	exon	$stInfo2[3]	$stInfo2[4]	$stInfo2[5]	$stInfo2[6]	$stInfo2[7]	$stAttri\n";
			print FH "$stInfo2[0]	$stInfo2[1]	CDS	$stInfo2[3]	$stInfo2[4]	$stInfo2[5]	$stInfo2[6]	$stInfo2[7]	$stAttri\n";
		}
	}
	print FH "\n"  if ($stInfo[0] =~ /[a-zA-Z]/); 
#	print "$stMerge10{$stKey}\n";
}

sub Merge
{
	my $stData1=$_[0];
	my $stData2=$_[1];
	my $stData3=$_[2];
	my @stPfamInfo=split(/,/,$_[3]);
	my $stType=$_[4];
	my %stFinalGff=%{$_[5]};
	my @stExceptInfo=split(/,/,$_[6]);
	my $stCol1_Name=$_[7];

	my %stExistPfam;
	my %stExistPEP;
	my @stGffInfo = ();
	my ($stID,$stSeq) = ("", "");
	my $nIdx1=0;

	open(DATA, "$stData3");
	
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		my @stInfo=split(/\t/,$stLine);
		next if($stExistPfam{$stInfo[0]} eq "No");
		for(my $j=0; $j<@stExceptInfo; $j++)
		{
			if($stInfo[4] eq "$stExceptInfo[$j]")
			{
				$stExistPfam{$stInfo[0]} = "No";
				last;
			}
		}
		next if($stExistPfam{$stInfo[0]} eq "No");
		next if($stExistPfam{$stInfo[0]} eq "Exist");
		for(my $i=0; $i<@stPfamInfo; $i++)
		{
			if($stInfo[4] eq "$stPfamInfo[$i]")
			{
				$stExistPfam{$stInfo[0]} = "Exist";
				last;
			}
		}
	}
	close DATA;
	
	open(DATA, "$stData2");
	
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		if($stLine =~ /^>([^\s]+)/)
		{
			if($stSeq ne "")
			{
				if((($stSeq =~ /^M/)||($stSeq =~ /\*$/))&&($stExistPfam{$stID} eq "Exist"))
				{
					if(($stType eq "Full")&&($stSeq =~ /^M.+\*$/))
					{
						$stExistPEP{$stID} = "OK";
						#print "$stData1	Full	$stID\n";
					}
					elsif(($stType eq "Partial")&&($stSeq !~ /^M.+\*$/))
					{
						$stExistPEP{$stID} = "OK";
						#print "$stData1	Partial	$stID\n";
					}
				}
				$stSeq = "";
			}
			$stID = $1;
		}
		else
		{
			$stSeq = $stSeq."$stLine";
		}
	}
	close DATA;
	if($stSeq ne "")
	{
		if((($stSeq =~ /^M/)||($stSeq =~ /\*$/))&&($stExistPfam{$stID} eq "Exist"))
		{
			if(($stType eq "Full")&&($stSeq =~ /^M.+\*$/))
			{
				$stExistPEP{$stID} = "OK";
						#print "$stData1	Full	$stID\n";
			}
			elsif(($stType eq "Partial")&&($stSeq !~ /^M.+\*$/))
			{
				$stExistPEP{$stID} = "OK";
						#print "$stData1	Partial	$stID\n";
			}
		}
		$stSeq = "";
	}
	
	open(DATA, "$stData1");

	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		if($stLine ne "")
		{
			push(@stGffInfo, "$stLine");
		}
		else
		{
			my @stInfo = split(/\t/,$stGffInfo[0]);
			$stInfo[8] =~ s/ID=//;
			foreach my $stKey (keys %stFinalGff)
			{
				my @stCheck=split(/_/,$stKey);
				next if(($stInfo[3] > $stCheck[2])||($stInfo[4] < $stCheck[1]));
				if(($stInfo[0] eq "$stCheck[0]")&&($stInfo[3] >= $stCheck[1])&&($stInfo[3] <= $stCheck[2]))
				{
					$nIdx1 = 1;
					last;
				}
				elsif(($stInfo[0] eq "$stCheck[0]")&&($stInfo[4] >= $stCheck[1])&&($stInfo[4] <= $stCheck[2]))
				{
					$nIdx1 = 1;
					last;
				}
				elsif(($stInfo[0] eq "$stCheck[0]")&&($stInfo[3] <= $stCheck[1])&&($stInfo[4] >= $stCheck[2]))
				{
					$nIdx1 = 1;
					last;
				}
			}
			if(($nIdx1 == 0)&&($stExistPEP{$stInfo[8]} eq "OK"))
			{
				for(my $i=0; $i<@stGffInfo; $i++)
				{
					my @stColInfo=split(/\t/,$stGffInfo[$i]);
					$stGffInfo[$i] =~ s/	$stColInfo[1]	/	$stCol1_Name	/;
					$stFinalGff{$stInfo[0]."_$stInfo[3]_$stInfo[4]"} = $stFinalGff{$stInfo[0]."_$stInfo[3]_$stInfo[4]"}."$stGffInfo[$i]\n";
				}
			}
			$nIdx1 = 0;
			@stGffInfo = ();
		}
	}
	if($#stGffInfo != -1)
	{
			my @stInfo = split(/\t/,$stGffInfo[0]);
			$stInfo[8] =~ s/ID=//;
			foreach my $stKey (keys %stFinalGff)
			{
				my @stCheck=split(/_/,$stKey);
				next if(($stInfo[3] > $stCheck[2])||($stInfo[4] < $stCheck[1]));
				if(($stInfo[0] eq "$stCheck[0]")&&($stInfo[3] >= $stCheck[1])&&($stInfo[3] <= $stCheck[2]))
				{
					$nIdx1 = 1;
					last;
				}
				elsif(($stInfo[0] eq "$stCheck[0]")&&($stInfo[4] >= $stCheck[1])&&($stInfo[4] <= $stCheck[2]))
				{
					$nIdx1 = 1;
					last;
				}
				elsif(($stInfo[0] eq "$stCheck[0]")&&($stInfo[3] <= $stCheck[1])&&($stInfo[4] >= $stCheck[2]))
				{
					$nIdx1 = 1;
					last;
				}
			}
			if(($nIdx1 == 0)&&($stExistPEP{$stInfo[8]} eq "OK"))
			{
				for(my $i=0; $i<@stGffInfo; $i++)
				{
					my @stColInfo=split(/\t/,$stGffInfo[$i]);
					$stGffInfo[$i] =~ s/	$stColInfo[1]	/	$stCol1_Name	/;
					$stFinalGff{$stInfo[0]."_$stInfo[3]_$stInfo[4]"} = $stFinalGff{$stInfo[0]."_$stInfo[3]_$stInfo[4]"}."$stGffInfo[$i]\n";
				}
			}
			$nIdx1 = 0;
	}
	close DATA;
	return %stFinalGff;
}
close(FH);

my $stGff3 = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.Merge.gff3";
my $stOut = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME";
my $stIDInfo = "$RUNNING_PATH/$OUTPUT_PREFIX.changeID.Info";
my $stGenome = "$TARGET_GENOME";

my $nCnt = 0;
my $nIdx = 0;
my %stList = {};
my ($stChr,$stTempChr,$stID);
open(DATA, "$stIDInfo");
while(my $stLine = <DATA>)
{
	chomp($stLine);
	my @stInfo = split /[\t]+/, $stLine;
	$stList{$stInfo[1]} = $stInfo[0];
}
close(DATA);
open(DATA, $stGff3);
open(GFF, ">$stOut.TGFam.gff3");
while(my $stLine = <DATA>)
{
	chomp($stLine);
	my @stInfo = split /[\t]+/, $stLine;
	if($stInfo[2] eq "gene")
	{
		$stChr = $stList{$stInfo[0]};
		if($stTempChr ne "")
		{
			if ($stTempChr ne $stChr)
			{
				$nCnt = 0;
			}
		}

		if ($stInfo[1] eq "PUBLIC")
		{
			$nIdx = 0;
		}
		else
		{
			$nIdx = 1;
			$nCnt++;
			$stID = "ID=$stOut.$stChr.$nCnt";
		}
		$stTempChr = $stChr;
	}

	if ($stInfo[0] ne "")
	{
		$stInfo[0] = $stChr;
		if($nIdx == 1)
		{
			$stInfo[$#stInfo] = $stID;
		}
		$stLine = join("	",@stInfo);
	}
	print GFF "$stLine\n";
}
close(DATA);
close(GFF);

my $stScaffold = "$TARGET_GENOME";
my $stGff = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.TGFam.gff3";
my $stOut = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.TGFam.CDS.fa";

my $stName = "";
my $stGCDS = "";
my $charStrand = "";
my %stScaffold = {};
my $nCnt = 0;
my $stSeq = "";

open(DATA, "$stScaffold");
while(my $stLine = <DATA>)
{
	chomp($stLine);
	if($stLine =~ /^>([^\s]+)/)
	{
		if($stSeq ne "")
		{
			$stScaffold{$stName} = uc($stSeq);
		}
		$stName = $1;
		$stSeq = "";
	}
	else
	{
		$stSeq = $stSeq.$stLine;
	}
}
close(DATA);

$stScaffold{$stName} = uc($stSeq);
$stSeq = "";

open(DATA, "$stGff");
open(OUT, ">$stOut");
while(my $stLine = <DATA>)
{
	chomp($stLine);

	my @stInfo = split /[\s\t]+/, $stLine;


	if($stInfo[2] =~ /gene/)
	{
		if($stInfo[8] =~ /ID=([^\s]+)/)
		{
			$stName = $1;
			$charStrand = $stInfo[6];
			$nCnt++;
		}
	}

	next if ($stInfo[2] =~ /mRNA/);
	next if ($stInfo[2] =~ /exon/);
	next if ($stInfo[2] =~ /gene/);

	if($stLine ne "")
	{

		my $nExonLen = $stInfo[4] - $stInfo[3] + 1;
		if($nExonLen == 0)
		{
			#			print "Error	$stName\n";
		}
		my $nExonStart = $stInfo[3] - 1;
		my $stScaffold = $stScaffold{$stInfo[0]};

		my $stExon = substr($stScaffold,$nExonStart,$nExonLen);
		$stExon = uc($stExon);
		$stGCDS = $stGCDS.$stExon;
	}
	else
	{
		if($charStrand eq "-")
		{
			$stGCDS = reverse($stGCDS);
			$stGCDS =~ tr/ACGT/TGCA/;
		}
		if(length($stGCDS)>2)
		{
			print OUT ">$stName\n$stGCDS\n";
		}
		$stGCDS = "";
	}
}
close(DATA);

if($stGCDS ne "")
{
	if($charStrand eq "-")
	{
		$stGCDS = reverse($stGCDS);
		$stGCDS =~ tr/ACGT/TGCA/;
	}
	if(length($stGCDS)>2)
	{
		print OUT ">$stName\n$stGCDS\n";
	}
}

#print "$nCnt\n";


my $stAssembly = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.TGFam.CDS.fa";
my $stCodon = "$TGFAM_SCRIPTS_PATH/CodonUsage";
my $stMergeOut = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.TGFam.PEP.fa";

my $stName = "";
my $stSeq = "";
my %stCodon = {};
my $nCut = 0;
my $stTemp = "tempaaaaa";

open(OUT, ">$stTemp");
print OUT ">AA";
close(OUT);

system ("cat $stAssembly $stTemp > NewAssembly.fa");

$stAssembly = "NewAssembly.fa";
open(DATA, "$stCodon");
while(my $stLine = <DATA>)
{
	chomp($stLine);
	my @stInfo = split /[\s\t]+/, $stLine;
	$stCodon{$stInfo[0]} = $stInfo[1];
}
close(DATA);

open(DATA, "$stAssembly");
open(OUT, ">$stMergeOut");

while(my $stLine = <DATA>)
{
	chomp($stLine);
	if($stLine =~ /^>/)
	{
		my @stTrans;
		my @stSubSeq;
		my $stTrans1;
		if($stName ne "")
		{
			for(my $i=0; $i<length($stSeq); $i= $i+3)
			{
				if($i+3<=length($stSeq))
				{
					$stSubSeq[0] = substr($stSeq,$i,3);
					$stSubSeq[1] = substr($stSeq,$i+1,3);
					$stSubSeq[2] = substr($stSeq,$i+2,3);
					if($stCodon{$stSubSeq[0]} eq "")
					{
						$stCodon{$stSubSeq[0]} = "X";
					}
					if($stCodon{$stSubSeq[1]} eq "")
					{
						$stCodon{$stSubSeq[1]} = "X";
					}
					if($stCodon{$stSubSeq[2]} eq "")
					{
						$stCodon{$stSubSeq[2]} = "X";
					}

					$stTrans[0] = $stTrans[0].$stCodon{$stSubSeq[0]};
					$stTrans[1] = $stTrans[1].$stCodon{$stSubSeq[1]};
					$stTrans[2] = $stTrans[2].$stCodon{$stSubSeq[2]};
				}
			}
			for(my $i=0; $i<3; $i++)
			{
				if ($stTrans[$i] !~ /\*./)
				{
					$stTrans1 = $stTrans[$i];
					last;
				}
			}
			if((length($stSeq) % 3 == 0 )&&($stTrans[0] !~ /\*./))
			{
				$stTrans1 = $stTrans[0];
			}

			if($stTrans1 !~ /\*./)
			{
				if(($stTrans1 =~ /^M/) && ($stTrans1 =~ /\*$/))
				{
					$stName = $stName."	FULL";
				}
				else
				{
					$stName = $stName."	Partial";
				}
				print OUT ">$stName\n$stTrans1\n";

			}
			else
			{
				print ">$stName\n$stTrans1\n";
			}
		}
		$stName = $stLine;
		$stName =~ s/>//g;
		$stSeq = "";
	}
	else
	{
		$stLine =~ tr/acgt/ACGT/;
		$stSeq = $stSeq.$stLine;
	}
}
close(DATA);
close(OUT);

my $stData = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.TGFam.PEP.fa";

open(FH2, ">$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.TGFam.PEP.fa.noStop");
open(DATA, "$stData");
while(my $stLine = <DATA>)
{
	chomp($stLine);
	if($stLine =~ /^>/)
	{
		print FH2 "$stLine\n";
	}
	else
	{
		$stLine =~ s/\*$//g;
		print FH2 "$stLine\n";
	}
}
close DATA;
close (FH2);

system("rm -rf $OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.TGFam.PEP.fa.noStop.tsv");
system("$IPRSCAN_PATH/interproscan.sh -appl pfam -i $OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.TGFam.PEP.fa.noStop -f tsv");


### system("perl /var2/TGFam/scripts/RunHmm.Generate.tsv.pl $HMMER_BIN_PATH /var2/script/PF00096.hmm $OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.TGFam.PEP.fa.noStop $OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.TGFam.PEP.fa.noStop.tsv 08-02-2019 $OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.TGFam");

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
my $now = sprintf("%04d-%02d-%02d", $year + 1900, $mon + 1, $mday);

my $stDate = $now;
my ($stDomain,$nIdx,$stPFam,$stGene, $nDomainLen);

if ($HMM_MATRIX_NAME eq "") {
}
else {

system("$HMMER_BIN_PATH/hmmsearch $HMM_MATRIX_NAME $OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.TGFam.PEP.fa.noStop > $OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.TGFam.search.out");

open(DATA, "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.TGFam.search.out");
open(OUT, ">$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.TGFam.SelfBuild.Hmm.out");
while(my $stLine = <DATA>)
{
	chomp($stLine);
	if ($stLine =~ /^Query:[\t\s]+([^\s]+)[\s\t]+\[M=([0-9]+)\]/)
	{
		$stDomain = $1;
		$nDomainLen = $2-1;
		$nIdx = 0;
	}
	elsif ($stLine =~ /^Accession:[\t\s]+([^\s]+)/)
	{
		$stPFam = $1;
	}
	elsif($stLine =~ />> ([^\s]+)/)
	{
		$stGene = $1;
	}
	elsif($stLine =~ /---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----/)
	{
		$nIdx++;
	}
	elsif (($nIdx == 1)&&($stLine ne ""))
	{
		my @stList = split /[\s\t]+/, $stLine;
		my $stTemp = join (",", @stList);
		next if ($stList[5]>$HMM_CUTOFF);
		print OUT "$stGene	Additional_Search	$nDomainLen	HMM	$stPFam	$stDomain	$stList[10]	$stList[11]	$stList[5]	T	$stDate\n";
		
	}
	elsif($stLine eq "")
	{
		$nIdx = 0;
	}
}
close(DATA);
close(OUT);

system("cat $OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.TGFam.SelfBuild.Hmm.out >> $OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.TGFam.PEP.fa.noStop.tsv");

}

system("rm -rf $OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.TGFam.PEP.fa.noStop");
system("mv $RUNNING_PATH/$MERGING_ANALYSIS_PATH/$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.TGFam.* $RUNNING_PATH/$FINAL_ANALYSIS_PATH");
system("rm -rf NewAssembly.fa");

print "\n########## 5.Generating_FinalGeneModel.pl is finished ##########\n\n";
