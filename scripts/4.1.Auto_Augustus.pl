use lib $ARGV[0];

BEGIN {
	require $ARGV[1];
		$ARGV[1]->import( qw(TGFAM_CONFIG) );
	require $ARGV[2];
		$ARGV[2]->import( qw(TGFAM_PATH) );
	}

use strict;

print "############### 4.1.Auto_Augustus.pl is started ################\n";

system("rm -rf $RUNNING_PATH/$AUGUSTUS_ANALYSIS_PATH/*");
system("cp -rf $RUNNING_PATH/$PROTEIN_MAPPING_ANALYSIS_PATH/temp_protein.$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME/$OUTPUT_PREFIX.RefPEP.Consensus.Repre.* $RUNNING_PATH/$MERGING_ANALYSIS_PATH"); 
system("cp -rf $RUNNING_PATH/$ISGAP_ANALYSIS_PATH/temp_ISGAP.$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME/$OUTPUT_PREFIX.Final.GeneModel.* $RUNNING_PATH/$MERGING_ANALYSIS_PATH"); 
system("cp -rf $RUNNING_PATH/$OUTPUT_PREFIX.geneIncluding*.fasta $RUNNING_PATH/$MERGING_ANALYSIS_PATH"); 

my $stPrefix = $OUTPUT_PREFIX;
my $stDomain = $REPRESENTATIVE_DOMAIN_NAME;
my $stPfamIDs = $TARGET_DOMAIN_ID;
my $charIdx = ""; # if you want to use data from PM and ISGAP without original PEP, you shoud use this option "N"

my @stPrefix = split /,/, $stPrefix;
my @stPfamID = split(/,/, $stPfamIDs);

for(my $i=0; $i<@stPrefix; $i++)
{
	#	print "$stPrefix[$i]\n";
	###	tsv file filtering by NB-ARC domain and sorting by e-value	###
	for(my $j=0; $j<@stPfamID; $j++)
	{
		if($j == 0)
		{
			
				my $stData = "$stPrefix[$i].Final.GeneModel.Whole.Repre.sort.PEP.fa.noStop.tsv";
				my $stPfam = "$stPfamID[$j]";
				my @NB=();
				my $num=0;

				open (FH, ">>$stPrefix[$i].Final.GeneModel.Whole.Repre.sort.PEP.fa.noStop.tsv.FilteredBy"."$stDomain"."domain.all") || die ("Can't open myfile");

				open(DATA, "$stData");
				while(my $stLine = <DATA>)
				{
					chomp($stLine);
					my @stList = split /[\t]/, $stLine;
					if(($stList[0] ne "")&&($stList[4] eq "$stPfam"))
					{
						push(@NB, "$stLine");	
					}
				}
				close DATA;

				my @sort= sort{(split /[\t]/, $a)[8] <=> (split /[\t]/, $b)[8]}@NB;
				for(my $i=0; $i<@sort; $i++) ### 30 or @sort
				{
					print FH "$sort[$i]\n";
				}
				close(FH);

				my $stData = "$stPrefix[$i].RefPEP.Consensus.Repre.PEP.fa.noStop.tsv";
				my $stPfam = "$stPfamID[$j]";
				my @NB=();
				my $num=0;

				open (FH1, ">>$stPrefix[$i].RefPEP.Consensus.Repre.PEP.fa.noStop.tsv.FilteredBy"."$stDomain"."domain.all") || die ("Can't open myfile");

				open(DATA, "$stData");
				while(my $stLine = <DATA>)
				{
					chomp($stLine);
					my @stList = split /[\t]/, $stLine;
					if(($stList[0] ne "")&&($stList[4] eq "$stPfam"))
					{
						push(@NB, "$stLine");	
					}
				}
				close DATA;

				my @sort= sort{(split /[\t]/, $a)[8] <=> (split /[\t]/, $b)[8]}@NB;
				for(my $i=0; $i<@sort; $i++) ### 30 or @sort
				{
					print FH1 "$sort[$i]\n";
				}
				close(FH1);
		}
		else
		{

				my $stData = "$stPrefix[$i].Final.GeneModel.Whole.Repre.sort.PEP.fa.noStop.tsv";
				my $stPfam = "$stPfamID[$j]";
				my @NB=();
				my $num=0;

				open (FH, ">>$stPrefix[$i].Final.GeneModel.Whole.Repre.sort.PEP.fa.noStop.tsv.FilteredBy"."$stDomain"."domain.all") || die ("Can't open myfile");

				open(DATA, "$stData");
				while(my $stLine = <DATA>)
				{
					chomp($stLine);
					my @stList = split /[\t]/, $stLine;
					if(($stList[0] ne "")&&($stList[4] eq "$stPfam"))
					{
						push(@NB, "$stLine");	
					}
				}
				close DATA;

				my @sort= sort{(split /[\t]/, $a)[8] <=> (split /[\t]/, $b)[8]}@NB;
				for(my $i=0; $i<@sort; $i++) ### 30 or @sort
				{
					print FH "$sort[$i]\n";
				}
				close(FH);

				my $stData = "$stPrefix[$i].RefPEP.Consensus.Repre.PEP.fa.noStop.tsv";
				my $stPfam = "$stPfamID[$j]";
				my @NB=();
				my $num=0;

				open (FH1, ">>$stPrefix[$i].RefPEP.Consensus.Repre.PEP.fa.noStop.tsv.FilteredBy"."$stDomain"."domain.all") || die ("Can't open myfile");

				open(DATA, "$stData");
				while(my $stLine = <DATA>)
				{
					chomp($stLine);
					my @stList = split /[\t]/, $stLine;
					if(($stList[0] ne "")&&($stList[4] eq "$stPfam"))
					{
						push(@NB, "$stLine");	
					}
				}
				close DATA;

				my @sort= sort{(split /[\t]/, $a)[8] <=> (split /[\t]/, $b)[8]}@NB;
				for(my $i=0; $i<@sort; $i++) ### 30 or @sort
				{
					print FH1 "$sort[$i]\n";
				}
				close(FH1);

		}
	}

	###	NB-ARC domain extraction	###

		my $stData = "$stPrefix[$i].Final.GeneModel.Whole.Repre.sort.PEP.fa.noStop.tsv.FilteredBy"."$stDomain"."domain.all"; # NB-ARC.tsv file such as 1_0.at_10.fasta.onlyNB-ARC.sorted.tsv
		my $stPEP = "$stPrefix[$i].Final.GeneModel.Whole.Repre.sort.PEP.fa";
		my %ID;
		my $nbID="";
		my $count=0;
		my $stSeq="";


		open (FH2, ">$stPrefix[$i].Final.GeneModel.Whole.Repre.sort.PEP.fa.filter") || die ("Can't open myfile");
		open(DATA, "$stData");
		while(my $stLine = <DATA>)
		{
				chomp($stLine);
				my @stList = split /[\t]/, $stLine;
				$ID{$stList[0]}=$ID{$stList[0]}."$stList[6]	$stList[7],";
		}
		close (DATA);

		open(DATA2, "$stPEP");
		while(my $stLine2 = <DATA2>)
		{
				chomp($stLine2);
				if($stLine2 =~ />([^\s]+)/)
				{
				if($stSeq ne "")
				{
						if($count == 1)
						{
								print FH2 ">$nbID\n$stSeq\n";
					}
					$count=0;
					$stSeq="";
				}
						$nbID= $1;
						if($ID{$nbID} ne "")
						{
								$count++;                       
						}
				}       
			else
			{
				$stSeq=$stSeq."$stLine2";
			}
		}
		if($stSeq ne "")
		{
				if($count == 1)
				{
						print FH2 ">$nbID\n$stSeq\n";
			}
			$count=0;
			$stSeq="";
		}
		close (DATA2);
		close(FH2);

		my $stData = "$stPrefix[$i].RefPEP.Consensus.Repre.PEP.fa.noStop.tsv.FilteredBy"."$stDomain"."domain.all"; # NB-ARC.tsv file such as 1_0.at_10.fasta.onlyNB-ARC.sorted.tsv
		my $stPEP = "$stPrefix[$i].RefPEP.Consensus.Repre.PEP.fa";
		my %ID;
		my $nbID="";
		my $count=0;
		my $stSeq="";


		open (FH3, ">$stPrefix[$i].RefPEP.Consensus.Repre.PEP.fa.filter") || die ("Can't open myfile");
		open(DATA, "$stData");
		while(my $stLine = <DATA>)
		{
				chomp($stLine);
				my @stList = split /[\t]/, $stLine;
				$ID{$stList[0]}=$ID{$stList[0]}."$stList[6]	$stList[7],";
		}
		close (DATA);

		open(DATA2, "$stPEP");
		while(my $stLine2 = <DATA2>)
		{
				chomp($stLine2);
				if($stLine2 =~ />([^\s]+)/)
				{
				if($stSeq ne "")
				{
						if($count == 1)
						{
								print FH3 ">$nbID\n$stSeq\n";
					}
					$count=0;
					$stSeq="";
				}
						$nbID= $1;
						if($ID{$nbID} ne "")
						{
								$count++;                       
						}
				}       
			else
			{
				$stSeq=$stSeq."$stLine2";
			}
		}
		if($stSeq ne "")
		{
				if($count == 1)
				{
						print FH3 ">$nbID\n$stSeq\n";
			}
			$count=0;
			$stSeq="";
		}
		close (DATA2);
		close(FH3);



	### filter gff, CDS ###

		my $stPEP = "$stPrefix[$i].Final.GeneModel.Whole.Repre.sort.PEP.fa.filter";
		my $stCDS = "$stPrefix[$i].Final.GeneModel.Whole.Repre.sort.CDS.fa";
		my $stGff = "$stPrefix[$i].Final.GeneModel.Whole.Repre.sort.gff3";

		my $stName;
		my %stList = {};
		my $stID = "";
		my $nCnt1=0;
		my $nCnt2=0;
		my $nIdx =0;
		my $nCnt = -1;

		open(DATA, "$stPEP");
		open(OUT, ">$stPEP.filter");
		while(my $stLine = <DATA>)
		{
				chomp($stLine);
				if ($stLine =~ /^>([^\s]+)/)
				{
						$stName = $1;
				}
			else
			{
				my $stSeq = $stLine;
				$stLine =~ s/[^*]//g;
				if (length($stLine)<2)
				{
							$stList{$stName}++;
								print OUT ">$stName\n$stSeq\n";
							
				}
			}
		}
		close(DATA);


		open(DATA, "$stCDS");
		open(OUT, ">$stCDS.filter");
		while(my $stLine = <DATA>)
		{
				chomp($stLine);
				if($stLine =~ /^>([^\s]+)/)
				{
						$stName = $1;
				}
				else
				{
						if($stList{$stName} ne "")
						{
								print OUT ">$stName\n$stLine\n";
						}
				}
		}
		close(DATA);
		close(OUT);

		open(DATA, "$stGff");
		open(OUT, ">$stGff.filter");
		while(my $stLine = <DATA>)
		{
				chomp($stLine);
				my @stInfo = split /[\s\t]+/, $stLine;
				if($stInfo[2] eq "gene")
				{
						$stID = $stInfo[$#stInfo];
						$stID =~ s/-Exonerate-.+//g;
						$stID =~ s/ID=//g;
						$nCnt++;
		#                $stID = "$stID"."_$nCnt";
						if($stList{$stID} ne "")
						{
								$nIdx = 1;
								$nCnt2++;
						}
						else
						{
								$nIdx = 0;
						}
						$nCnt1++;
				}
				if($nIdx == 1)
				{
						print OUT "$stLine\n";
				}
		}
		close(DATA);
		close(OUT);

		#		print "$nCnt1   $nCnt2\n";

		my $stPEP = "$stPrefix[$i].RefPEP.Consensus.Repre.PEP.fa.filter";
		my $stCDS = "$stPrefix[$i].RefPEP.Consensus.Repre.CDS.fa";
		my $stGff = "$stPrefix[$i].RefPEP.Consensus.Repre.gff3";

		my $stName;
		my %stList = {};
		my $stID = "";
		my $nCnt1=0;
		my $nCnt2=0;
		my $nIdx =0;
		my $nCnt = -1;

		open(DATA, "$stPEP");
		open(OUT, ">$stPEP.filter");
		while(my $stLine = <DATA>)
		{
				chomp($stLine);
				if ($stLine =~ /^>([^\s]+)/)
				{
						$stName = $1;
				}
			else
			{
				my $stSeq = $stLine;
				$stLine =~ s/[^*]//g;
				if (length($stLine)<2)
				{
							$stList{$stName}++;
								print OUT ">$stName\n$stSeq\n";
							
				}
			}
		}
		close(DATA);


		open(DATA, "$stCDS");
		open(OUT, ">$stCDS.filter");
		while(my $stLine = <DATA>)
		{
				chomp($stLine);
				if($stLine =~ /^>([^\s]+)/)
				{
						$stName = $1;
				}
				else
				{
						if($stList{$stName} ne "")
						{
								print OUT ">$stName\n$stLine\n";
						}
				}
		}
		close(DATA);
		close(OUT);

		open(DATA, "$stGff");
		open(OUT, ">$stGff.filter");
		while(my $stLine = <DATA>)
		{
				chomp($stLine);
				my @stInfo = split /[\s\t]+/, $stLine;
				if($stInfo[2] eq "gene")
				{
						$stID = $stInfo[$#stInfo];
						$stID =~ s/-Exonerate-.+//g;
						$stID =~ s/ID=//g;
						$nCnt++;
		#                $stID = "$stID"."_$nCnt";
						if($stList{$stID} ne "")
						{
								$nIdx = 1;
								$nCnt2++;
						}
						else
						{
								$nIdx = 0;
						}
						$nCnt1++;
				}
				if($nIdx == 1)
				{
						print OUT "$stLine\n";
				}
		}
		close(DATA);
		close(OUT);



	### Full, Partial divide 

		my $stData = "$stPrefix[$i].Final.GeneModel.Whole.Repre.sort.gff3.filter";


		open (FH4, ">$stPrefix[$i].Final.GeneModel.Whole.Repre.sort.gff3.filter.Info") || die ("Can't open myfile");
		open(DATA, "$stData");
		while(my $stLine = <DATA>)
		{
				chomp($stLine);
				my @stList = split /[\s]+/, $stLine;
				if($stList[2] eq "gene")
		#	if($stList[2] eq "mRNA")
				{
					print FH4 "$stLine\n";       
			}

		}
		close DATA;
		close(FH4);

		my $stData = "$stPrefix[$i].RefPEP.Consensus.Repre.gff3.filter";


		open (FH5, ">$stPrefix[$i].RefPEP.Consensus.Repre.gff3.filter.Info") || die ("Can't open myfile");
		open(DATA, "$stData");
		while(my $stLine = <DATA>)
		{
				chomp($stLine);
				my @stList = split /[\s]+/, $stLine;
				if($stList[2] eq "gene")
		#	if($stList[2] eq "mRNA")
				{
					print FH5 "$stLine\n";       
			}

		}
		close DATA;
		close(FH5);

			my $stData = "$stPrefix[$i].Final.GeneModel.Whole.Repre.sort.gff3.filter.Info"; ### final Info
			my $stPEP = "$stPrefix[$i].Final.GeneModel.Whole.Repre.sort.PEP.fa.filter";
			my $stID="";
			my %IDSeq;

			open(DATA, "$stPEP");
			open Full, ">", "$stPEP.Full";
			open Partial, ">", "$stPEP.Partial";
			while(my $stLine = <DATA>)
			{
				chomp($stLine);
				if($stLine =~ /^>([^\s]+)/)
				{
					$stID = $1;		
				}
				else
				{
					$IDSeq{$stID} = $stLine;	
					if(($stLine =~ /^M/)&&($stLine =~ /\*$/))
					{
						print Full ">$stID\n$stLine\n";
					}
					else
					{
						print Partial ">$stID\n$stLine\n";
					}
				}
			}
			close DATA;
			close Full;
			close Partial;

			open(DATA, "$stData");
			open FH, ">", "$stData.Full";
			open FH2, ">", "$stData.Partial";
			while(my $stLine = <DATA>)
			{
				chomp($stLine);
				my @stList = split /[\t]/, $stLine;
				$stList[8] =~ s/ID=//g;
				my $seq = $IDSeq{$stList[8]};

				if(($seq =~ /^M/)&&($seq =~ /\*$/))
				{
					print FH "$stLine\n";
				}
				else
				{
					print FH2 "$stLine\n";
				}
			}
			close DATA;
			close FH;
			close FH2;

			my $stData = "$stPrefix[$i].RefPEP.Consensus.Repre.gff3.filter.Info"; ### final Info
			my $stPEP = "$stPrefix[$i].RefPEP.Consensus.Repre.PEP.fa.filter";
			my $stID="";
			my %IDSeq;

			open(DATA, "$stPEP");
			open Full, ">", "$stPEP.Full";
			open Partial, ">", "$stPEP.Partial";
			while(my $stLine = <DATA>)
			{
				chomp($stLine);
				if($stLine =~ /^>([^\s]+)/)
				{
					$stID = $1;		
				}
				else
				{
					$IDSeq{$stID} = $stLine;	
					if(($stLine =~ /^M/)&&($stLine =~ /\*$/))
					{
						print Full ">$stID\n$stLine\n";
					}
					else
					{
						print Partial ">$stID\n$stLine\n";
					}
				}
			}
			close DATA;
			close Full;
			close Partial;

			open(DATA, "$stData");
			open FH, ">", "$stData.Full";
			open FH2, ">", "$stData.Partial";
			while(my $stLine = <DATA>)
			{
				chomp($stLine);
				my @stList = split /[\t]/, $stLine;
				$stList[8] =~ s/ID=//g;
				my $seq = $IDSeq{$stList[8]};

				if(($seq =~ /^M/)&&($seq =~ /\*$/))
				{
					print FH "$stLine\n";
				}
				else
				{
					print FH2 "$stLine\n";
				}
			}
			close DATA;
			close FH;
			close FH2;

	### cat Full type ###

	if($charIdx eq "N")
	{
		sleep 3;
		system("cat $stPrefix[$i].Final.GeneModel.Whole.Repre.sort.PEP.fa.filter.Full $stPrefix[$i].RefPEP.Consensus.Repre.PEP.fa.filter.Full > $stPrefix[$i].TrainingSet.ForAugustus.fasta");
	}
	else
	{
		sleep 3;
		system("cat $stPrefix[$i].geneIncluding*.fasta $stPrefix[$i].Final.GeneModel.Whole.Repre.sort.PEP.fa.filter.Full $stPrefix[$i].RefPEP.Consensus.Repre.PEP.fa.filter.Full > $stPrefix[$i].TrainingSet.ForAugustus.fasta");
	}

	### non redundant ###

		my $stData="$stPrefix[$i].TrainingSet.ForAugustus.fasta";	###	Chi.Final.GeneModel.Whole.PEP.fa


		open (FH6, ">$stPrefix[$i].TrainingSet.ForAugustus.NR.fasta") || die ("Can't open myfile");
		open(DATA, "$stData");

		my %stList;
		my $stID="";

		while(my $stLine = <DATA>)
		{
			chomp($stLine);
			if($stLine =~ /^>[^\s]+/)
			{
				$stID=$stLine;
			}
			else
			{
				$stList{$stLine}="$stID";
			}
		}
		close DATA;

		foreach my $stKey (keys %stList)
		{
			print FH6 "$stList{$stKey}\n$stKey\n";
		}
		close(FH6);

}

print "\n############### 4.1.Auto_Augustus.pl is finished ###############\n\n";
