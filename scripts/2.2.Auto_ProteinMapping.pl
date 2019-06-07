use lib $ARGV[0];

BEGIN {
	require $ARGV[1];
		$ARGV[1]->import( qw(TGFAM_CONFIG) );
	require $ARGV[2];
		$ARGV[2]->import( qw(TGFAM_PATH) );
	}

use strict;

print "############ 2.2.Auto_ProteinMapping.pl is started #############\n";

my @stGenome = glob("*.genome.fa");
my @stProtein = glob("*.protein.fa");

system("rm -rf $OUTPUT_PREFIX.RefPEP.Protein.sorted.gff3");
system("rm -rf $OUTPUT_PREFIX.RefPEP.Protein.sorted.PosiModi.gff3");
system("rm -rf Error.list");
system("rm -rf Error.Extend.list");
system("rm -rf Error.Consensus.list");
system("rm -rf $OUTPUT_PREFIX.RefPEP.Consensus.Repre.gff3");
system("rm -rf $OUTPUT_PREFIX.RefPEP.Consensus.Repre.PEP.fa.noStop");

for(my $i=0; $i<@stGenome; $i++)
{
	my $stOut = $stGenome[$i];
	$stOut =~ s/genome.fa/exonerate.out/g;
	system("$EXONERATE_BIN_PATH/exonerate --model protein2genome --percent 50 --showtargetgff yes --showquerygff yes -t $stGenome[$i] -q $stProtein[$i] --querytype protein --targettype dna --maxintron $MAX_INTRON_LENGTH --forcegtag yes > $stOut");

	#	print "$i th $stGenome[$i] is finished\n";
}

my $stData = "$OUTPUT_PREFIX.RefPEP.Protein.out"; #### "$stPrefix.RefPEP.Protein.out";
my $stMaskedGenome = "$RUNNING_PATH/$OUTPUT_PREFIX.tempID.Masked.Except$REPRESENTATIVE_DOMAIN_NAME.fasta.Main.fa";
my $stPrefix = $OUTPUT_PREFIX;         ### *'s name in *.RefPEP.Protein.out ex)ATH, CAN
my $stDomainName = $REPRESENTATIVE_DOMAIN_NAME;
my $stOut=$stData;
$stOut =~ s/out/gff3/;

sleep 3;
system("cat *.exonerate.out > $OUTPUT_PREFIX.RefPEP.Protein.out");
system("rm -rf *.genome.fa *.protein.fa *.exonerate.out");

	my $stData = "$OUTPUT_PREFIX.RefPEP.Protein.out"; #Target + Query
	my $stList = "$RUNNING_PATH/$PROTEIN_MAPPING_ANALYSIS_PATH/$stPrefix.$stDomainName.MatchedEvidenceProtein.list";
	my $fCoverage = "90";
	my $stWeightName = "REF";
	my $stQueryPrefix = "tempCh[0-9]+";
	my $stTargetID = "";
	my $nIdx = 0;
	my $nFrame = 0;
	my $nMatchedLen = 0;
	my $stID = "";
	my $nCnt = 0;
	my %stList = {};
	my $nGene = 0;
	my $nExtractGene = 0;
	my $nGeneID = 0;
	my $stAttribute = "";

	open(DATA, "$stList");
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		my @stInfo = split /[\s\t]+/, $stLine;
		$stList{$stInfo[1]} = $stInfo[$#stInfo];
		
	}
	close(DATA);

	open(DATA, "$stData");
	open(OUT, ">$stOut");
	while(my $stLine =<DATA>)
	{
		chomp($stLine);
		if($stLine =~ /Query\: ([^\s]+)/)
		{
			$stID = $1;
		}
		elsif($stLine =~ /Query range\: ([0-9]+) \-\> ([0-9]+)/)
		{
			$nMatchedLen = $2 -$1 + 1;
			if(($nMatchedLen/$stList{$stID})*100 > $fCoverage)
			{
				$nCnt = 1;
			}
			else
			{
				$nCnt = 0;
			}
			$nGene++;
		}
		if($stLine =~ /Target\: ([^\s]+)/)
		{
			$stTargetID = $1;
		}
		if($stLine =~ /Target range\: ([0-9]+) \-\> ([0-9]+)/)
		{
	#		$nMatchedLen = abs($2 -$1) + 1;
	#		if(($nMatchedLen/$stList{$stTargetID})*100 > $fCoverage)
	#		{
	#			$nCnt = 1;
	#		}
			$nIdx = 1;
		}
		if($stLine =~ /vulgar/)
		{
			$nIdx = 0;
		}

		if($nIdx == 1)
		{
			if($stLine =~ /\#/)
			{
				$nCnt = 0;
				$nFrame++;
			}
		}

		if($nCnt == 1)
		{
			if($stLine =~ /^$stQueryPrefix/)
			{
				my @stInfo = split /[\t]+/, $stLine;
	#			$stInfo[3] = $stInfo[3] - 1;
	#			$stInfo[5] = 1;
				$stInfo[1] = $stWeightName;
				$stLine = join("	", @stInfo);
				if($stInfo[2] eq "gene")
				{
					$stInfo[8] =~ s/gene_id ([0-9]+) .+/ID=$stID;$stInfo[0]-Exonerate-prediction-$1/g;
					$stAttribute = $stInfo[8];
					$nGeneID = $stAttribute;
					$nGeneID =~ s/.+Exonerate-prediction-//g;
				}
				$stInfo[8] = $stAttribute;
				$stLine = join("	", @stInfo);

				next if ($nGeneID>2);

				if($stLine =~ /similarity/)
				{
					$nExtractGene++;
					print OUT "\n";
				}
				elsif(($stLine =~ /gene/)||($stLine =~ /cds/)||($stLine =~ /exon/)||($stLine =~ /mRNA/))
				{
					print OUT "$stLine\n";
				}
			}
		}
		
		
	}
	close(DATA);
	close(OUT);

	#	print "$nGene	$nExtractGene	$nFrame\n";


#system("perl ./SortGff.pl $stOut > $stPrefix.RefPEP.Protein.sorted.gff3"); ### vali -> by eyes

	my @stGene;
	my @stCDS;
	my @stInfo;
	my $stLine;
	my $nIdx1=0;
	my $stGff = $stOut;

open (FH, ">$stPrefix.RefPEP.Protein.sorted.gff3") || die ("Can't open myfile");
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
				print FH "$stGene[$i]\n";
			}
			for(my $i=0; $i<@stCDS; $i++)
			{
				my @stTemp = split /[\s\t]+/, $stCDS[$i];
				$stTemp[2] = "exon";	
				my $stExon = join ("	", @stTemp);
				print FH "$stExon\n";
				print FH "$stCDS[$i]\n";
			}
			print FH "$stLine\n";
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
				print FH "$stGene[$i]\n";
			}
			for(my $i=0; $i<@stCDS; $i++)
			{
				my @stTemp = split /[\s\t]+/, $stCDS[$i];
				$stTemp[2] = "exon";	
				my $stExon = join ("	", @stTemp);
				print FH "$stExon\n";
				print FH "$stCDS[$i]\n";
			}
			print FH "\n";
			@stCDS = ();
			@stGene = ();
	}

close(FH);

#system("perl ./Essential_gff_position_change.pl $stPrefix.RefPEP.Protein.sorted.gff3 > $stPrefix.RefPEP.Protein.sorted.PosiModi.gff3");

	my $stGFF = "$stPrefix.RefPEP.Protein.sorted.gff3";

open (FH2, ">$stPrefix.RefPEP.Protein.sorted.PosiModi.gff3") || die ("Can't open myfile");
	open(DATA, "$stGFF");
	while(my $stLine = <DATA>)
	{
			chomp($stLine);
			if($stLine ne "")
			{
					my @stList = split /[\s]+/, $stLine;
					$stList[0] =~ /([^\s]+)\_([0-9]+)\_[0-9]+/;
					my $posi = $2;
					my $modiCh = $1;
					my $str = $stList[3] + $posi - 1;
					my $end = $stList[4] + $posi - 1;
			my $ID = $stList[8];
			if($ID =~ /ID=[^\s]+\;tempCh[0-9]+(\_[0-9]+\_[0-9]+)\-Exonerate\-prediction\-[0-9]+/)
			{
				$ID =~ s/$1//g;
			}
			else
			{
				print FH2 "error	$stLine\n";
			}
					print FH2 "$modiCh	$stList[1]	$stList[2]	$str	$end	$stList[5]	$stList[6]	$stList[7]	$ID\n";
			}
			else
			{
					print FH2 "$stLine\n";
			}
	}
	close DATA;
close(FH2);

###	initial CDS, PEP	###

	my $stScaffold = "$RUNNING_PATH/$OUTPUT_PREFIX.tempID.Masked.Except$REPRESENTATIVE_DOMAIN_NAME.fasta.Main.fa";
	my $stGff = "$stPrefix.RefPEP.Protein.sorted.PosiModi.gff3";
	my $stOut = "$stPrefix.RefPEP.Mapped.CDS.fa";

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
			if($stName ne "")
			{
				$stScaffold{$stName} = $stSeq;
			}
			$stName = $1;
			#		print "$stName\n";
			$stSeq = "";
		}
		else
		{
			$stSeq = $stSeq.$stLine;
		}
	}
	close(DATA);
	$stScaffold{$stName} = $stSeq;


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
				$stName =~ s/-Exonerate-prediction.+//g;
				$stName = $stName."_$nCnt";
				$nCnt++;
			}
		}

		next if ($stInfo[2] =~ /mRNA/);
		next if ($stInfo[2] =~ /exon/);
		next if ($stInfo[2] =~ /gene/);

		if($stLine ne "")
		{

			my $nExonLen = $stInfo[4] - $stInfo[3] + 1;
			my $nExonStart = $stInfo[3] - 1;
			my $stScaffold = $stScaffold{$stInfo[0]};

			my $stExon = substr($stScaffold,$nExonStart,$nExonLen);
			$stExon = uc($stExon);
			
	#		if($charStrand eq "-")
	#		{
	#			$stExon = reverse($stExon);
	#			$stExon =~ tr/ACGT/TGCA/;
	#			$stGCDS = $stExon.$stGCDS;
	#		}
	#		else
	#		{
				$stGCDS = $stGCDS.$stExon;
	#		}
		}
		else
		{
			if($charStrand eq "-")
			{
				$stGCDS = reverse($stGCDS);
				$stGCDS =~ tr/ACGT/TGCA/;
			}
			if($stGCDS ne "")
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
		if($stGCDS ne "")
		{
			print OUT ">$stName\n$stGCDS\n";
		}
	}
	#	print "$nCnt\n";


	my $stAssembly = "$stPrefix.RefPEP.Mapped.CDS.fa";
	my $stCodon = "$TGFAM_SCRIPTS_PATH/CodonUsage";
	my $stOut = "$stPrefix.RefPEP.Mapped";

	my $stName = "";
	my $stSeq = "";
	my %stCodon = {};
	my $nCut = 20;
	my $stTemp = "tempaaaaa";


	open(OUT, ">$stTemp");
	print OUT ">AA";
	close(OUT);

	sleep 3;
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
	open(OUT, ">$stOut.PEP.fa");
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		if($stLine =~ /^>/)
		{
			my $stTrans1 = "";
			if($stName ne "")
			{
				for(my $i=0; $i<length($stSeq); $i= $i+3)
				{
					if($i+3<=length($stSeq))
					{
						my $stSubSeq = substr($stSeq,$i,3);
						if($stCodon{$stSubSeq} eq "")
						{
							$stCodon{$stSubSeq} = "X";
						}
						$stTrans1 = $stTrans1.$stCodon{$stSubSeq};
					}
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
	#				if($stName =~ /Total/)
	#				{
	#					if(($stName =~ /FULL/)||(length($stTrans1)>300))
	#					{
	#						print OUT ">$stName\n$stTrans1\n";
	#					}
	#				}
	#				else
	#				{
						print OUT ">$stName\n$stTrans1\n";
	#				}

					if($stTrans1 eq "")
					{
		#				print "$stName\n";
					}

				}
				else
				{
					open (FH3, ">>Error.list") || die ("Can't open myfile");
					print FH3 "$stName\n$stTrans1\n";
					close(FH3);
				}
			}
			$stName = $stLine;
			$stName =~ s/>//g;
			$stSeq = "";
		}
		else
		{
			$stSeq = $stSeq.$stLine;
		}
	}
	close(DATA);
	close(OUT);

	my $stPEP = "$stPrefix.RefPEP.Mapped.PEP.fa";
	my $stCDS = "$stPrefix.RefPEP.Mapped.CDS.fa";
	my $stGff = "$stPrefix.RefPEP.Protein.sorted.PosiModi.gff3";

	my $stName;
	my %stList = {};
	my $stID = "";
	my $nCnt1=0;
	my $nCnt2=0;
	my $nIdx =0;
	my $nCnt = -1;

	open(DATA, "$stPEP");
	open(OUT, "$stPEP.filter");
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
			$stID = "$stID"."_$nCnt";
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

	#print "$nCnt1	$nCnt2\n";


####     Extension 3 base pairs       ###
	my $stGff = "$stPrefix.RefPEP.Protein.sorted.PosiModi.gff3.filter";
	my $stGenome = "$RUNNING_PATH/$OUTPUT_PREFIX.tempID.Masked.Except$REPRESENTATIVE_DOMAIN_NAME.fasta.Main.fa";
	my $stOut = "$stPrefix.RefPEP";


	my $stStart = "ATG";
	my @stStop = ("TAA","TAG","TGA");
	my %stScaffold = {};
	my ($stName,$stSeq,$charStrand,$nGeneStart,$nGeneEnd,$stGffInfo,$nCnt,$stGCDS);
	my $stNextCodon="";

	open(DATA, "$stGenome");
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		if($stLine =~ /^>([^\s]+)/)
		{
			if($stName ne "")
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


	open(DATA, "$stGff");
	open(GFF, ">$stOut.Extend.gff3");
	while(my $stLine = <DATA>)
	{
		chomp($stLine);

		my @stInfo = split /[\s\t]+/, $stLine;
		$stGffInfo = $stGffInfo.$stLine."\n";

		if($stInfo[2] =~ /gene/)
		{
			if($stInfo[8] =~ /ID=([^\s]+)/)
			{
				$stName = $1;
				$charStrand = $stInfo[6];
				$stName =~ s/Exonerate-prediction.+//g;
				$stName = $stName."Score_$stInfo[5]";
				$nGeneStart = $stInfo[3];
				$nGeneEnd = $stInfo[4];
				$stSeq = $stScaffold{$stInfo[0]};
				$nCnt++;
			}
		}

		next if ($stInfo[2] =~ /mRNA/);
		next if ($stInfo[2] =~ /exon/);
		next if ($stInfo[2] =~ /gene/);

		if($stLine ne "")
		{
			my $nExonLen = $stInfo[4] - $stInfo[3] + 1;
			my $nExonStart = $stInfo[3] - 1;
			if($charStrand eq "-")
			{
				$stNextCodon = substr($stSeq, $stInfo[3]-3, 3);
			}
			else
			{
				$stNextCodon = substr($stSeq, $stInfo[4], 3);
			}
			my $stExon = substr($stSeq,$nExonStart,$nExonLen);

			$stGCDS = $stGCDS.$stExon;
		}
		else
		{
			if($charStrand eq "-")
			{
				$stGCDS = reverse($stGCDS);
				$stGCDS =~ tr/ACGT/TGCA/;

				$stNextCodon = reverse($stNextCodon);
				$stNextCodon =~ tr/ACGT/TGCA/;
			}
			my $stLastCodon = substr($stGCDS,length($stGCDS)-3,3);
	#		my $stNextCodon = substr($stSeq, $nGeneEnd, 3);
				
			if(($stLastCodon eq "$stStop[0]")||($stLastCodon eq "$stStop[1]")||($stLastCodon eq "$stStop[2]"))
			{
				#print "error	$stLine\n";
			}
			elsif(($stNextCodon eq "$stStop[0]")||($stNextCodon eq "$stStop[1]")||($stNextCodon eq "$stStop[2]"))
			{
				my $stInformation = "";
				if($charStrand eq "-")
				{
					my $nModiGeneStart = $nGeneStart - 3;
					$stGffInfo =~ s/	$nGeneStart	/	$nModiGeneStart	/g;
				}
				else
				{
					my $nModiGeneEnd = $nGeneEnd + 3;
					$stGffInfo =~ s/	$nGeneEnd	/	$nModiGeneEnd	/g;
				}
			}
			print GFF "$stGffInfo";
			$stGCDS = "";
			$stGffInfo = "";
			$stNextCodon="";
		}
	}
	close(DATA);
	close(GFF);

	my $stScaffold = "$RUNNING_PATH/$OUTPUT_PREFIX.tempID.Masked.Except$REPRESENTATIVE_DOMAIN_NAME.fasta.Main.fa";
	my $stGff = "$stPrefix.RefPEP.Extend.gff3";
	my $stOut = "$stPrefix.RefPEP.Extend.CDS.fa";

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
			if($stName ne "")
			{
				$stScaffold{$stName} = $stSeq;
			}
			$stName = $1;
			#print "$stName\n";
			$stSeq = "";
		}
		else
		{
			$stSeq = $stSeq.$stLine;
		}
	}
	close(DATA);
	$stScaffold{$stName} = $stSeq;


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
				$stName =~ s/-Exonerate-prediction.+//g;
				$stName = $stName."_$nCnt";
				$nCnt++;
			}
		}

		next if ($stInfo[2] =~ /mRNA/);
		next if ($stInfo[2] =~ /exon/);
		next if ($stInfo[2] =~ /gene/);

		if($stLine ne "")
		{

			my $nExonLen = $stInfo[4] - $stInfo[3] + 1;
			my $nExonStart = $stInfo[3] - 1;
			my $stScaffold = $stScaffold{$stInfo[0]};

			my $stExon = substr($stScaffold,$nExonStart,$nExonLen);
			$stExon = uc($stExon);
			
	#		if($charStrand eq "-")
	#		{
	#			$stExon = reverse($stExon);
	#			$stExon =~ tr/ACGT/TGCA/;
	#			$stGCDS = $stExon.$stGCDS;
	#		}
	#		else
	#		{
				$stGCDS = $stGCDS.$stExon;
	#		}
		}
		else
		{
			if($charStrand eq "-")
			{
				$stGCDS = reverse($stGCDS);
				$stGCDS =~ tr/ACGT/TGCA/;
			}
			if($stGCDS ne "")
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
		if($stGCDS ne "")
		{
			print OUT ">$stName\n$stGCDS\n";
		}
	}

	#	print "$nCnt\n";

	my $stAssembly = "$stPrefix.RefPEP.Extend.CDS.fa";
	my $stCodon = "$TGFAM_SCRIPTS_PATH/CodonUsage";
	my $stOut = "$stPrefix.RefPEP.Extend";

	my $stName = "";
	my $stSeq = "";
	my %stCodon = {};
	my $nCut = 20;
	my $stTemp = "tempaaaaa";

	open(OUT, ">$stTemp");
	print OUT ">AA";
	close(OUT);

	sleep 3;
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
	open(OUT, ">$stOut.PEP.fa");
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		if($stLine =~ /^>/)
		{
			my $stTrans1 = "";
			if($stName ne "")
			{
				for(my $i=0; $i<length($stSeq); $i= $i+3)
				{
					if($i+3<=length($stSeq))
					{
						my $stSubSeq = substr($stSeq,$i,3);
						if($stCodon{$stSubSeq} eq "")
						{
							$stCodon{$stSubSeq} = "X";
						}
						$stTrans1 = $stTrans1.$stCodon{$stSubSeq};
					}
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
	#				if($stName =~ /Total/)
	#				{
	#					if(($stName =~ /FULL/)||(length($stTrans1)>300))
	#					{
	#						print OUT ">$stName\n$stTrans1\n";
	#					}
	#				}
	#				else
	#				{
						print OUT ">$stName\n$stTrans1\n";
	#				}

					if($stTrans1 eq "")
					{
		#				print "$stName\n";
					}

				}
				else
				{
					open (FH4, ">>Error.Extend.list") || die ("Can't open myfile");
					print FH4 "$stName\n$stTrans1\n";
					close (FH4);
				}
			}
			$stName = $stLine;
			$stName =~ s/>//g;
			$stSeq = "";
		}
		else
		{
			$stSeq = $stSeq.$stLine;
		}
	}
	close(DATA);
	close(OUT);

####     Consensus gene models   ###

	my $stGff = "$stPrefix.RefPEP.Extend.gff3";
	my $stCDS = "$stPrefix.RefPEP.Extend.CDS.fa";
	my $stOut = "$stPrefix.RefPEP.Consensus";

	my %nStart = {};
	my %nEnd = {};
	my %stList = {};
	my %stResult = {};
	my %charStrand = {};
	my %stGffInfo = {};
	my @stData =();
	my @stResult = ();
	my @stFilter = ();

	my $stName="";
	my $nGene = 1;
	my $nDiff = 0;
	my $stResult = "";
	my @nBack = ();
	my $stGffID = "";
	my $stPos = "";
	my %stPos = {};
	my $nIdx = 0;
	my @nPos = ();
	my $stIdx = "";
	my $nCnt = 0;
	my %nScore = {};

	open(DATA, "$stGff");
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		my @stInfo = split /[\s\t]+/, $stLine;
		if ($stInfo[2] eq "gene")
		{
			$stInfo[$#stInfo] =~ s/-Exonerate.+//g;
			$stInfo[$#stInfo] =~ s/ID=//g;

			$stGffID = $stInfo[$#stInfo]."_$nCnt";
			
			$nScore{$stGffID} = $stInfo[5];
			$nStart{$stGffID} = $stInfo[3];
			$nEnd{$stGffID} = $stInfo[4];
			$charStrand{$stGffID} = $stInfo[6];
			$nCnt++;
		}
		elsif($stInfo[2] eq "cds")
		{
			push(@nPos, "$stInfo[3]	$stInfo[4]");
		}
		if($stInfo[1] eq "")
		{
			for(my $i=0; $i<@nPos; $i++)
			{
				my @stTemp = split /[\s\t]+/, $nPos[$i];
				for(my $j=$stTemp[0]; $j<=$stTemp[1]; $j++)
				{
					$stPos = $stPos.$j.",";
				}
			}
			@nPos = ();
			if($stPos{$stGffID} eq "")
			{
				$stPos{$stGffID} = $stPos;
			}
			$stPos = "";
		}
	}
	close(DATA);

	$nCnt = 0;

	open(DATA, "$stCDS");
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		if($stLine =~ />/)
		{
			$stName = $stLine;
			$stName =~ s/>//g;
			$stList{$stName."_$nCnt"}++;
			$nCnt++;
		}
		else
		{
			my $nLen = length($stLine);
			$stLine = uc($stLine);
	#		if($stList{$stName}<2)
	#		{
				push(@stData, "$stName	$nStart{$stName}	$nEnd{$stName}	$nLen	$stLine	$stPos{$stName}	$charStrand{$stName}");
	#		}
		}
	}
	close(DATA);

	@stData = sort
	{
			($a =~ /^.+;tempCh([0-9]+)\_[0-9]+	/)[0] <=> ($b =~ /^.+;tempCh([0-9]+)\_[0-9]+	/)[0]||
			($a =~ /^.+;tempCh[0-9]+\_[0-9]+	[0-9]+	[0-9]+	[0-9]+	[^\s]+	[^\s]+	([+-])/)[0] cmp ($b =~ /^.+;tempCh[0-9]+\_[0-9]+	[0-9]+	[0-9]+	[0-9]+	[^\s]+	[^\s]+	([+-])/)[0]||
			($a =~ /^.+;tempCh[0-9]+\_[0-9]+	([0-9]+)/)[0] <=> ($b =~ /^.+;tempCh[0-9]+\_[0-9]+	([0-9]+)/)[0]||
			($b =~ /^.+;tempCh[0-9]+\_[0-9]+	[0-9]+	([0-9]+)/)[0] <=> ($a =~ /^.+;tempCh[0-9]+\_[0-9]+	[0-9]+	([0-9]+)/)[0]
		
	}@stData;

	open(OUT, ">$stOut.Putative.SortBy.HighlyConserve.out");
	for(my $i=0; $i<@stData; $i++)
	{
		@nBack = ();
		$nIdx = 0;
		next if ($stData[$i] eq "");
		my @stInfo1 = split /[\s\t]+/, $stData[$i];

		my ($stID,$stSeq,$nScore,$nCumCnt,$stPos,$charStrand) = ($stInfo1[0],$stInfo1[4],$nScore{$stInfo1[0]},1,$stInfo1[5],$stInfo1[6]);
		my $stID1 = $stInfo1[0];

		$stSeq = uc($stSeq);
		$stID1 =~ s/.+;(.+)\_[0-9]+/$1/g;
		my $stProtein = $stInfo1[0];
		$stProtein =~ s/(.+);.+/$1/g;

		for(my $j=$i+1; $j<@stData; $j++)
		{
			next if ($stData[$j] eq "");
			my @stInfo2 = split /[\s\t]+/, $stData[$j];
			$stInfo2[4] = uc($stInfo2[4]);
			my $stID2 = $stInfo2[0];

			$stID2 =~ s/.+;(tempCh[0-9]+)\_[0-9]+/$1/g;
			my $nScore2 = $nScore{$stInfo2[0]};
			my $stProtein2 = $stInfo2[0];
			$stProtein2 =~ s/(.+);.+/$1/g;

			if($stID1 ne $stID2)
			{
				$nGene = 1;
				last;
			}
			else
			{
				if($charStrand ne $stInfo2[6])
				{
					if($stInfo1[2] <$stInfo2[1])
					{
						$nGene++;
						last;
					}
				}
				else
				{
					if ($stPos =~ /$stInfo2[5]/)
					{
						if(($stSeq =~ /^ATG/)&&(($stSeq =~ /TA[AG]$/)||($stSeq =~ /TGA$/)))
						{
	#						$stID = $stID."FULL";
							$nScore = $nScore + $nScore2;
							$stProtein = $stProtein.",$stProtein2";
							$stData[$j] = "";
						}
						elsif(($stInfo2[4] =~ /^ATG/) &&(($stInfo2[4]=~/TA[AG]$/)||($stInfo2[4]=~ /TGA$/))){}
						else
						{
							$stData[$j] = "";
							$nScore = $nScore + $nScore2;
							$stProtein = $stProtein.",$stProtein2";
						}
						$nCumCnt++;
					}
					elsif($stInfo2[5] =~ /$stPos/)
					{
						my $nTempStart = $stInfo1[1];
						my $nTempEnd = $stInfo1[2];

						$stInfo1[1] = $stInfo2[1];
						$stInfo1[2] = $stInfo2[2];

						if(($stInfo2[4] =~ /^ATG/) &&(($stInfo2[4]=~/TA[AG]$/)||($stInfo2[4]=~ /TGA$/)))
						{
							$stSeq = uc($stInfo2[4]);
	#						$stID = $stInfo2[0]."FULL";
							$stPos = $stInfo2[5];
							$stInfo1[1] = $nTempStart;
							$stInfo1[2] = $nTempEnd;
							$nScore = $nScore + $nScore2;
							$stProtein = $stProtein.",$stProtein2";
							$stData[$j] = "";
						}
						elsif(($stSeq =~ /^ATG/)&&(($stSeq =~ /TA[AG]$/)||($stSeq =~ /TGA$/)))
						{
	#						$stID = $stID."FULL";
						}
						else
						{
							$stSeq = uc($stInfo2[4]);
							$stID = $stInfo2[0];
							$stPos = $stInfo2[5];
							$stInfo1[1] = $nTempStart;
							$stInfo1[2] = $nTempEnd;
							$nScore = $nScore + $nScore2;
							$stProtein = $stProtein.",$stProtein2";
							$stData[$j] = "";
						}
						$nCumCnt++;
					}
					else
					{
						if($stInfo1[2] <$stInfo2[1])
						{
							$nGene++;
							last;
						}
					}
				}
				my @stTemp = split /,/, $stPos;
				@stTemp = sort {$a <=> $b}@stTemp;
			}
		}
		$nDiff = 0;

		my @stTemp = split /,/, $stPos;
		my @nProtein = split /,/, $stProtein;
		my $nProtein = $#nProtein + 1;
		my $nCDSLen = $#stTemp;

		@stTemp = sort{$a <=> $b}@stTemp;
		my $nStart = $stTemp[0];
		my $nEnd = $stTemp[$#stTemp];

		my $stTemp = "$stID1	$nGene	$stID	$nProtein	$nStart	$nEnd	$nCDSLen	$nScore	$stProtein	$stPos	$charStrand";
		print OUT "$stTemp\n";
	}
	close(OUT);

	my $stData = "$stPrefix.RefPEP.Consensus.Putative.SortBy.HighlyConserve.out";
	my $stPrefix = "$stPrefix";
	my $nOption = "1"; #1 : Step1, 2: Step2
	my $stOut = "$stPrefix.RefPEP.Consensus";

	my @stNewData;
	my $nLoci = 1;
	my @stTemp = ();
	my $nCnt = 0;
	my $nGeneID = 10;

	open(DATA, "$stData");
	my @stNewData = <DATA>;
	chomp(@stNewData);
	close(DATA);

	if($nOption == 1)
	{
		&GffForStep1(\@stNewData,$stOut);
	}
	elsif($nOption == 2)
	{
		&GffForStep2(\@stNewData,$stOut);
	}
	else
	{
		exit;
	}

	sub GffForStep1
	{
		my @stNewData = @{$_[0]};
		my $stOut = $_[1];

		open(OUT, ">$stOut.gff3");
		for(my $i=0; $i<@stNewData; $i++)
		{
			my @stInfo = split /[\s\t]+/, $stNewData[$i];
			my @stInfo2 = split /[\s\t]+/, $stNewData[$i+1];
			my @stPos = split /,/, $stInfo[9];
			my $stID = $stInfo[0];
			my $nTemp = 0;
			$stID =~ s/_cov[0-9]+//g;

			print OUT "$stInfo[0]	$stPrefix	gene	$stInfo[4]	$stInfo[5]	$stInfo[7]	$stInfo[10]	.	ID=$stID.$nGeneID\n";
			print OUT "$stInfo[0]	$stPrefix	mRNA	$stInfo[4]	$stInfo[5]	.	$stInfo[10]	.	ID=$stID.$nGeneID;Parents=$stInfo[8]\n";

			my $nExonStart = $stPos[0];
			my $nExonEnd = $stPos[1];

			for(my $j=0; $j<$#stPos; $j++)
			{
				if(abs($stPos[$j]-$stPos[$j+1])>1)
				{
					$nTemp++;
					if($nTemp == 2)
					{
						$nExonEnd = $nExonStart;
					}
					print OUT "$stInfo[0]	$stPrefix	exon	$nExonStart	$nExonEnd	.	$stInfo[10]	.	Parents=$stInfo[8]\n";
					print OUT "$stInfo[0]	$stPrefix	cds	$nExonStart	$nExonEnd	.	$stInfo[10]	.	Parents=$stInfo[8]\n";
					$nExonStart = $stPos[$j+1];
				}
				else
				{
					$nExonEnd = $stPos[$j+1];
					$nTemp = 0;
				}
			}
			print OUT "$stInfo[0]	$stPrefix	exon	$nExonStart	$nExonEnd	.	$stInfo[10]	.	Parents=$stInfo[8]\n";
			print OUT "$stInfo[0]	$stPrefix	cds	$nExonStart	$nExonEnd	.	$stInfo[10]	.	Parents=$stInfo[8]\n\n";

			if($stInfo[0] ne $stInfo2[0])
			{
				$nGeneID = 10;
			}
			else
			{
				$nGeneID = $nGeneID + 10;
			}

		}
		close(OUT);
	}

	sub GffForStep2
	{
		my @stNewData = @{$_[0]};
		my $stOut = $_[1];

		open(OUT, ">$stOut.gff3");
		for(my $i=0; $i<@stNewData; $i++)
		{
			my @stInfo = split /[\s\t]+/, $stNewData[$i];
			my @stInfo2 = split /[\s\t]+/, $stNewData[$i+1];
			my @stPos = split /,/, $stInfo[10];
			my $stID = $stInfo[1];
			my $nLoci = $stInfo[13];
			$nLoci =~ s/Loci_//g;
			my $nTemp = 0;
			$stID =~ s/_cov[0-9]+//g;

			if($nLoci>1)
			{
				$nGeneID = $nGeneID - 10;
			}

			print OUT "$stInfo[1]	$stPrefix	gene	$stInfo[5]	$stInfo[6]	$stInfo[8]	$stInfo[11]	.	ID=$stID.$nGeneID.$nLoci\n";
			print OUT "$stInfo[1]	$stPrefix	mRNA	$stInfo[5]	$stInfo[6]	.	$stInfo[11]	.	ID=$stID.$nGeneID;Parents=$stInfo[9]\n";

			my $nExonStart = $stPos[0];
			my $nExonEnd = $stPos[1];

			for(my $j=0; $j<$#stPos; $j++)
			{
				if(abs($stPos[$j]-$stPos[$j+1])>1)
				{
					$nTemp++;
					if($nTemp == 2)
					{
						$nExonEnd = $nExonStart;
					}
					print OUT "$stInfo[1]	$stPrefix	exon	$nExonStart	$nExonEnd	.	$stInfo[11]	.	Parents=$stInfo[9]\n";
					print OUT "$stInfo[1]	$stPrefix	cds	$nExonStart	$nExonEnd	.	$stInfo[11]	.	Parents=$stInfo[9]\n";
					$nExonStart = $stPos[$j+1];
				}
				else
				{
					$nExonEnd = $stPos[$j+1];
					$nTemp = 0;
				}
			}
			print OUT "$stInfo[1]	$stPrefix	exon	$nExonStart	$nExonEnd	.	$stInfo[11]	.	Parents=$stInfo[9]\n";
			print OUT "$stInfo[1]	$stPrefix	cds	$nExonStart	$nExonEnd	.	$stInfo[11]	.	Parents=$stInfo[9]\n\n";

			if($stInfo[1] ne $stInfo2[1])
			{
				$nGeneID = 10;
			}
			else
			{
				$nGeneID = $nGeneID + 10;
			}

		}
		close(OUT);
	}

	my $stScaffold = "$RUNNING_PATH/$OUTPUT_PREFIX.tempID.Masked.Except$REPRESENTATIVE_DOMAIN_NAME.fasta.Main.fa";
	my $stGff = "$stPrefix.RefPEP.Consensus.gff3";
	my $stOut = "$stPrefix.RefPEP.Consensus.CDS.fa";

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
			if($stName ne "")
			{
				$stScaffold{$stName} = $stSeq;
			}
			$stName = $1;
			#			print "$stName\n";
			$stSeq = "";
		}
		else
		{
			$stSeq = $stSeq.$stLine;
		}
	}
	close(DATA);
	$stScaffold{$stName} = $stSeq;


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
				$stName =~ s/-Exonerate-prediction.+//g;
				$stName = $stName."_$nCnt";
				$nCnt++;
			}
		}

		next if ($stInfo[2] =~ /mRNA/);
		next if ($stInfo[2] =~ /exon/);
		next if ($stInfo[2] =~ /gene/);

		if($stLine ne "")
		{

			my $nExonLen = $stInfo[4] - $stInfo[3] + 1;
			my $nExonStart = $stInfo[3] - 1;
			my $stScaffold = $stScaffold{$stInfo[0]};

			my $stExon = substr($stScaffold,$nExonStart,$nExonLen);
			$stExon = uc($stExon);
			
	#		if($charStrand eq "-")
	#		{
	#			$stExon = reverse($stExon);
	#			$stExon =~ tr/ACGT/TGCA/;
	#			$stGCDS = $stExon.$stGCDS;
	#		}
	#		else
	#		{
				$stGCDS = $stGCDS.$stExon;
	#		}
		}
		else
		{
			if($charStrand eq "-")
			{
				$stGCDS = reverse($stGCDS);
				$stGCDS =~ tr/ACGT/TGCA/;
			}
			if($stGCDS ne "")
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
		if($stGCDS ne "")
		{
			print OUT ">$stName\n$stGCDS\n";
		}
	}

	#	print "$nCnt\n";

	my $stAssembly = "$stPrefix.RefPEP.Consensus.CDS.fa";
	my $stCodon = "$TGFAM_SCRIPTS_PATH/CodonUsage";
	my $stOut = "$stPrefix.RefPEP.Consensus";

	my $stName = "";
	my $stSeq = "";
	my %stCodon = {};
	my $nCut = 20;
	my $stTemp = "tempaaaaa";

	open(OUT, ">$stTemp");
	print OUT ">AA";
	close(OUT);

	sleep 3;
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
	open(OUT, ">$stOut.PEP.fa");
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		if($stLine =~ /^>/)
		{
			my $stTrans1 = "";
			if($stName ne "")
			{
				for(my $i=0; $i<length($stSeq); $i= $i+3)
				{
					if($i+3<=length($stSeq))
					{
						my $stSubSeq = substr($stSeq,$i,3);
						if($stCodon{$stSubSeq} eq "")
						{
							$stCodon{$stSubSeq} = "X";
						}
						$stTrans1 = $stTrans1.$stCodon{$stSubSeq};
					}
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
	#				if($stName =~ /Total/)
	#				{
	#					if(($stName =~ /FULL/)||(length($stTrans1)>300))
	#					{
	#						print OUT ">$stName\n$stTrans1\n";
	#					}
	#				}
	#				else
	#				{
						print OUT ">$stName\n$stTrans1\n";
	#				}

					if($stTrans1 eq "")
					{
		#				print "$stName\n";
					}

				}
				else
				{
					
					open (FH5, ">>Error.Consensus.list") || die ("Can't open myfile");
					print FH5 "$stName\n$stTrans1\n";
					close(FH5);

				}
			}
			$stName = $stLine;
			$stName =~ s/>//g;
			$stSeq = "";
		}
		else
		{
			$stSeq = $stSeq.$stLine;
		}
	}
	close(DATA);
	close(OUT);



###	Representative extractiion	###
#system("perl ./ExtractGeneInfo.pl $stPrefix.RefPEP.Consensus.gff3 > $stPrefix.RefPEP.Consensus.gff3.Info");


	my $stData = "$stPrefix.RefPEP.Consensus.gff3";


open (FH6, ">$stPrefix.RefPEP.Consensus.gff3.Info") || die ("Can't open myfile");
	open(DATA, "$stData");
	while(my $stLine = <DATA>)
	{
			chomp($stLine);
			my @stList = split /[\s]+/, $stLine;
			if($stList[2] eq "gene")
	#	if($stList[2] eq "mRNA")
			{
				print FH6"$stLine\n";       
		}

	}
	close DATA;
close (FH6);

	my $stFilterPEP="$stPrefix.RefPEP.Consensus.PEP.fa";
	my $stPositionInfo="$stPrefix.RefPEP.Consensus.gff3.Info";
	my $stOut="$stPrefix.Representative.Whole.Info";

	open(DATA, "$stFilterPEP");

	my %stList;
	my $nIdx1=0;
	my $nIdx2=0;

	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		if($stLine =~ /^>([^\s]+)	([^\s]+)/)
		{
			my $stID=$1;
			my $stType=$2;
			$stID =~ s/\_[0-9]+//g;
			$stList{$stID}="	$stType";
			$nIdx1++;
		}
	}
	close DATA;

	#	print "Total Full,Partial: $nIdx1\n";

	open(DATA, "$stPositionInfo");

	my @stPos=();

	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		my @stInfo=split(/\t/,$stLine);
		$stInfo[8] =~ s/^ID=//g;
		if($stList{$stInfo[8]} =~ /^	[^\s]+/)
		{
			if($stList{$stInfo[8]} =~ /[0-9]+/) ### ???ERROR check. duplicated geneID
			{
				#				print "$stLine\n";
			}
			$stList{$stInfo[8]} = "$stInfo[5]".$stList{$stInfo[8]};
			push(@stPos,"$stLine");
		}
	}
	close DATA;
	$nIdx1=@stPos;
	@stPos=sort{	($a =~ /^tempCh([0-9]+)	/)[0] <=> ($b =~ /^tempCh([0-9]+)	/)[0] ||
			($a =~ /	gene	([0-9]+)/)[0] <=> ($b =~ /	gene	([0-9]+)/)[0]
			}@stPos;
			#	print "Total Matched Transcripts: $nIdx1\n";

	open FH, ">", "$stOut";
	open FH1, ">", "OverlapGene2.out";
	if($#stPos == 0)
	{
		my @stInfo1=split(/\t/,$stPos[0]);
		$stInfo1[8] =~ s/^ID=//g;
		my @stScoreType1=split(/\t/,$stList{$stInfo1[8]});
		print FH "$stPos[0];$stScoreType1[1];\n";
	}
	for(my $i=0; $i<@stPos-1; $i++)
	{
		my @stInfo1=split(/\t/,$stPos[$i]);
		my @stInfo2=split(/\t/,$stPos[$i+1]);
		$stInfo1[8] =~ s/^ID=//g;
		$stInfo2[8] =~ s/^ID=//g;
		if(($stInfo1[0] eq "$stInfo2[0]")&&($stInfo1[3] <= $stInfo2[3])&&($stInfo1[4] >= $stInfo2[3]))
		{
			$nIdx2++;
			print FH1 "$stPos[$i]\n$stPos[$i+1]\n\n";
			my @stScoreType1=split(/\t/,$stList{$stInfo1[8]});
			my @stScoreType2=split(/\t/,$stList{$stInfo2[8]});
			if(($stScoreType1[1] =~ /Full/)&&($stScoreType2[1] !~ /Full/))
			{
				$stPos[$i+1]=$stPos[$i];
				if($i == $#stPos-1)
				{
					#$stPos[$i+1] =~ s/_[^\s]+$//g;
					$stPos[$i+1] = $stPos[$i+1].";$stScoreType1[1];";
					print FH "$stPos[$i+1]\n";
				}
			}
			elsif(($stScoreType1[1] !~ /Full/)&&($stScoreType2[1] =~ /Full/))
			{
				if($i == $#stPos-1)
				{
					#$stPos[$i+1] =~ s/_[^\s]+$//g;
					$stPos[$i+1] = $stPos[$i+1].";$stScoreType2[1];";
					print FH "$stPos[$i+1]\n";
				}
			}
			else
			{
				if($stScoreType1[0] >= $stScoreType2[0])
				{
					$stPos[$i+1]=$stPos[$i];
				}
				if($i == $#stPos-1)
				{
					#$stPos[$i+1] =~ s/_[^\s]+$//g;
					$stPos[$i+1] = $stPos[$i+1].";$stScoreType2[1];";
					print FH "$stPos[$i+1]\n";
				}
			}
		}
		else
		{
			my @stScoreType1=split(/\t/,$stList{$stInfo1[8]});
			#$stPos[$i] =~ s/_[^\s]+$//g;
			$stPos[$i] = $stPos[$i].";$stScoreType1[1];";
			print FH "$stPos[$i]\n";
			if($i == $#stPos-1)
			{
				my @stScoreType2=split(/\t/,$stList{$stInfo2[8]});
				#$stPos[$i+1] =~ s/_[^\s]+$//g;
				$stPos[$i+1] = $stPos[$i+1].";$stScoreType2[1];";
				print FH "$stPos[$i+1]\n";
			}
		}
	}
	#	print "Overlap Transcripts: $nIdx2\n";

close (FH);
close (FH1);

	my $stGff1="$stPrefix.RefPEP.Consensus.gff3";	###	Chi.Final.GeneModel.Whole.sort.gff3
	my $stInfoData="$stPrefix.Representative.Whole.Info";	###	tem.Whole

	my %stList;
	my %dupl;
	my $stKey="";
	my $num=0;

	open(DATA, "$stGff1");
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		if($stLine ne "")
		{
			my @stInfo=split(/\t/,$stLine);
			if($stInfo[2] eq "gene")
			{		
				$stKey=$stLine;
				if($dupl{$stKey} eq "OK")
				{
					$num=0;
					next;
				}
				else
				{
					$num=0;
					$stList{$stKey}=$stList{$stKey}."$stLine\^";
					$num++;
				}
				$dupl{$stKey}="OK";
				
			}
			elsif($num == 1)
			{
				$stList{$stKey}=$stList{$stKey}."$stLine\^";
			}
		}
	}
	close DATA;

	my @stTotal=();

open (FH7, ">$stPrefix.RefPEP.Consensus.Repre.gff3") || die ("Can't open myfile");
	open(DATA, "$stInfoData");

	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		push(@stTotal,"$stLine");
	}
	close DATA;

	@stTotal=sort{	($a =~ /^tempCh([0-9]+)	/)[0] <=> ($b =~ /^tempCh([0-9]+)	/)[0] ||
			($a =~ /	gene	([0-9]+)/)[0] <=> ($b =~ /	gene	([0-9]+)/)[0]
			}@stTotal;
	for(my $i=0; $i<@stTotal; $i++)
	{
		my $stGenePos=$stTotal[$i];
		my $stType="";
		if($stGenePos =~ /;([^\s]+;)$/)
		{
			$stType=$1;
		}
		$stGenePos =~ s/;[^\s]+;$//g;
		my @stInfo=split(/\^/,$stList{$stGenePos});
		for(my $j=0; $j<@stInfo; $j++)
		{
			print FH7 "$stInfo[$j]$stType\n";
		}
		print FH7 "\n";
	}

close (FH7);

	my $stScaffold = "$RUNNING_PATH/$OUTPUT_PREFIX.tempID.Masked.Except$REPRESENTATIVE_DOMAIN_NAME.fasta.Main.fa";
	my $stGff = "$stPrefix.RefPEP.Consensus.Repre.gff3";
	my $stOut = "$stPrefix.RefPEP.Consensus.Repre.CDS.fa";

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
				print "Error	$stName\n";
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

			print OUT ">$stName\n$stGCDS\n";
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

		print OUT ">$stName\n$stGCDS\n";
	}

	#	print "$nCnt\n";

	my $stAssembly = "$stPrefix.RefPEP.Consensus.Repre.CDS.fa";
	my $stCodon = "$TGFAM_SCRIPTS_PATH/CodonUsage";
	my $stOut = "$stPrefix.RefPEP.Consensus.Repre";

	my $stName = "";
	my $stSeq = "";
	my %stCodon = {};
	my $nCut = 20;
	my $stTemp = "tempaaaaa";

	open(OUT, ">$stTemp");
	print OUT ">AA";
	close(OUT);

	sleep 3;
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
	open(OUT, ">$stOut.PEP.fa");
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		if($stLine =~ /^>/)
		{
			my $stTrans1 = "";
			if($stName ne "")
			{
				for(my $i=0; $i<length($stSeq); $i= $i+3)
				{
					if($i+3<=length($stSeq))
					{
						my $stSubSeq = substr($stSeq,$i,3);
						if($stCodon{$stSubSeq} eq "")
						{
							$stCodon{$stSubSeq} = "X";
						}
						$stTrans1 = $stTrans1.$stCodon{$stSubSeq};
					}
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
	#				if($stName =~ /Total/)
	#				{
	#					if(($stName =~ /FULL/)||(length($stTrans1)>300))
	#					{
	#						print OUT ">$stName\n$stTrans1\n";
	#					}
	#				}
	#				else
	#				{
						print OUT ">$stName\n$stTrans1\n";
	#				}

					if($stTrans1 eq "")
					{
		#				print "$stName\n";
					}

				}
				else
				{
					#print "$stName\n$stTrans1\n";
				}
			}
			$stName = $stLine;
			$stName =~ s/>//g;
			$stSeq = "";
		}
		else
		{
			$stSeq = $stSeq.$stLine;
		}
	}
	close(DATA);
	close(OUT);


my $stData = "$stPrefix.RefPEP.Consensus.Repre.PEP.fa";

open (FH8, ">$stPrefix.RefPEP.Consensus.Repre.PEP.fa.noStop") || die ("Can't open myfile");
open(DATA, "$stData");
while(my $stLine = <DATA>)
{
	chomp($stLine);
	if($stLine =~ /^>([^\s]+)/)
	{
		
		print FH8">$1\n";
	}
	else
	{
		$stLine =~ s/\*$//g;
		print FH8"$stLine\n";
	}
}
close DATA;
close (FH8);

system("rm -rf $stPrefix.RefPEP.Consensus.Repre.PEP.fa.noStop.tsv");
system("$IPRSCAN_PATH/interproscan.sh -appl pfam -i $stPrefix.RefPEP.Consensus.Repre.PEP.fa.noStop -f tsv");

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(); 
my $now = sprintf("%04d-%02d-%02d", $year + 1900, $mon + 1, $mday); 
 
my $stDate = $now; 
my ($stDomain,$nIdx,$stPFam,$stGene, $nDomainLen);

if ($HMM_MATRIX_NAME eq "") {
}
else {

system("$HMMER_BIN_PATH/hmmsearch $HMM_MATRIX_NAME $OUTPUT_PREFIX.RefPEP.Consensus.Repre.PEP.fa.noStop > $OUTPUT_PREFIX.RefPEP.Consensus.Repre.search.out");

open(DATA, "$OUTPUT_PREFIX.RefPEP.Consensus.Repre.search.out");
open(OUT, ">$OUTPUT_PREFIX.RefPEP.Consensus.Repre.SelfBuild.Hmm.out");
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

system("cat $OUTPUT_PREFIX.RefPEP.Consensus.Repre.SelfBuild.Hmm.out >> $OUTPUT_PREFIX.RefPEP.Consensus.Repre.PEP.fa.noStop.tsv");

}


system("rm -rf NewAssembly.fa");

print "\n############ 2.2.Auto_ProteinMapping.pl is finished ############\n\n";
