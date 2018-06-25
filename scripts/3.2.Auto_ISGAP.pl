use lib $ARGV[0];

BEGIN {
	require $ARGV[1];
		$ARGV[1]->import( qw(TGFAM_CONFIG) );
	require $ARGV[2];
		$ARGV[2]->import( qw(TGFAM_PATH) );
	}

use strict;

my @stGenome = glob("*.genome.fa");
my @stProtein = glob("*.protein.fa");

system("rm -rf Error.list");
system("rm -rf Error.Extend.list");
system("rm -rf Error.Final.list");
system("rm -rf Error.repre.list");
for(my $i=0; $i<@stGenome; $i++)
{
	my $stOut = $stGenome[$i];
	$stOut =~ s/genome.fa/exonerate.out/g;
	system("$EXONERATE_BIN_PATH/exonerate --model protein2genome --percent 50 --showtargetgff yes --showquerygff yes -t $stGenome[$i] -q $stProtein[$i] --querytype protein --targettype dna --maxintron $MAX_INTRON_LENGTH --forcegtag yes > $stOut");

	#	print "$i th $stGenome[$i] is finished\n";
}
my $stData = "$OUTPUT_PREFIX.RefPEP.Protein.out"; #### "$stPrefix.RefPEP.Protein.out";
my $stMaskedGenome = "$RUNNING_PATH/$OUTPUT_PREFIX.tempID.Masked.Except$REPRESENTATIVE_DOMAIN_NAME.fasta.Main.fa";  ### masking genome
my $stPrefix = $OUTPUT_PREFIX;         ### *'s name in *.RefPEP.Protein.out ex)ATH, CAN
my $nCPU = $THREADS;
my $stOut=$stData;
$stOut =~ s/out/gff3/;

sleep 3;

system("cat *.exonerate.out > $stData");
system("rm -rf *.genome.fa *.protein.fa *.exonerate.out");
	my $stData = "$OUTPUT_PREFIX.RefPEP.Protein.out"; #Target + Query
	my $stList = "$RUNNING_PATH/$ISGAP_ANALYSIS_PATH/$stPrefix.MatchedEvidenceProtein.list";
	my $fCoverage = "70";
	my $stWeightName = "REF";
	my $stQueryPrefix = "[0-9]+";
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
		$stInfo[1] =~ s/Transcript_//g;
		$stInfo[1] =~ s/\//./g;
		$stInfo[1] =~ s/Confidence_[\.0-9]+_//g;
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

	my $stGff = $stOut;
	my @stGene;
	my @stCDS;
	my $nIdx1=0;

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

###	initial CDS, PEP	###
	my $stScaffold = "$RUNNING_PATH/$ISGAP_ANALYSIS_PATH/$stPrefix.transcripts.fa";
	my $stGff = "$stPrefix.RefPEP.Protein.sorted.gff3";
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
					open (FH2, ">>Error.list") || die ("Can't open myfile");
					print FH2 "$stName\n$stTrans1\n";
					close(FH2);
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
	my $stGff = "$stPrefix.RefPEP.Protein.sorted.gff3";

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

#	print "$nCnt1	$nCnt2\n";


####     Extension       ###
	my $stGff = "$stPrefix.RefPEP.Protein.sorted.gff3.filter";
	my $stGenome = "$RUNNING_PATH/$ISGAP_ANALYSIS_PATH/$stPrefix.transcripts.fa";
	my $stOut = "$stPrefix.RefPEP";


	my $stStart = "ATG";
	my @stStop = ("TAA","TAG","TGA");
	my %stScaffold = {};
	my ($stName,$stSeq,$charStrand,$nGeneStart,$nGeneEnd,$stGffInfo,$nCnt,$stGCDS);

	open(DATA, "$stGenome");
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		$stLine =~ s/Transcript_//g;
		$stLine =~ s/\//./g;
		$stLine =~ s/Confidence_[0-9\.]+_//g;
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
	open(OUT, ">$stOut.Extend.Info");
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

			my $stExon = substr($stSeq,$nExonStart,$nExonLen);
			$stExon = uc($stExon);
			$stGCDS = $stGCDS.$stExon;
			
			if($charStrand eq "-")
			{
	#			$stGCDS = $stExon.$stGCDS;
			}
			else
			{
	#			$stGCDS = $stGCDS.$stExon;
			}
		}
		else
		{
			if($charStrand eq "-")
			{
				$stGCDS = reverse($stGCDS);
				$stGCDS =~ tr/ACGT/TGCA/;
			}

			my $stLastCodon = substr($stGCDS,length($stGCDS)-3,3);
			my $stStartCodon = substr($stGCDS,0,3);

			if(($stStartCodon eq "ATG")&&(($stLastCodon eq "$stStop[0]")||($stLastCodon eq "$stStop[1]")||($stLastCodon eq "$stStop[2]")))
			{
				print OUT "FULL	No_Change	$stName	$nGeneStart	$nGeneEnd	$nGeneStart	$nGeneEnd\n";
			}
			elsif($stGCDS =~ /^ATG/)
			{
				my $stInformation = "";
				if($charStrand eq "+")
				{
					my $nModiGeneEnd = 0;
					for(my $i=$nGeneEnd;$i<=length($stSeq); $i=$i+3)
					{
						my $stCodon = substr($stSeq,$i,3);
						$nModiGeneEnd = $i+3;

						if($nModiGeneEnd>length($stSeq))
						{
							last;
						}

						if(($stCodon eq "$stStop[0]")||($stCodon eq "$stStop[1]")||($stCodon eq "$stStop[2]"))
						{
							$stGffInfo =~ s/	$nGeneEnd	/	$nModiGeneEnd	/g;
							$stInformation = "FULL	Add_Stop	$stName	$nGeneStart	$nGeneEnd	$nGeneStart	$nModiGeneEnd\n";
							print OUT "$stInformation\n";
							last;
						}
					}
					if($stInformation eq "")
					{
						$nModiGeneEnd = $nModiGeneEnd - 3;
						$stGffInfo =~ s/	$nGeneEnd	/	$nModiGeneEnd	/g;
						$stInformation = "Partial	Extend_End	$stName	$nGeneStart	$nGeneEnd	$nGeneStart	$nModiGeneEnd\n";
						print OUT "$stInformation\n";
					}
				}
				else
				{
					my $nModiGeneStart = 0;
					for(my $i=$nGeneStart; $i>=0; $i=$i-3)
					{
						$nModiGeneStart = $i-3;
						if($nModiGeneStart<1)
						{
							last;
						}
						my $stCodon = substr($stSeq,$nModiGeneStart-1,3);
						$stCodon = reverse($stCodon);
						$stCodon =~ tr/ACGT/TGCA/;

						if(($stCodon eq "$stStop[0]")||($stCodon eq "$stStop[1]")||($stCodon eq "$stStop[2]"))
						{
							$stGffInfo =~ s/	$nGeneStart	/	$nModiGeneStart	/g;
							$stInformation = "FULL	Add_Stop	$stName	$nGeneStart	$nGeneEnd	$nModiGeneStart	$nGeneEnd\n";
							print OUT "$stInformation\n";
							last;
						}

					}
					if($stInformation eq "")
					{
						$nModiGeneStart = $nModiGeneStart + 3;
						$stGffInfo =~ s/	$nGeneStart	/	$nModiGeneStart	/g;
						$stInformation = "Partial	Extend_End	$stName	$nGeneStart	$nGeneEnd	$nModiGeneStart	$nGeneEnd\n";
						print OUT "$stInformation\n";
					}

				}
			}
			elsif(($stLastCodon eq "$stStop[0]")||($stLastCodon eq "$stStop[1]")||($stLastCodon eq "$stStop[2]"))
			{
				my $stInformation = "";
				if($charStrand eq "-")
				{
					my $nInitialGeneEnd = $nGeneEnd;
					my $nModiGeneEnd = 0;
					for(my $i=$nGeneEnd;$i<=length($stSeq); $i=$i+3)
					{
						my $stCodon = substr($stSeq,$i,3);
						$stCodon = reverse($stCodon);
						$stCodon =~ tr/ACGT/TGCA/;
						$nModiGeneEnd = $i+3;

						if($nModiGeneEnd>length($stSeq))
						{
							last;
						}
						if(($stCodon eq "$stStop[0]")||($stCodon eq "$stStop[1]")||($stCodon eq "$stStop[2]"))
						{
							last;
						}


						if($stCodon eq "ATG")
						{
							$stGffInfo =~ s/	$nGeneEnd	/	$nModiGeneEnd	/g;
							$nGeneEnd = $nModiGeneEnd;
							$stInformation = "FULL	Add_Start	$stName	$nGeneStart	$nInitialGeneEnd	$nGeneStart	$nModiGeneEnd\n";
	#						last;
						}
					}
					if($stInformation ne "")
					{
						print OUT "$stInformation\n";
					}
					else
					{
						$nModiGeneEnd = $nModiGeneEnd - 3;
						$stGffInfo =~ s/	$nGeneEnd	/	$nModiGeneEnd	/g;
						$stInformation = "Partial	Extend_Start	$stName	$nGeneStart	$nGeneEnd	$nGeneStart	$nModiGeneEnd\n";
						print OUT "$stInformation\n";
					}
				}
				else
				{
					my $nInitialGeneStart = $nGeneStart;
					my $nModiGeneStart = 0;
					for(my $i=$nGeneStart; $i>=0; $i=$i-3)
					{
						$nModiGeneStart = $i-3;
						if($nModiGeneStart<1)
						{
							last;
						}
						my $stCodon = substr($stSeq,$nModiGeneStart-1,3);

						if(($stCodon eq "$stStop[0]")||($stCodon eq "$stStop[1]")||($stCodon eq "$stStop[2]"))
						{
							last;
						}


						if($stCodon eq "ATG")
						{
							$stGffInfo =~ s/	$nGeneStart	/	$nModiGeneStart	/g;
							$nGeneStart = $nModiGeneStart;
							$stInformation = "FULL	Add_Start	$stName	$nInitialGeneStart	$nGeneEnd	$nModiGeneStart	$nGeneEnd\n";
	#						last;
						}

					}
					if ($stInformation ne "")
					{
						print OUT "$stInformation\n";
					}
					else
					{
						$nModiGeneStart = $nModiGeneStart + 3;
						$stGffInfo =~ s/	$nGeneStart	/	$nModiGeneStart	/g;
						$stInformation = "Partial	Extend_Start	$stName	$nGeneStart	$nGeneEnd	$nModiGeneStart	$nGeneEnd\n";
						print OUT "$stInformation\n";
					}

				}
			}
			else
			{
				my ($nModiGeneStart,$nModiGeneEnd,$nEIdx,$nSIdx,$nInitialStart,$nInitialEnd,$stInformation) =(0,0,0,0,$nGeneStart,$nGeneEnd,"");
				if($charStrand eq "+")
				{
					for(my $i=$nGeneEnd;$i<=length($stSeq); $i=$i+3)
					{
						my $stCodon = substr($stSeq,$i,3);
						$nModiGeneEnd = $i+3;

						if($nModiGeneEnd>length($stSeq))
						{
							last;
						}

						if(($stCodon eq "$stStop[0]")||($stCodon eq "$stStop[1]")||($stCodon eq "$stStop[2]"))
						{
							$stGffInfo =~ s/	$nGeneEnd	/	$nModiGeneEnd	/g;
							$nEIdx = 1;
							last;
						}
					}
					for(my $i=$nGeneStart; $i>=0; $i=$i-3)
					{
						$nModiGeneStart = $i-3;
						if($nModiGeneStart<1)
						{
							last;
						}
						my $stCodon = substr($stSeq,$nModiGeneStart-1,3);
						if(($stCodon eq "$stStop[0]")||($stCodon eq "$stStop[1]")||($stCodon eq "$stStop[2]"))
						{
							last;
						}


						if($stCodon eq "ATG")
						{
							$stGffInfo =~ s/	$nGeneStart	/	$nModiGeneStart	/g;
							$nGeneStart = $nModiGeneStart;
							$nSIdx = 1;
	#						last;
						}

					}
					$nModiGeneStart = $nGeneStart;

					if (($nEIdx == 1)&&($nSIdx!=1))
					{
						print OUT "Partial	Add_Stop	$stName	$nInitialStart	$nInitialEnd	$nInitialStart	$nModiGeneEnd\n";
					}
					elsif(($nSIdx == 1)&&($nEIdx!=1))
					{
						print OUT "Partial	Add_Start	$stName	$nInitialStart	$nInitialEnd	$nModiGeneStart	$nInitialEnd\n";
					}
					elsif(($nSIdx == 1)&&($nEIdx==1))
					{
						print OUT "FULL	Add_Start_Stop	$stName	$nInitialStart	$nInitialEnd	$nModiGeneStart	$nModiGeneEnd\n";
					}
					else
					{
	#					print OUT "Partial	No_change	$stName	$nGeneStart	$nGeneEnd	$nGeneStart	$nGeneEnd\n";
					}

				}
				else
				{
					for(my $i=$nGeneStart; $i>=0; $i=$i-3)
					{
						$nModiGeneStart = $i-3;
						if($nModiGeneStart<1)
						{
							last;
						}
						my $stCodon = substr($stSeq,$nModiGeneStart-1,3);
						$stCodon = reverse($stCodon);
						$stCodon =~ tr/ACGT/TGCA/;

						if(($stCodon eq "$stStop[0]")||($stCodon eq "$stStop[1]")||($stCodon eq "$stStop[2]"))
						{
							$stGffInfo =~ s/	$nGeneStart	/	$nModiGeneStart	/g;
							$nEIdx = 1;
							last;
						}

					}
					for(my $i=$nGeneEnd;$i<=length($stSeq); $i=$i+3)
					{
						my $stCodon = substr($stSeq,$i,3);
						$stCodon = reverse($stCodon);
						$stCodon =~ tr/ACGT/TGCA/;
						$nModiGeneEnd = $i+3;

						if(($stCodon eq "$stStop[0]")||($stCodon eq "$stStop[1]")||($stCodon eq "$stStop[2]"))
						{
							last;
						}


						if($nModiGeneEnd>length($stSeq))
						{
							last;
						}


						if($stCodon eq "ATG")
						{
							$stGffInfo =~ s/	$nGeneEnd	/	$nModiGeneEnd	/g;
							$nGeneEnd = $nModiGeneEnd;
							$nSIdx = 1;
	#						last;
						}
					}
					$nModiGeneEnd = $nGeneEnd;

					if (($nEIdx == 1)&&($nSIdx!=1))
					{
						print OUT "Partial	Add_Stop	$stName	$nInitialStart	$nInitialEnd	$nModiGeneStart	$nInitialEnd\n";
					}
					elsif(($nSIdx == 1)&&($nEIdx!=1))
					{
						print OUT "Partial	Add_Start	$stName	$nInitialStart	$nInitialEnd	$nInitialStart	$nModiGeneEnd\n";
					}
					elsif(($nSIdx == 1)&&($nEIdx==1))
					{
						print OUT "FULL	Add_Start_Stop	$stName	$nInitialStart	$nInitialEnd	$nModiGeneStart	$nModiGeneEnd\n";
					}
					else
					{
	#					print OUT "Partial	No_change	$stName	$nGeneStart	$nGeneEnd	$nGeneStart	$nGeneEnd\n";
					}
				}

			}
			print GFF "$stGffInfo";
			$stGCDS = "";
			$stGffInfo = "";
		}
	}
	close(DATA);
	close(OUT);
	close(GFF);

	my $stScaffold = "$RUNNING_PATH/$ISGAP_ANALYSIS_PATH/$stPrefix.transcripts.fa";
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
					open (FH3, ">>Error.Extend.list") || die ("Can't open myfile");
					print FH3 "$stName\n$stTrans1\n";
					close (FH3);
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
			$stLine = uc($stLine);
			my $nLen = length($stLine);
	#		if($stList{$stName}<2)
	#		{
				push(@stData, "$stName	$nStart{$stName}	$nEnd{$stName}	$nLen	$stLine	$stPos{$stName}	$charStrand{$stName}");
	#		}
		}
	}
	close(DATA);

	@stData = sort
	{
			($a =~ /^.+;([^\s]+)\_[0-9]+	/)[0] cmp ($b =~ /^.+;([^\s]+)\_[0-9]+	/)[0]||
			($a =~ /^.+;[^\s]+\_[0-9]+	[0-9]+	[0-9]+	[0-9]+	[^\s]+	[^\s]+	([+-])/)[0] cmp ($b =~ /^.+;[^\s]+\_[0-9]+	[0-9]+	[0-9]+	[0-9]+	[^\s]+	[^\s]+	([+-])/)[0]||
			($a =~ /^.+;[^\s]+\_[0-9]+	([0-9]+)/)[0] <=> ($b =~ /^.+;[^\s]+\_[0-9]+	([0-9]+)/)[0]||
			($b =~ /^.+;[^\s]+\_[0-9]+	[0-9]+	([0-9]+)/)[0] <=> ($a =~ /^.+;[^\s]+\_[0-9]+	[0-9]+	([0-9]+)/)[0]
		
	}@stData;

	foreach(@stData)
	{
	#	print "$_\n";
	}

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
			my $stID2 = $stInfo2[0];

			$stID2 =~ s/.+;(.+)\_[0-9]+/$1/g;
			my $nScore2 = $nScore{$stInfo2[0]};
			my $stProtein2 = $stInfo2[0];
			$stProtein2 =~ s/(.+);.+/$1/g;
			$stInfo2[4] = uc($stInfo2[4]);

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

my $stData = "$OUTPUT_PREFIX.RefPEP.Consensus.Putative.SortBy.HighlyConserve.out";
my $stPrefix = $OUTPUT_PREFIX;
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

	my $stPrefix = $OUTPUT_PREFIX;
	my $stScaffold = "$RUNNING_PATH/$ISGAP_ANALYSIS_PATH/$stPrefix.transcripts.fa";
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
			# print "$stName\n";
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

	my $stPrefix = $OUTPUT_PREFIX;
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
				open (FH4, ">$stPrefix.transcripts.Info") || die ("Can't open myfile");
					print FH4 "$stName\n$stTrans1\n";
				close(FH4)
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

#### Representative form/AS form in Transcript level ####
	my $stPrefix = $OUTPUT_PREFIX;
	my $stFasta="$RUNNING_PATH/$ISGAP_ANALYSIS_PATH/$stPrefix.transcripts.fa";

	open(DATA, "$stFasta");

	my @stPos=();

	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		if(($stLine =~ /^>([0-9]+) ([^\s]+) ([^\s]+)([\+\-\.]) ([0-9]+)\-([0-9]+)$/)||($stLine =~ /^>([0-9]+) ([^\s]+) ([^\s]+)([\+\-\.]) ([0-9]+).+\-([0-9]+)$/))
		{
			my $stID=$1;
			my $stASID=$2;
			my $stScaffold=$3;
			my $stStrand=$4;
			my $nStr=$5;
			my $nEnd=$6;
			if($stASID =~ /^[A-Z]+_([0-9]+\.[0-9]+)$/)
			{
				my $stChange=$1;
				$stChange =~ s/^[0]+//g;
				$stASID="CUFF.$stChange";
			}
			push(@stPos, "$stScaffold	cufflinks	gene	$nStr	$nEnd	.	$stStrand	.	ID=$stID;AS_ID=$stASID");
		}
	}
	close DATA;

	@stPos = sort{($a =~ /^tempCh([0-9]+)	/)[0] <=> ($b =~ /^tempCh([0-9]+)/)[0] ||
			  ($a =~ /([0-9]+)	[0-9]+	\.	/)[0] <=> ($b =~ /([0-9]+)	[0-9]+	\.	/)[0]
			}@stPos;

open (FH5, ">$stPrefix.transcripts.Info") || die ("Can't open myfile");
	for(my $i=0; $i<@stPos; $i++)
	{
		print FH5 "$stPos[$i]\n";
	}
close(FH5);

	my $stPrefix = $OUTPUT_PREFIX;
	my $stGff = "$stPrefix.RefPEP.Consensus.gff3";
	my $stPEP = "$stPrefix.RefPEP.Consensus.PEP.fa";
	my $stTranscriptInfo = "$stPrefix.transcripts.Info";

	my %stType = {};
	my %nLoci = {};
	my @stData;
	my %stFullLoci = {};
	my $nTemp = 0;
	my @stResult;
	my %charStrand = {};

open (FH6, ">$stPrefix.Loci.Info") || die ("Can't open myfile");
	open(DATA, "$stPEP");
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		if ($stLine =~ /^>([^\s]+)[\s\t]+([^\s]+)/)
		{
			my $stName = $1;
			my $stType = $2;
			$stName =~ s/\_[0-9]+//g;
			$stType{$stName} = $stType;
			my $stLoci = $stName;
			$stLoci =~ s/\.[0-9]+//g;
			$nLoci{$stLoci}++;
			if($stType eq "FULL")
			{
				$stFullLoci{$stLoci}++;
			}
		
		}
	}
	close(DATA);

	open(DATA, "$stTranscriptInfo");
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		my @stInfo = split /[\t]+/, $stLine;
		my $stID=$stInfo[$#stInfo];
		$stID =~ s/^ID=([0-9]+).+/$1/g;
		$charStrand{$stID} = $stInfo[6];
	}
	close(DATA);

	open (DATA, "$stGff");
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		my @stInfo = split /[\s\t]+/, $stLine;
		next if ($stInfo[2] ne "gene");
		$nTemp++;
		my $stID = $stInfo[$#stInfo];
		$stID =~ s/ID=//g;

		if((($charStrand{$stInfo[0]} =~ /[+-]/)&&($charStrand{$stInfo[0]} eq "$stInfo[6]"))||($charStrand{$stInfo[0]} !~ /[+-]/))
		{
			if($stFullLoci{$stInfo[0]} eq "")
			{
				push(@stData, "$stID	$stInfo[5]	$stInfo[3]	$stInfo[4]");
			}
			else
			{
				if($stType{$stID} eq "FULL")
				{
					push(@stData, "$stID	$stInfo[5]	$stInfo[3]	$stInfo[4]");
				}
			}
		}
		if($nLoci{$stInfo[0]} == $nTemp)
		{
			print FH5 "$nTemp	$stInfo[0]	$nLoci{$stInfo[0]}\n";
			my @stTemp = &SortAndExtractGeneModel (\@stData,$stFullLoci{$stInfo[0]});
			push(@stResult,@stTemp);

			$nTemp=0;
			@stData = ();
		}

	}
	close(DATA);

	for(my $i=0; $i<@stResult; $i++)
	{
		print FH6 "$stResult[$i]\n";
	}

	sub SortAndExtractGeneModel
	{
		my @stData = @{$_[0]};
		my $nFull = $_[1];
		my @stResult;
		my $nTemp = 0;


		@stData = sort
		{
				($b =~ /^[^\t]+[\t]+([0-9]+)[\t]+[0-9]+/)[0] <=> ($a =~ /^[^\t]+[\t]+([0-9]+)[\t]+[0-9]+/)[0]
		}@stData;

		if($nFull eq "")
		{
			my @stInfo = split /[\t]+/,$stData[0];
			if($stInfo[1]>=1000)
			{
				push(@stResult, "$stData[0]	Best_Partial");
			}
		}
		else
		{
			my $nStart = 0;
			my $nEnd = 0;
			my $nStart1=0;
			my $nEnd1=0;
			for(my $i=0; $i<3; $i++)
			{
				if($stData[$i] eq "")
				{
					last;
				}
				my @stInfo = split /[\t]+/, $stData[$i];

	#			next if ($stInfo[1]<=1000);

				if($i==0)
				{
					$nStart = $stInfo[2];
					$nEnd = $stInfo[3];
					push(@stResult, "$stData[$i]	Best_Full");
				}
				else
				{
					if(($nStart>$stInfo[3])||($nEnd<$stInfo[2]))
					{
						if ($nTemp != 1)
						{
							push(@stResult, "$stData[$i]	Representative");
							$nStart1 = $stInfo[2];
							$nEnd1 = $stInfo[3];
						}
						$nTemp = 1;
					}
					else
					{
						if(($nStart1<=$stInfo[3])&&($nStart1>=$stInfo[2]))
						{
							push(@stResult, "$stData[$i]	Repre_AS");
						}
						elsif(($nEnd1<=$stInfo[3])&&($nEnd1>=$stInfo[2]))
						{
							push(@stResult, "$stData[$i]	Repre_AS");
						}
						else
						{
							push(@stResult, "$stData[$i]	AS");
						}
					}
				}
			}
		}
		return (@stResult);
	}

close(FH6);

	my $stPrefix = $OUTPUT_PREFIX;
	my $stTranscriptInfo = "$stPrefix.transcripts.Info";
	my $stLociInfo = "$stPrefix.Loci.Info";
	my $stPEP = "$stPrefix.RefPEP.Consensus.PEP.fa";
	my $stGenome = "$stMaskedGenome";
	my $stOut = "$stPrefix.GeneModel";

	my ($stName,$stSeq,$nIdx);
	my %stList = {};
	my %stGenome = {};
	my %stType ={};
	my $stID;

	open(DATA, "$stGenome");
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		if ($stLine =~ /^>([^\s]+)/)
		{
			if ($stSeq ne "")
			{
				$stGenome{$stName} = $stSeq;
			}
			$stName = $1;
			$stSeq = "";
		}
		else
		{
			$stSeq = $stSeq.$stLine;
		}
	}
	$stGenome{$stName} = $stSeq;
	close(DATA);

	open(DATA, "$stLociInfo");
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		my @stInfo = split /[\t]+/, $stLine;
		my $stID = $stInfo[0];
		$stID =~ s/\..+//g;
		$stList{$stID}++;
		$stType{$stInfo[0]} = $stInfo[4];
		
	}
	close(DATA);

	open(DATA, "$stTranscriptInfo");
	open(OUT, ">$stOut.transcript.genomic.fa");
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		my @stInfo = split /[\t]+/, $stLine;
		my $stID = $stInfo[$#stInfo];
		$stID =~ s/ID=([0-9]+);.+/$1/g;
		next if ($stList{$stID} eq "");
		my $nLen = $stInfo[4] - $stInfo[3] + 1;

		$stSeq = substr($stGenome{$stInfo[0]},$stInfo[3]-1,$nLen);
		$stSeq = uc($stSeq);
		print OUT ">$stID\n$stSeq\n";
	}
	close(DATA);
	close(OUT);

	open(DATA, "$stPEP");
	open(BEST, ">$stOut.Repre.PEP.fa");
	open(OUT, ">$stOut.Whole.PEP.fa");
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		if ($stLine =~ />([^\t]+)/)
		{
			$stID = $1;
			$stID =~ s/\_[0-9]+//g;
			my $stLoci = $stID;
			$stLoci =~ s/\.[0-9]+//g;

			if(($stType{$stID} =~/Best/)||($stType{$stID} =~ /Representative/))
			{
				$nIdx = 1;
			}
			elsif($stType{$stID} =~ /AS/)
			{
				$nIdx = 2;
			}
			else
			{
				$nIdx = 0;
			}
		}
		else
		{
			if($nIdx == 1)
			{
				print BEST ">$stID"."_$stType{$stID}\n$stLine\n";
				print OUT ">$stID"."_$stType{$stID}\n$stLine\n";
			}
			elsif($nIdx == 2)
			{
				print OUT ">$stID"."_$stType{$stID}\n$stLine\n";
			}
		}
	}
	close(DATA);
	close(BEST);
	close(OUT);

	my $stPrefix = $OUTPUT_PREFIX;
	my $stData="$stPrefix.GeneModel.Whole.PEP.fa";

open (FH7, ">$stPrefix.GeneModel.list") || die ("Can't open myfile");
	open(DATA, "$stData");

	my ($stID,$stSeq)=("","");

	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		if($stLine =~ />([^\s]+)/)
		{
			if($stSeq ne "")
			{
				my $nLen=length($stSeq);
				print FH7 "  $stID             $nLen\n";
				$stSeq="";
			}
			$stID=$1;
		}
		else
		{
			$stSeq=$stSeq."$stLine";
		}
	}
	if($stSeq ne "")
	{
		my $nLen=length($stSeq);
		print FH7 "  $stID             $nLen\n";
		$stSeq="";
	}
	close DATA;
close(FH7);

	my $stPrefix = $OUTPUT_PREFIX;
	my $stTranscript = "$stPrefix.GeneModel.transcript.genomic.fa";
	my $stProtein = "$stPrefix.GeneModel.Whole.PEP.fa";
	my $nCPU = $THREADS;

	my @stTranscript = ();

	my %stInfo = {};
	my %stList = {};
	my %stProtein = {};
	my ($stName,$stID);
	my $nCnt = 0;
	open(DATA, "$stProtein");
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		if($stLine =~ /^>([^\s]+)/)
		{
			$stID = $stLine;
		}
		else
		{
			my $stLoci = $stID;
			$stLoci =~ s/>([0-9]+)\..+/$1/g;
			$stProtein{$stLoci} = $stProtein{$stLoci}."$stID\n$stLine\n";
		}
		
	}
	close(DATA);

	open(DATA, "$stTranscript");
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		if ($stLine =~ /^>/)
		{
			$stName = $stLine;
			$nCnt++;
		}
		else
		{
			open(TRANS, ">$stName.trans.seq");
			open(PROTEIN, ">$stName.protein.seq");

			print TRANS "$stName\n$stLine";
			$stName =~ s/>//g;
			print PROTEIN "$stProtein{$stName}";
			close(TRANS);
			close(PROTEIN);
			if($nCnt<$nCPU)
			{
				system("$EXONERATE_BIN_PATH/exonerate --model protein2genome --percent 70 --showtargetgff yes --showquerygff yes -t $stName.trans.seq -q $stName.protein.seq --querytype protein --targettype dna --maxintron $MAX_INTRON_LENGTH --forcegtag yes > $stName.proteinmapping.out &");
			#	print "$nCnt\n";
			}
			else
			{
				system("$EXONERATE_BIN_PATH/exonerate --model protein2genome --percent 70 --showtargetgff yes --showquerygff yes -t $stName.trans.seq -q $stName.protein.seq --querytype protein --targettype dna --maxintron $MAX_INTRON_LENGTH --forcegtag yes > $stName.proteinmapping.out");
				$nCnt = 0;
			}
		}
	}
	close(DATA);
sleep 3;
system("rm -rf NewAssembly.fa");

my $stData = "$OUTPUT_PREFIX.$REPRESENTATIVE_DOMAIN_NAME.ProteinInGenome.out"; ### *.ProteinInGenome.out
my $stMaskedGenome = "$RUNNING_PATH/$OUTPUT_PREFIX.tempID.Masked.Except$REPRESENTATIVE_DOMAIN_NAME.fasta";
my $stPrefix = $OUTPUT_PREFIX;

sleep 3;
system("cat *.proteinmapping.out > $stData");
system("rm -rf *.proteinmapping.out");
system("rm -rf *.seq");

	my $stData = "$OUTPUT_PREFIX.$REPRESENTATIVE_DOMAIN_NAME.ProteinInGenome.out"; #Target + Query
	my $stList = "$stPrefix.GeneModel.list";
	my $fCoverage = "99.5";
	my $stWeightName = "ISGAP";
	my $stQueryPrefix = "[0-9]+";
	my $stTranscriptInfo = "$stPrefix.transcripts.Info";
	my $stOut = "$stPrefix.Final.GeneModel.Whole.gff3";

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
	my %stScaffold = {};
	my %nStart = {};
	my $nTemp = 0;

	open(DATA, "$stTranscriptInfo");
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		my @stInfo = split /[\t]+/, $stLine;
		my $stID = $stInfo[$#stInfo];
		$stID =~ s/ID=([0-9]+);.+/$1/g;
		$stScaffold{$stID} = $stInfo[0];
		$nStart{$stID} = $stInfo[3];
	}
	close(DATA);

	open(DATA, "$stList");
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		my @stInfo = split /[\s\t]+/, $stLine;
		$stInfo[1] =~ s/Transcript_//g;
		$stInfo[1] =~ s/\//./g;
		$stInfo[1] =~ s/Confidence_[\.0-9]+_//g;
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
			if(($nMatchedLen/$stList{$stID})*100 >= $fCoverage)
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
				$stInfo[1] = $stWeightName;
				$stLine = join("	", @stInfo);
				if($stInfo[2] eq "gene")
				{
					my $stTemp = $stScaffold{$stInfo[0]};
					$stInfo[8] =~ s/gene_id ([0-9]+) .+/ID=$stID.$1;/g;
					$stAttribute = $stInfo[8];
					$nGeneID = $stAttribute;
					$nExtractGene++;
				}
				$stInfo[8] = $stAttribute;
				$stInfo[3] = $stInfo[3]+ $nStart{$stInfo[0]} -1;
				$stInfo[4] = $stInfo[4]+ $nStart{$stInfo[0]} -1;
				$stInfo[0] = $stScaffold{$stInfo[0]};
				$stLine = join("	", @stInfo);

				if($stLine =~ /similarity/)
				{
					$nTemp++;
					if($nTemp == 2)
					{
						print OUT "\n";
						$nTemp=0;
					}
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
	#	print "$nExtractGene\n";

	my $stPrefix = $OUTPUT_PREFIX;
	my $stGff = "$stPrefix.Final.GeneModel.Whole.gff3";

	my @stGene;
	my @stCDS;
	my $nIdx1=0;

open (FH, ">$stPrefix.Final.GeneModel.Whole.sorted.gff3") || die ("Can't open myfile");
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

	my $stPrefix = $OUTPUT_PREFIX;
	my $stScaffold = "$stMaskedGenome";
	my $stGff = "$stPrefix.Final.GeneModel.Whole.sorted.gff3";
	my $stOut = "$stPrefix.Final.GeneModel.Whole.CDS.fa";

	my $stName = "";
	my $stGCDS = "";
	my $charStrand = "";
	my %stScaffold = {};
	my $nCnt = 0;
	my $stSeq = "";
	my $nCnt = 0;


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
	my $nStr=0;
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
			$nStr=$stInfo[3];
		}

		next if ($stInfo[2] =~ /mRNA/);
		next if ($stInfo[2] =~ /exon/);
		next if ($stInfo[2] =~ /gene/);

		next if ($stName eq "");

		if($stInfo[5] ne "")
		{
			my $nExonLen = $stInfo[4] - $stInfo[3] + 1;
			my $nExonStart = $stInfo[3] - 1;
			my $stScaffold = $stScaffold{$stInfo[0]};

			my $stExon = substr($stScaffold,$nExonStart,$nExonLen);
			
	#		if($charStrand eq "-")
	#		{
	#			$stExon = reverse($stExon);
	#			$stExon =~ tr/ACGT/TGCA/;
	#		}
			$stGCDS = $stGCDS.$stExon;
		}
		else
		{
			if($charStrand eq "-")
			{
				$stGCDS = reverse($stGCDS);
				$stGCDS =~ tr/ACGT/TGCA/;
			}
			if(length($stGCDS)%3>0)
			{
				#				print "$stName\n";
			}
			if($stGCDS ne "")
			{
				print OUT ">$stName	$nStr\n$stGCDS\n";
			}
			$stGCDS = "";
			$stName = "";
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
		if(length($stGCDS)%3>0)
		{
			#			print "$stName\n";
		}
		print OUT ">$stName	$nStr\n$stGCDS\n";
	}

	#	print "$nCnt\n";

	my $stPrefix = $OUTPUT_PREFIX;
	my $stAssembly = "$stPrefix.Final.GeneModel.Whole.CDS.fa";
	my $stCodon = "$TGFAM_SCRIPTS_PATH/CodonUsage";
	my $stOut = "$stPrefix.Final.GeneModel.Whole";

	my $stName = "";
	my $stSeq = "";
	my %stCodon = {};
	my $nCut = 20;
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
					open (FH2, ">>Error.Final.list") || die ("Can't open myfile");
					print FH2 "$stName\n$stTrans1\n";
					close(FH2);
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

	my $stPrefix = $OUTPUT_PREFIX;
	my $stPEP = "$stPrefix.Final.GeneModel.Whole.PEP.fa";
	my $stCDS = "$stPrefix.Final.GeneModel.Whole.CDS.fa";
	my $stGff = "$stPrefix.Final.GeneModel.Whole.sorted.gff3";

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

#	print "$nCnt1   $nCnt2\n";


### Extraction of representative form in genome ###
sleep 3;
system("grep \"gene\" $stPrefix.Final.GeneModel.Whole.sorted.gff3.filter > $stPrefix.Final.GeneModel.Whole.Info");

	my $stPrefix = $OUTPUT_PREFIX;
	my $stFilterPEP="$stPrefix.Final.GeneModel.Whole.PEP.fa";
	my $stPositionInfo="$stPrefix.Final.GeneModel.Whole.Info";
	my $stOut="$stPrefix.Representative.Whole.Info";

	open(DATA, "$stFilterPEP");

	my %stList;
	my $nIdx1=0;
	my $nIdx2=0;

	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		if($stLine =~ /^>([^\s]+)	([^\s]+)	([^\s]+)/)
		{
			my $stID=$1;
			my $stPosition=$2;
			my $stType=$3;
			$stList{$stID."_$stPosition"} = "	$stType";
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
		if($stList{$stInfo[8]."_$stInfo[3]"} =~ /^	[^\s]+/)
		{
			if($stList{$stInfo[8]."_$stInfo[3]"} =~ /[0-9]+/)
			{
	#			print "$stLine\n";
			}
			$stList{$stInfo[8]."_$stInfo[3]"} = "$stInfo[5]".$stList{$stInfo[8]."_$stInfo[3]"};
			push(@stPos,"$stLine");
		}
	}
	close DATA;
	$nIdx1=@stPos;
	@stPos=sort{	($a =~ /^tempCh([0-9]+)	/)[0] <=> ($b =~ /^tempCh([0-9]+)	/)[0] ||
	#		($a =~ /cov[0-9]+\.([0-9]+)/)[0] <=> ($b =~ /cov[0-9]+\.([0-9]+)/)[0] ||
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
			my @stScoreType1=split(/\t/,$stList{$stInfo1[8]."_$stInfo1[3]"});
			my @stScoreType2=split(/\t/,$stList{$stInfo2[8]."_$stInfo2[3]"});
			if(($stScoreType1[1] =~ /Full/)&&($stScoreType2[1] !~ /Full/))
			{
				$stPos[$i+1]=$stPos[$i];
				if($i == $#stPos-1)
				{
					#$stPos[$i+1] =~ s/_[^\s]+$//g;
					$stPos[$i+1] = $stPos[$i+1]."$stScoreType1[1];";
					print FH "$stPos[$i+1]\n";
				}
			}
			elsif(($stScoreType1[1] !~ /Full/)&&($stScoreType2[1] =~ /Full/))
			{
				if($i == $#stPos-1)
				{
					#$stPos[$i+1] =~ s/_[^\s]+$//g;
					$stPos[$i+1] = $stPos[$i+1]."$stScoreType2[1];";
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
					$stPos[$i+1] = $stPos[$i+1]."$stScoreType2[1];";
					print FH "$stPos[$i+1]\n";
				}
			}
		}
		else
		{
			my @stScoreType1=split(/\t/,$stList{$stInfo1[8]."_$stInfo1[3]"});
			#$stPos[$i] =~ s/_[^\s]+$//g;
			$stPos[$i] = $stPos[$i]."$stScoreType1[1];";
			print FH "$stPos[$i]\n";
			if($i == $#stPos-1)
			{
				my @stScoreType2=split(/\t/,$stList{$stInfo2[8]."_$stInfo2[3]"});
				#$stPos[$i+1] =~ s/_[^\s]+$//g;
				$stPos[$i+1] = $stPos[$i+1]."$stScoreType2[1];";
				print FH "$stPos[$i+1]\n";
			}
		}
	}
	#	print "Overlap Transcripts: $nIdx2\n";

close(FH);
close(FH1);

	my $stPrefix = $OUTPUT_PREFIX;
	my $stData="$stPrefix.Representative.Whole.Info";



open (FH3, ">$stPrefix.Full_candidate.Whole.Info") || die ("Can't open myfile");
	open(DATA, "$stData");

	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		if(($stLine =~ /Best_Full/)&&($stLine =~ /Partial/))
		{
			print FH3 "$stLine\n";
		}
		elsif(($stLine =~ /Representative/)&&($stLine =~ /Partial/))
		{
			print FH3 "$stLine\n";
		}
		elsif(($stLine =~ /_AS/)&&($stLine =~ /Partial/))
		{
			print FH3 "$stLine\n";
		}
	}
	close DATA;
close(FH3);

	my $stPrefix = $OUTPUT_PREFIX;
	my $stGff1="$stPrefix.Final.GeneModel.Whole.sorted.gff3.filter";	###	Chi.Final.GeneModel.Whole.sort.gff3
	my $stInfoData="$stPrefix.Representative.Whole.Info";	###	tem.Whole

	my %stList;
	my %dupl;
	my $stKey="";
	my $num=0;

open (FH4, ">$stPrefix.Final.GeneModel.Whole.Repre.sort.gff3") || die ("Can't open myfile");
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
					$stList{$stKey}=$stList{$stKey}."$stLine\|";
					$num++;
				}
				$dupl{$stKey}="OK";
				
			}
			elsif($num == 1)
			{
				$stList{$stKey}=$stList{$stKey}."$stLine\|";
			}
		}
	}
	close DATA;

	my @stTotal=();

	open(DATA, "$stInfoData");

	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		push(@stTotal,"$stLine");
	}
	close DATA;

	@stTotal=sort{	($a =~ /^tempCh([0-9]+)	/)[0] <=> ($b =~ /^tempCh([0-9]+)	/)[0] ||
	#		($a =~ /cov[0-9]+\.([0-9]+)/)[0] <=> ($b =~ /cov[0-9]+\.([0-9]+)/)[0] ||
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
		$stGenePos =~ s/;[^\s]+;$/;/g;
		my @stInfo=split(/\|/,$stList{$stGenePos});
		for(my $j=0; $j<@stInfo; $j++)
		{
			print FH4 "$stInfo[$j]$stType\n";
		}
		print FH4 "\n";
	}
close(FH4);

	my $stPrefix = $OUTPUT_PREFIX;
	my $stScaffold = $stMaskedGenome;
	my $stGff = "$stPrefix.Final.GeneModel.Whole.Repre.sort.gff3";
	my $stOut = "$stPrefix.Final.GeneModel.Whole.Repre.sort.CDS.fa";

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
			#	print "Error	$stName\n";
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

	my $stPrefix = $OUTPUT_PREFIX;
	my $stAssembly = "$stPrefix.Final.GeneModel.Whole.Repre.sort.CDS.fa";
	my $stCodon = "$TGFAM_SCRIPTS_PATH/CodonUsage";
	my $stOut = "$stPrefix.Final.GeneModel.Whole.Repre.sort";

	my $stName = "";
	my $stSeq = "";
	my %stCodon = {};
	my $nCut = 20;
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
					open (FH5, ">>Error.repre.list") || die ("Can't open myfile");
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



#### non-redundant ###

	my $stPrefix = $OUTPUT_PREFIX;
	my $stData="$stPrefix.Final.GeneModel.Whole.PEP.fa";	###	Chi.Final.GeneModel.Whole.PEP.fa


open (FH6, ">$stPrefix.Final.GeneModel.Whole.NR.PEP.fa") || die ("Can't open myfile");
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

###	tsv file	###

	my $stPrefix = $OUTPUT_PREFIX;
	my $stData = "$stPrefix.Final.GeneModel.Whole.Repre.sort.PEP.fa";

open (FH7, ">$stPrefix.Final.GeneModel.Whole.Repre.sort.PEP.fa.noStop") || die ("Can't open myfile");
	open(DATA, "$stData");
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		if($stLine =~ /^>([^\s]+)/)
		{
			
			print FH7 ">$1\n";
		}
		else
		{
			$stLine =~ s/\*$//g;
			print FH7 "$stLine\n";
		}
	}
	close DATA;
close(FH7);

sleep 3;
system("rm -rf $stPrefix.Final.GeneModel.Whole.Repre.sort.PEP.fa.noStop.tsv");
system("$IPRSCAN_PATH/interproscan.sh -appl pfam -i $stPrefix.Final.GeneModel.Whole.Repre.sort.PEP.fa.noStop -f tsv");
system("rm -rf NewAssembly.fa");
