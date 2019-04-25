use lib $ARGV[0];

BEGIN {
	require $ARGV[1];
		$ARGV[1]->import( qw(TGFAM_CONFIG) );
	require $ARGV[2];
		$ARGV[2]->import( qw(TGFAM_PATH) );
	}

use strict;

print "############### 4.3.Auto_Augustus.pl is started ################\n";

system("cp -rf $RUNNING_PATH/$ISGAP_ANALYSIS_PATH/temp_ISGAP.$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME/$OUTPUT_PREFIX.transcripts.Info $OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.transcripts.Info"); 


my $stData = "$RUNNING_PATH/$ISGAP_ANALYSIS_PATH\/$OUTPUT_PREFIX.transcripts.fa"; ### transcript
my $stPrefix = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME"; ### Prefix_gene ex) ATH_NLR, ATH_PPR, ATH_MADS

###	Augustus	###
system("$AUGUSTUS_BIN_PATH/augustus --species=$stPrefix --genemodel=partial --gff3=on --noInFrameStop=true $stData > $stPrefix.GeneInfo"); 

###	gff, CDS, PEP	###

	my $stData = "$stPrefix.GeneInfo";
	my @stData = ();
	my @stTemp = ();


open (FH, ">$stPrefix.ISGAP.augustus.gff3") || die ("Can't open myfile");
	open(DATA, "$stData");
	while(my $stLine = <DATA>)
	{
		chomp($stLine);

		if(($stData[0] ne "")&&($stLine =~ /^#/))
		{
			for(my $i=0; $i<@stData; $i++)
			{
				my @stInfo = split /[\s\t]+/,$stData[$i];
				if($stInfo[2] eq "stop_codon")
				{
					if($stInfo[6] eq "-")
					{
						my @stTemp = split /[\s\t]+/, $stData[$i+1];
						$stTemp[3] = $stInfo[3];
						$stData[$i+1] = join('	',@stTemp);
					}
					else
					{
						my @stTemp = split /[\s\t]+/, $stData[$i-1];
						$stTemp[4] = $stInfo[4];
						$stData[$i-1] = join('	',@stTemp);
					}
				}
				
			}
			for(my $i=0; $i<@stData; $i++)
			{
				next if ($stData[$i] =~ /stop_codon/);
				$stData[$i] =~ s/transcript/mRNA/g;
				print FH "$stData[$i]\n";
			}
			print FH "\n";
			@stData = ();
		}

		next if ($stLine =~ /^#/);
		next if ($stLine =~ /start_codon/);
		next if ($stLine =~ /intron/);
		push(@stData, "$stLine");

		


	}
	close(DATA);
close(FH);

	my $stScaffold = "$RUNNING_PATH/$ISGAP_ANALYSIS_PATH\/$OUTPUT_PREFIX.transcripts.fa";
	my $stGff = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.ISGAP.augustus.gff3";
	my $stOut = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.ISGAP.augustus.CDS.fa";

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
				#		print "Error	$stName\n";
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

	my $stAssembly = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.ISGAP.augustus.CDS.fa";
	my $stCodon = "$TGFAM_SCRIPTS_PATH/CodonUsage";
	my $stOut = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.ISGAP.augustus";

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

	open (FH2, ">Error.temp_$stPrefix.ISGAP.augustus") || die ("Can't open myfile");
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
					print FH2 "$stName\n$stTrans1\n";
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
	close(FH2);

	my $stPEP = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.ISGAP.augustus.PEP.fa";
	my $stCDS = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.ISGAP.augustus.CDS.fa";
	my $stGff = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.ISGAP.augustus.gff3";

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
system("rm -rf NewAssembly.fa");
	#	print "$nCnt1   $nCnt2\n";


my $stFinalPEP = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.ISGAP.augustus.PEP.fa.filter"; ### *_MADS.PEP.fa
my $stFinalGff = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.ISGAP.augustus.gff3.filter"; ### *_MADS.augustus.gff3.filter
my $stData2 = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.AugustusInGenome.out"; ### *.AugustusInGenome.out
my $stMaskedGenome = "$RUNNING_PATH/$OUTPUT_PREFIX.tempID.Masked.Except$REPRESENTATIVE_DOMAIN_NAME.fasta.Main.fa";

	my $stData2="$stFinalPEP";

	open (FH, ">$stPrefix.Augustus.list") || die ("Can't open myfile");
	open(DATA, "$stData2");

	my ($stID,$stSeq)=("","");

	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		if($stLine =~ />([^\s]+)/)
		{
			if($stSeq ne "")
			{
				my $nLen=length($stSeq);
				print FH "  $stID             $nLen\n";
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
		print FH "  $stID             $nLen\n";
		$stSeq="";
	}
	close DATA;
close(FH);

	my $stGff = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.ISGAP.augustus.gff3.filter"; ### *.augustus.gff3.filter
	my $stGenome = "$RUNNING_PATH/$OUTPUT_PREFIX.tempID.Masked.Except$REPRESENTATIVE_DOMAIN_NAME.fasta.Main.fa"; ### masked genome
	my $stTranscriptInfo = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.transcripts.Info"; ### *.transcripts.Info
	my $stOut = $stTranscriptInfo;
	$stOut =~ s/\.transcripts\.Info//;
	my %stList;
	my %stGenome;
	my $stName="";
	my $stSeq = "";

	open(DATA, "$stGff");
	while(my $stLine = <DATA>)
	{
		my @stInfo = split /[\t]/, $stLine;
		if($stInfo[2] eq "gene")
		{
			$stList{$stInfo[0]}++;
		}	
	}
	close DATA;

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

	open(OUT, ">$stOut.transcript.genomic.fa");
	open(DATA, "$stTranscriptInfo");
	while(my $stLine = <DATA>)
	{
			chomp($stLine);
			my @stInfo = split /[\t]+/, $stLine;
			my $stID = $stInfo[$#stInfo];
			$stID =~ s/ID=([0-9]+);.+/$1/g;
			next if ($stList{$stID} eq "");
			my $nLen = $stInfo[4] - $stInfo[3] + 1;

			$stSeq = substr($stGenome{$stInfo[0]},$stInfo[3]-1,$nLen);
			print OUT ">$stID\n$stSeq\n";
	}
	close(DATA);
	close(OUT);

	my $stGff = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.ISGAP.augustus.gff3.filter"; ### *.augustus.gff3.filter
	my $stTranscript = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.transcript.genomic.fa"; ### *.transcript.genomic.fa
	my $stProtein = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.ISGAP.augustus.PEP.fa.filter"; ### *_MADS.PEP.fa
	my $nCPU = 20;

	my @stTranscript = ();

	my %stInfo = {};
	my %stList = {};
	my %stProtein = {};
	my %nLoci = {};
	my ($stName,$stID);
	my $nCnt = 0;

	open(DATA, "$stGff");
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		my @stInfo = split /[\t]/, $stLine;
		$stInfo[8] =~ s/ID=//g;
		$nLoci{$stInfo[8]} = "$stInfo[0]";
	}
	close DATA;

	open(DATA, "$stProtein");
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		if($stLine =~ /^>([^\s]+)/)
		{
			$stID = $1;
		}
		else
		{
			my $stLoci = $nLoci{$stID};
			$stProtein{$stLoci} = $stProtein{$stLoci}.">$stID\n$stLine\n";
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
				system("$EXONERATE_BIN_PATH/exonerate --model protein2genome --percent 70 --showtargetgff yes --showquerygff yes -t $stName.trans.seq -q $stName.protein.seq --querytype protein --targettype dna --maxintron 30000 --forcegtag yes > $stName.proteinmapping.out &");
				#				print "$nCnt\n";
			}
			else
			{
				system("$EXONERATE_BIN_PATH/exonerate --model protein2genome --percent 70 --showtargetgff yes --showquerygff yes -t $stName.trans.seq -q $stName.protein.seq --querytype protein --targettype dna --maxintron 30000 --forcegtag yes > $stName.proteinmapping.out &");
				$nCnt = 0;
			}
		}
	}
	close(DATA);

system("cat *.proteinmapping.out > $OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.AugustusInGenome.out");
system("rm -rf *.proteinmapping.out");
system("rm -rf *.seq");

	my $stData2 = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.AugustusInGenome.out"; #Target + Query
	my $stList = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.Augustus.list";
	my $fCoverage = "99.5";
	my $stWeightName = "ISGAPaugu";
	my $stQueryPrefix = "[0-9]+";
	my $stTranscriptInfo = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.transcripts.Info";
	my $stOut = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.ISGAP.Augustus.Whole.gff3";

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
		$stList{$stInfo[1]} = $stInfo[$#stInfo];
		
	}
	close(DATA);

	open(DATA, "$stData2");
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
					if($nTemp == 1)
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

	my $stGff = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.ISGAP.Augustus.Whole.gff3";

	my @stGene;
	my @stCDS;
	my $nIdx1=0;

open (FH2, ">$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.ISGAP.Augustus.Whole.sorted.gff3") || die ("Can't open myfile");
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
				print FH2 "$stGene[$i]\n";
			}
			for(my $i=0; $i<@stCDS; $i++)
			{
				my @stTemp = split /[\s\t]+/, $stCDS[$i];
				$stTemp[2] = "exon";	
				my $stExon = join ("	", @stTemp);
				print FH2 "$stExon\n";
				print FH2 "$stCDS[$i]\n";
			}
			print FH2 "$stLine\n";
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
				print FH2 "$stGene[$i]\n";
			}
			for(my $i=0; $i<@stCDS; $i++)
			{
				my @stTemp = split /[\s\t]+/, $stCDS[$i];
				$stTemp[2] = "exon";	
				my $stExon = join ("	", @stTemp);
				print FH2 "$stExon\n";
				print FH2 "$stCDS[$i]\n";
			}
			print FH2 "\n";
			@stCDS = ();
			@stGene = ();
	}

close(FH2);

	my $stScaffold = "$RUNNING_PATH/$OUTPUT_PREFIX.tempID.Masked.Except$REPRESENTATIVE_DOMAIN_NAME.fasta.Main.fa";
	my $stGff = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.ISGAP.Augustus.Whole.sorted.gff3";
	my $stOut = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.ISGAP.Augustus.Whole.CDS.fa";

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

	my $stAssembly = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.ISGAP.Augustus.Whole.CDS.fa";
	my $stCodon = "$TGFAM_SCRIPTS_PATH/CodonUsage";
	my $stOut = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.ISGAP.Augustus.Whole";

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


open (FH3, ">Error.ISGAP.Augustus.list") || die ("Can't open myfile");
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
					print FH3 "$stName\n$stTrans1\n";
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
close(FH3);

	my $stPEP = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.ISGAP.Augustus.Whole.PEP.fa";
	my $stCDS = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.ISGAP.Augustus.Whole.CDS.fa";
	my $stGff = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.ISGAP.Augustus.Whole.sorted.gff3";

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
system("grep \"gene\" $stPrefix.ISGAP.Augustus.Whole.sorted.gff3.filter > $stPrefix.ISGAP.Augustus.Whole.Info");

	my $stFilterPEP="$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.ISGAP.Augustus.Whole.PEP.fa";
	my $stPositionInfo="$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.ISGAP.Augustus.Whole.Info";
	my $stOut="$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.ISGAP.Augustus.Representative.Whole.Info";

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
				#				print "$stLine\n";
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

	my $stGff1="$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.ISGAP.Augustus.Whole.sorted.gff3.filter";	###	Chi.Final.GeneModel.Whole.sort.gff3
	my $stInfoData="$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.ISGAP.Augustus.Representative.Whole.Info";	###	tem.Whole

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

open (FH4, ">$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.ISGAP.Augustus.Whole.Repre.sort.gff3") || die ("Can't open myfile");
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

	my $stScaffold = "$RUNNING_PATH/$OUTPUT_PREFIX.tempID.Masked.Except$REPRESENTATIVE_DOMAIN_NAME.fasta.Main.fa";
	my $stGff = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.ISGAP.Augustus.Whole.Repre.sort.gff3";
	my $stOut = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.ISGAP.Augustus.Whole.Repre.sort.CDS.fa";

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
				#				print "Error	$stName\n";
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

	close(OUT);

	my $stAssembly = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.ISGAP.Augustus.Whole.Repre.sort.CDS.fa";
	my $stCodon = "$TGFAM_SCRIPTS_PATH/CodonUsage";
	my $stOut = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.ISGAP.Augustus.Whole.Repre.sort";

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


open (FH5, ">Error.repre.list") || die ("Can't open myfile");
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
					print FH5 "$stName\n$stTrans1\n";
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
close(FH5);

	my $stData2 = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.ISGAP.Augustus.Whole.Repre.sort.PEP.fa";

open (FH6, ">$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.ISGAP.Augustus.Whole.Repre.sort.PEP.fa.noStop") || die ("Can't open myfile");
	open(DATA, "$stData2");
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		if($stLine =~ /^>([^\s]+)/)
		{
			
			print FH6 ">$1\n";
		}
		else
		{
			$stLine =~ s/\*$//g;
			print FH6 "$stLine\n";
		}
	}
	close DATA;
close(FH6);

sleep 3;
system("rm -rf $stPrefix.ISGAP.Augustus.Whole.Repre.sort.PEP.fa.noStop.tsv");
system("$IPRSCAN_PATH/interproscan.sh -appl pfam -i $stPrefix.ISGAP.Augustus.Whole.Repre.sort.PEP.fa.noStop -f tsv");

### system("perl /var2/TGFam/scripts/RunHmm.Generate.tsv.pl $HMMER_BIN_PATH /var2/script/PF00096.hmm $stPrefix.ISGAP.Augustus.Whole.Repre.sort.PEP.fa.noStop $stPrefix.ISGAP.Augustus.Whole.Repre.sort.PEP.fa.noStop.tsv 08-02-2019 $stPrefix.ISGAP.Augustus.Whole.Repre.sort");

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
my $now = sprintf("%04d-%02d-%02d", $year + 1900, $mon + 1, $mday);

my $stDate = $now;
my ($stDomain,$nIdx,$stPFam,$stGene, $nDomainLen);

if ($HMM_MATRIX_NAME eq "") {
}
else {

system("$HMMER_BIN_PATH/hmmsearch $HMM_MATRIX_NAME $OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.ISGAP.Augustus.Whole.Repre.sort.PEP.fa.noStop > $OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.ISGAP.Augustus.Whole.Repre.sort.search.out");

open(DATA, "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.ISGAP.Augustus.Whole.Repre.sort.search.out");
open(OUT, ">$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.ISGAP.Augustus.Whole.Repre.sort.SelfBuild.Hmm.out");
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

system("cat $OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.ISGAP.Augustus.Whole.Repre.sort.SelfBuild.Hmm.out >> $OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.ISGAP.Augustus.Whole.Repre.sort.PEP.fa.noStop.tsv");

}

system("cp $stPrefix.ISGAP.Augustus.Whole.Repre.sort.*.fa *.tsv $stPrefix.ISGAP.Augustus.Whole.Repre.sort.gff3 $RUNNING_PATH/$MERGING_ANALYSIS_PATH");

system("cp *.ISGAP.Augustus.Whole.Repre.sort.* $RUNNING_PATH/$MERGING_ANALYSIS_PATH");
system("rm -rf NewAssembly.fa");

print "\n############### 4.3.Auto_Augustus.pl is finished ###############\n\n";
