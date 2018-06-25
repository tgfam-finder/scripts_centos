use lib $ARGV[0];

BEGIN {
	require $ARGV[1];
		$ARGV[1]->import( qw(TGFAM_CONFIG) );
	require $ARGV[2];
		$ARGV[2]->import( qw(TGFAM_PATH) );
	}

use strict;

my $stData = "$RUNNING_PATH/$OUTPUT_PREFIX.tempID.Masked.Except$REPRESENTATIVE_DOMAIN_NAME.fasta.Main.fa"; ### masking genome
my $stGene = "$RUNNING_PATH/$MERGING_ANALYSIS_PATH/$OUTPUT_PREFIX.TrainingSet.ForAugustus.NR.fasta"; ### Prefix.TrainingSet.ForAugustus.NR.fasta
my $stPrefix = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME"; ### Prefix_gene ex) ATH_NLR, ATH_PPR, ATH_MADS


###     Augustus training   ###
#######################################################
my $stGenome = "$RUNNING_PATH/$PROTEIN_MAPPING_ANALYSIS_PATH/$OUTPUT_PREFIX.Genome.tempID.Masked.Except$REPRESENTATIVE_DOMAIN_NAME.fasta.validation.fa";
#######################################################
my $stFULLPEP = "$RUNNING_PATH/$MERGING_ANALYSIS_PATH/$OUTPUT_PREFIX.TrainingSet.ForAugustus.NR.fasta";
my $stAugustusPATH = "$AUGUSTUS_PATH";
my $stScipioPATH = "$SCIPIO_PATH";
my $stBlatPATH = "$BLAT_PATH";
my $stSpecies = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME";
my $stOut = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME";

system("rm -rf $RUNNING_PATH/$AUGUSTUS_ANALYSIS_PATH/*");
system("rm -rf $stAugustusPATH/config/species/$stSpecies");
#######################################################
system ("perl $stScipioPATH/scipio.1.4.1.pl --blat_output=$stOut.prot.vs.genome.psl --blat_bin=$BLAT_PATH/blat $stGenome $stFULLPEP > $stOut.yaml");
#######################################################
system ("cat $stOut.yaml | $stScipioPATH/yaml2gff.1.4.pl > $stOut.scipiogff");
system("$stAugustusPATH/scripts/scipiogff2gff.pl --in=$stOut.scipiogff --out=$stOut.gff");
system("cat $stOut.yaml | $stScipioPATH/yaml2log.1.4.pl > $stOut.log");
system("$stAugustusPATH/scripts/gff2gbSmallDNA.pl $stOut.gff $stGenome 1000 $stOut.raw.gb");
system("$stAugustusPATH/scripts/new_species.pl --species=$stSpecies --AUGUSTUS_CONFIG_PATH=$stAugustusPATH/config");
system("$stAugustusPATH/bin/etraining --species=$stSpecies --stopCodonExcludedFromCDS=true $stOut.raw.gb 2> train.err");
system("ls -ort $stAugustusPATH/config/species/$stSpecies/");
system("cat train.err | perl -pe 's/.*in sequence (\\S+): .*/\$1/' > badgenes.lst\n"); # badgene
system("perl $stAugustusPATH/scripts/filterGenes.pl badgenes.lst $stOut.raw.gb > $stOut.gb");
system("$stAugustusPATH/bin/etraining --species=$stSpecies --stopCodonExcludedFromCDS=true $stOut.gb");

###	Augustus	###
#######################################################
system("$AUGUSTUS_BIN_PATH/augustus --AUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH --species=$stPrefix --genemodel=partial --gff3=on --noInFrameStop=true $stGenome > $stPrefix.GeneInfo"); 
#######################################################

###	gff, CDS, PEP	###

	my $stData = "$stPrefix.GeneInfo";
	my @stData = ();
	my @stTemp = ();

	open (FH1, ">$stPrefix.augustus.gff3") || die ("Can't open myfile");
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
				#######################################################
				my @stTempInfo = split /[\t]+/, $stData[$i];
				my @stPosInfo = split /[\_]/, $stTempInfo[0];
				$stTempInfo[0] = $stPosInfo[0];
				$stTempInfo[3] = $stPosInfo[1] + $stTempInfo[3] - 1;
				$stTempInfo[4] = $stPosInfo[1] + $stTempInfo[4] - 1;
				$stData[$i] = join("	", @stTempInfo);
				#######################################################
				print FH1"$stData[$i]\n";
			}
			print FH1"\n";
			@stData = ();
		}

		next if ($stLine =~ /^#/);
		next if ($stLine =~ /start_codon/);
		next if ($stLine =~ /intron/);
		push(@stData, "$stLine");
	}
	close(DATA);
	close(FH1);

	my $stScaffold = "$RUNNING_PATH/$OUTPUT_PREFIX.tempID.Masked.Except$REPRESENTATIVE_DOMAIN_NAME.fasta.Main.fa";
	my $stGff = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.augustus.gff3";
	my $stOut = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.augustus.CDS.fa";

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
		my @stPosInfo = split /\_/, $stInfo[0];

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
	close(OUT);

	my $stAssembly = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.augustus.CDS.fa";
	my $stCodon = "$TGFAM_SCRIPTS_PATH/CodonUsage";
	my $stOut = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.augustus";

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

	open (FH2, ">temp_$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME") || die ("Can't open myfile");
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

	my $stPEP = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.augustus.PEP.fa";
	my $stCDS = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.augustus.CDS.fa";
	my $stGff = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.augustus.gff3";

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

###	tsv file	###
	my $stData = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.augustus.PEP.fa.filter";

	open (FH3, ">$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.augustus.PEP.fa.filter.noStop") || die ("Can't open myfile");
	open(DATA, "$stData");
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		if($stLine =~ /^>([^\s]+)/)
		{
			
			print FH3 ">$1\n";
		}
		else
		{
			$stLine =~ s/\*$//g;
			print FH3 "$stLine\n";
		}
	}
	close DATA;
	close(FH3);

sleep 3;
system("rm -rf  $stPrefix.augustus.PEP.fa.filter.noStop.tsv");
system("$IPRSCAN_PATH/interproscan.sh -appl pfam -i $stPrefix.augustus.PEP.fa.filter.noStop -f tsv");
system("cp *.filter *.tsv $RUNNING_PATH/$MERGING_ANALYSIS_PATH");
system("rm -rf NewAssembly.fa");
