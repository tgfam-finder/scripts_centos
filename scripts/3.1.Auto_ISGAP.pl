use lib $ARGV[0];

BEGIN {
	require $ARGV[1];
		$ARGV[1]->import( qw(TGFAM_CONFIG) );
	require $ARGV[2];
		$ARGV[2]->import( qw(TGFAM_PATH) );
	}

use strict;

print "################# 3.1.Auto_ISGAP.pl is started #################\n";

system("$BOWTIE_PATH/bowtie2-build $RUNNING_PATH/$OUTPUT_PREFIX.tempID.Masked.Except$REPRESENTATIVE_DOMAIN_NAME.fasta.Main.fa $OUTPUT_PREFIX");

system("perl $TGFAM_SCRIPTS_PATH/TopHat_Cufflink_Pipeline.pl $ARGV[0] $ARGV[1] $ARGV[2] -num-threads $THREADS -Bowtie2_Index $OUTPUT_PREFIX -For $RNASEQ_FORWARD_PATH -Rev $RNASEQ_REVERSE_PATH");

sleep 3;
system("cp -rf $BLAST_DB_NAME* $RUNNING_PATH/$ISGAP_ANALYSIS_PATH");
sleep 3;
system("cp -rf $RESOURCE_PROTEIN $RUNNING_PATH/$ISGAP_ANALYSIS_PATH");


my $stRemoveRedundant = "$RUNNING_PATH/$OUTPUT_PREFIX.Final.$REPRESENTATIVE_DOMAIN_NAME\domain.search.out.Extended.RemoveRedundant"; ### *.NB-ARCdomain.search.out.Extended.RemoveRedundant
my $stGTF = "Cufflinks_0.out/transcripts.gtf"; ### transcripts.gtf
my $stMaskedGenome = "$RUNNING_PATH/$OUTPUT_PREFIX.tempID.Masked.Except$REPRESENTATIVE_DOMAIN_NAME.fasta.Main.fa";
my $outputPrefix = $OUTPUT_PREFIX; ### ATH, GMA.....
my $stSourcePEPindex = $BLAST_DB_NAME; ### 33species_P450_PEP.fasta,33species_NB-ARC_PEP.fasta... (Blast index file)
my $stSourcePEPfasta = $RESOURCE_PROTEIN; ### 33species_P450_PEP.fasta,33species_NB-ARC_PEP.fasta... (Fasta file)
my $nCPU = $THREADS;
my $stDomainName = $REPRESENTATIVE_DOMAIN_NAME;

	my $stData = "$RUNNING_PATH/$OUTPUT_PREFIX.Final.$REPRESENTATIVE_DOMAIN_NAME\domain.search.out.Extended.RemoveRedundant";
	my $stGTF = "Cufflinks_0.out/transcripts.gtf";
	my %IDPosi;
	my $str="";
	my $end="";
	my $chr="";
	my $num =0;

open (FH1, ">$outputPrefix.transcript.filtered.gft") || die ("Can't open myfile");
	open(DATA, "$stData");
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		my @stList = split /[\t]/, $stLine;

		$IDPosi{$stList[0]}=$IDPosi{$stList[0]}."$stList[1]	$stList[2],";	
	}
	close DATA;

	open(DATA, "$stGTF");
	open FH, ">", "$stGTF.Info.forValidation";
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		my @stList = split /[\t]/, $stLine;
		
		if($stList[2] =~ /transcript/)
		{
			$num=0;
			$str = $stList[3];
			$end = $stList[4];		
			$chr = $stList[0];
			
			$IDPosi{$chr} =~ s/,$//g;
			my @stPosi = split /,/, $IDPosi{$chr};
			for(my $i=0; $i<@stPosi; $i++)
			{
				my @StrEnd = split /[\t]/, $stPosi[$i];
				if(($str >= $StrEnd[0])&&($end <= $StrEnd[1]))
				{
					print FH1"$stLine\n";
					print FH "$chr	$StrEnd[0]	$StrEnd[1]	$stLine\n";
					$num++;
					next;
				}
			}
		}
		elsif($num == 1)
		{
			print FH1 "$stLine\n";
		}
	}
	close DATA;
	close FH;
close(FH1);

	my $stData = "$stGTF.Info.forValidation"; ### transcripts.gtf.Info.forValidation

	open(DATA, "$stData");
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		my @stList = split /[\t]/, $stLine;
		if(($stList[1] > $stList[6])||($stList[2] < $stList[7]))
		{
			#			print "error_gtf_filtering	$stLine\n";
		}
	}
	close DATA;

system("$TOPHAT_PATH/gtf_to_fasta $outputPrefix.transcript.filtered.gft $stMaskedGenome $outputPrefix.transcripts.fa");


	my $stDB=$BLAST_DB_NAME;
	my $nCPU=$THREADS;
	my $stInput="$outputPrefix.transcripts.fa";
	my $stOut="$outputPrefix.transcripts.vs.RefPEP.out";
#	my @stData=glob("$stInput.Number*");

	system("$BLAST_BIN_PATH/blastx -db $stDB -query $stInput -out $stOut -evalue 1e-5 -outfmt 7 -max_target_seqs 200 -num_threads $nCPU"); ### unefficient !!!

	my $stBlast = "$outputPrefix.transcripts.vs.RefPEP.out";
	my $stFasta = $RESOURCE_PROTEIN;

	my @stBlast = split /,/, $stBlast;
	my @stFasta = split /,/, $stFasta;

	my $stName;

	my %stName = {};

open (FH2, ">$outputPrefix.MatchedEvidenceProtein.list") || die ("Can't open myfile");
	for(my $i=0; $i<@stBlast; $i++)
	{
		open(DATA, "$stBlast[$i]");
		while(my $stLine = <DATA>)
		{
			chomp($stLine);
			next if($stLine =~ /^#/);
			my @stInfo = split /[\s\t]+/, $stLine;
			$stName{$stInfo[1]}++;
		}
		close(DATA);
	}


	for(my $i=0; $i<@stFasta; $i++)
	{
		my $stSeq;
		open(DATA, "$stFasta[$i]");
		while(my $stLine = <DATA>)
		{
			chomp($stLine);
			if($stLine =~ /^>([^\s]+)/)
			{
				if ($stSeq ne "")
				{
					print FH2" $stName	".length($stSeq)."\n";
				}
				$stName = $1;
				$stSeq = "";
			}
			else
			{
				$stSeq = $stSeq.$stLine;
			}
		}
		print FH2" $stName	".length($stSeq)."\n";
		close(DATA);
	}

close(FH2);

unless (-d "temp_ISGAP.$outputPrefix\_$stDomainName")
{
	system("mkdir temp_ISGAP.$outputPrefix\_$stDomainName");
}

	my $stBlast = "$outputPrefix.transcripts.vs.RefPEP.out";
	my $stPEP = $RESOURCE_PROTEIN;
	my $stGenome = "$outputPrefix.transcripts.fa";

	my %stGenome = {};
	my %stProtein = {};
	my %stList = {};
	my %stProteinList = {};
	my ($stName,$stSeq,$stProtein) = "";
	my $nIdx = 0;
	my $nCnt = 0;
	my $nMatchedSeq = 0;
	my @stGenome = ();	###	MS_Modification
	my @stProtein = ();	###	MS_Modification

	open(DATA, "$stGenome");
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		if($stLine =~ /^>([^\s]+)/)
		{
			if($stSeq ne "")
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

	open(DATA, "$stPEP");
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		if($stLine =~ /^>([^\s]+)/)
		{
			if($stProtein ne "")
			{
				$stProtein{$stName} = $stProtein;
			}
			$stName = $1;
			$stProtein = "";
		}
		else
		{
			$stProtein = $stProtein.$stLine;
		}
	}
	close(DATA);
	$stProtein{$stName} = $stProtein;

	open(DATA, "$stBlast");
	while (my $stLine = <DATA>)
	{
		chomp($stLine);
		if($stLine =~ /^# Query: ([^\s]+)/)
		{
			$stName = $1;
			$nCnt = 0;
		}
		elsif($stLine =~ /^# ([0-9]+) hits found/)
		{
			if($1>0)
			{
				$nIdx = 1;
				$nMatchedSeq++;
			}
			else
			{
				$nIdx = 0;
			}
		}
		if($nIdx == 1)
		{
			if($stLine =~ /^[^#]/)
			{
				$nCnt++;
				my @stInfo = split /[\s\t]+/, $stLine;
				$stInfo[1] =~ s/#.+//g;
				if($nCnt == 1)
				{
					open(OUT1, ">$stName.genome.fa");
					open(OUT2, ">$stName.protein.fa");
					push(@stGenome, "$stName.genome.fa");		###		MS_Modification
					push(@stProtein, "$stName.protein.fa");		###		MS_Modification
					print OUT1 ">$stName\n$stGenome{$stName}\n";
				}
				if($stProteinList{$stInfo[1]}==0)
				{
					if($stProtein{$stInfo[1]} !~ /[Uu]/)
					{
						next if ($stProteinList{$stProtein{$stInfo[1]}}>1);
						print OUT2 ">$stInfo[1]\n$stProtein{$stInfo[1]}\n";
					}
				}
				$stProteinList{$stInfo[1]}++;
				$stProteinList{$stProtein{$stInfo[1]}}++;
			}
			else
			{
				$nCnt = 0;
				%stProteinList = {};
				close(OUT1);
				close(OUT2);
			}
		}
	}
	close(DATA);

	#print "$#stGenome will be aligned using proteins\n";

sleep 3;
system("mv *.genome.fa temp_ISGAP.$outputPrefix\_$stDomainName");
sleep 3;
system("mv *.protein.fa temp_ISGAP.$outputPrefix\_$stDomainName");

print "\n################# 3.1.Auto_ISGAP.pl is finished ################\n\n";
