use lib $ARGV[0];

BEGIN {
	require $ARGV[1];
		$ARGV[1]->import( qw(TGFAM_CONFIG) );
	require $ARGV[2];
		$ARGV[2]->import( qw(TGFAM_PATH) );
	}

use strict;

print "############ 2.1.Auto_ProteinMapping.pl is started #############\n";

system("cp -rf $BLAST_DB_NAME* $RUNNING_PATH/$PROTEIN_MAPPING_ANALYSIS_PATH");
system("cp -rf $RESOURCE_PROTEIN $RUNNING_PATH/$PROTEIN_MAPPING_ANALYSIS_PATH;");

my $stDomainName = $REPRESENTATIVE_DOMAIN_NAME; ### NB-ARC, SRF, AP2, PPR
my $outputPrefix = $OUTPUT_PREFIX;
my $stSourcePEPindex = $BLAST_DB_NAME; ### 33species_P450_PEP.fasta,33species_NB-ARC_PEP.fasta... (Blast index file)
my $stSourcePEPfasta = $RESOURCE_PROTEIN; ### 33species_P450_PEP.fasta,33species_NB-ARC_PEP.fasta... (Fasta file)
my $nCPU = $THREADS;

	my $stData = "$RUNNING_PATH/$outputPrefix.tempID.Masked.Except$stDomainName.fasta";
	my $stRemoveRedundant = "$RUNNING_PATH/$outputPrefix.tempID.Masked.Except$stDomainName.fasta.validation";
	my %IDSeq;
	my $stID="";

open (FH, ">$outputPrefix.Genome.tempID.Masked.Except$stDomainName.fasta.validation.fa") || die ("Can't open myfile");
	open(DATA, "$stData");
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		if($stLine =~ /^>([^\s]+)/)	
		{
			$stID = $1;
		}		
		else
		{	
			$IDSeq{$stID}="$stLine";
		}
	}
	close DATA;

	open(DATA, "$stRemoveRedundant");
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		my @stList = split /[\t]/, $stLine;
		$stLine =~ s/[\t]/\_/g;
		my $stSeq = $IDSeq{$stList[0]};
		my $length = $stList[2] - $stList[1] +1;
		my $start = $stList[1] - 1;
		my $fasta = substr($stSeq, $start, $length);
		print FH ">$stLine\n$fasta\n";
	}
	close DATA;
close(FH);

	my $stDB=$BLAST_DB_NAME;
	my $nCPU=$THREADS;
	my $stInput="$outputPrefix.Genome.tempID.Masked.Except$stDomainName.fasta.validation.fa";
	my $stOut="$outputPrefix.genome.vs.RefPEP.out";

	system("$BLAST_BIN_PATH/blastx -db $stDB -query $stInput -out $stOut -evalue 1e-5 -outfmt 7 -max_target_seqs 200 -num_threads $nCPU"); ### unefficient !!!

open (FH1, ">$outputPrefix.$stDomainName.MatchedEvidenceProtein.list") || die ("Can't open myfile");
	my $stBlast = "$outputPrefix.genome.vs.RefPEP.out";
	my $stFasta = $RESOURCE_PROTEIN;

	my @stBlast = split /,/, $stBlast;
	my @stFasta = split /,/, $stFasta;

	my $stName;

	my %stName = {};

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
					print FH1 " $stName	".length($stSeq)."\n";
				}
				$stName = $1;
				$stSeq = "";
			}
			else
			{
				$stSeq = $stSeq.$stLine;
			}
		}
		print FH1 " $stName	".length($stSeq)."\n";
		close(DATA);
	}
close(FH1);

unless(-d "temp_protein.$outputPrefix\_$stDomainName")
{
	system("mkdir temp_protein.$outputPrefix\_$stDomainName");
}

my $stBlast = "$OUTPUT_PREFIX.genome.vs.RefPEP.out";
my $stPEP = $RESOURCE_PROTEIN;
my $stGenome = "$OUTPUT_PREFIX.Genome.tempID.Masked.Except$stDomainName.fasta.validation.fa";

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
system("mv *.genome.fa temp_protein.$outputPrefix\_$stDomainName");
sleep 3;
system("mv *.protein.fa temp_protein.$outputPrefix\_$stDomainName");

print "\n############ 2.1.Auto_ProteinMapping.pl is finished ############\n\n";
