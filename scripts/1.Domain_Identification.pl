use lib $ARGV[0];

BEGIN {
	require $ARGV[1];
		$ARGV[1]->import( qw(TGFAM_CONFIG) );
	require $ARGV[2];
		$ARGV[2]->import( qw(TGFAM_PATH) );
	}

use strict;

print "############ 1.Domain_Identification.pl is started #############\n";

my $stData = $TARGET_GENOME; ### genome fasta file (tempID)
my $stGff3 = $GFF3_OF_TARGET_GENOME;
my $stTSV = "$RUNNING_PATH/$OUTPUT_PREFIX.tsv.addition";
my $stTSVOutput = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME";
my $stPEP = $PROTEINS_FOR_DOMAIN_IDENTIFICATION;
my $stPrefix = $OUTPUT_PREFIX;
my $stPfamIDs = $TARGET_DOMAIN_ID;
my $stDomains = $TARGET_DOMAIN_NAME;
my $stRepreDom = $REPRESENTATIVE_DOMAIN_NAME;
my $nPerc = "100";	###	percent

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
my $now = sprintf("%04d-%02d-%02d", $year + 1900, $mon + 1, $mday); 

my $stDate = $now;
my ($stDomain,$nIdx,$stPFam,$stGene, $nDomainLen);

my @stPfamID=split(/,/,$stPfamIDs);
my @stDomain=split(/,/,$stDomains);

my $stOut = $stData;
$stOut =~ s/\.fna//;
$stOut =~ s/\.fasta//;
$stOut =~ s/\.fa//;
my $stMerge = "";

if ($HMM_MATRIX_NAME eq "") {
}
else {

system("$HMMER_BIN_PATH/hmmsearch $HMM_MATRIX_NAME $PROTEINS_FOR_DOMAIN_IDENTIFICATION > $RUNNING_PATH/$OUTPUT_PREFIX.tsv.addition.search.out");

open(DATA, "$RUNNING_PATH/$OUTPUT_PREFIX.tsv.addition.search.out");
open(OUT, ">$RUNNING_PATH/$OUTPUT_PREFIX.addition.SelfBuild.Hmm.out");
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

}

system("cat $TSV_FOR_DOMAIN_IDENTIFICATION $OUTPUT_PREFIX.addition.SelfBuild.Hmm.out > $RUNNING_PATH/$OUTPUT_PREFIX.tsv.addition");

for(my $i=0; $i<@stPfamID; $i++)
{
###	Make HMM training set
	
	my $stData = $stTSV;
	my $cut = $nPerc;
	my $stPfam = $stPfamID[$i];
	my %NB;
	my $num=0;

open FH, ">", "$RUNNING_PATH/$stTSVOutput.FilteredBy$stDomain[$i]domain.$nPerc" || die ("Can't open myfile");
	open(DATA, "$stData");
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		my @stList = split /[\t]/, $stLine;
		if(($stList[0] ne "")&&($stList[4] eq "$stPfam"))	###	PF00931=NB-ARC, PF00319=MADS
		{
			$NB{$stLine} = "OK"; 	
		}
	}
	close DATA;

	my @nTotal=keys(%NB);
	my $nPerc=int(($#nTotal+1)*$cut/100);

	foreach my $nkey (keys %NB) ### 30 or @sort
	{
		$num++;
		if($num <= $nPerc)
		{
			print FH "$nkey\n";
		}
		else
		{
			last;
		}
	}
close FH;
#system("cp -rf $stTSV.FilteredBy$stDomain[$i]domain.$cut $RUNNING_PATH/");

#	$stData = undef;
	my $nPerc ="100";	###	percent
	my $stData = $stTSV;
	my $stPfam = $stPfamID[$i];
	my @NB=();
	my $num=0;

open FH1, ">", "$RUNNING_PATH/$stTSVOutput.FilteredBy$stDomain[$i]domain.all" || die ("Can't open myfile");
	open(DATA1, "$stData");
	while(my $stLine = <DATA1>)
	{
		chomp($stLine);
		my @stList = split /[\t]/, $stLine;
		if(($stList[0] ne "")&&($stList[4] eq "$stPfam"))
		{
			push(@NB, "$stLine");	
		}
	}
	close DATA1;

	my @sort= sort{(split /[\t]/, $a)[8] <=> (split /[\t]/, $b)[8]}@NB;
	for(my $i=0; $i<@sort; $i++) ### 30 or @sort
	{
		print FH1 "$sort[$i]\n";
	}
close FH1;

#system("cp -rf $stTSV.FilteredBy$stDomain[$i]domain.all $RUNNING_PATH/");

my $stData = "$RUNNING_PATH/$stTSVOutput.FilteredBy$stDomain[$i]domain.$nPerc"; # NB-ARC.tsv file such as 1_0.at_10.fasta.onlyNB-ARC.sorted.tsv
my $stPEP = $PROTEINS_FOR_DOMAIN_IDENTIFICATION;
my %ID;
my $nbID="";
my $count=0;
my $str="";
my $end="";
my $seq="";
my @str=();
my @end=();
my $stSeq="";


open FH2, ">", "$stPrefix.$stDomain[$i]domain.$nPerc.fasta.temp" || die ("Can't open myfile");
	open(DATA, $stData);
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		my @stList = split /[\t]/, $stLine;
		$ID{$stList[0]}=$ID{$stList[0]}."$stList[6]	$stList[7],";	
	}
	close (DATA);

	open(DATA, $stPEP);
	while(my $stLine = <DATA>)
	{
			chomp($stLine);
			if($stLine =~ />([^\s]+)/)
			{
			if($stSeq ne "")
			{
					if($count == 1)
					{
					if($str eq "")
					{
						for(my $j=0; $j<@str; $j++)
						{
							$seq= substr($stSeq, $str[$j]-1, $end[$j]-$str[$j]+1);
							print FH2 ">$nbID\.$j\n$seq\n";
						}
						@str=();
						@end=();
					}
					else
					{
						$seq= substr($stSeq, $str-1, $end-$str+1);                
						print FH2 ">$nbID\.0\n$seq\n";
						$str="";
						$end="";
					}
				}
				$stSeq="";
				$count=0;
			}
					$nbID = $1;
					if($ID{$nbID} ne "")
					{
							$ID{$nbID} =~ s/,$//g;

				if($ID{$nbID} =~ /,/)
				{
					my @devide = split /,/, $ID{$nbID};
					for(my $i=0; $i<@devide; $i++)
					{
						my @StrEnd = split /[\t]/, $devide[$i];
						push(@str, "$StrEnd[0]");
						push(@end, "$StrEnd[1]");	
					}	
					$count++;
				}
				else
				{
					my @StrEnd = split /[\t]/, $ID{$nbID};
					$str= $StrEnd[0];
					$end= $StrEnd[1];
								$count++;
				}
					}
			}       
		else
		{
			$stSeq=$stSeq."$stLine";
		}
	}
	if($stSeq ne "")
	{
			if($count == 1)
			{
			if($str eq "")
			{
				for(my $j=0; $j<@str; $j++)
				{
					$seq= substr($stSeq, $str[$j]-1, $end[$j]-$str[$j]+1);
					print FH2 ">$nbID\.$j\n$seq\n";
				}
				@str=();
				@end=();
			}
			else
			{
				$seq= substr($stSeq, $str-1, $end-$str+1);                
				print FH2 ">$nbID\.0\n$seq\n";
				$str="";
				$end="";
			}
		}
		$stSeq="";
		$count=0;
	}
	close DATA;
close FH2;

open FH3, ">$stPrefix."."$stDomain[$i]"."domain.$nPerc.fasta" || die ("Can't open myfile");
	my $stData="$stPrefix.$stDomain[$i]domain.$nPerc.fasta.temp";	###	Chi.Final.GeneModel.Whole.PEP.fa

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
		print FH3 "$stList{$stKey}\n$stKey\n";
	}
close FH3;


	my $stData = "$RUNNING_PATH/$stTSVOutput.FilteredBy$stDomain[$i]domain.all"; # NB-ARC.tsv file such as 1_0.at_10.fasta.onlyNB-ARC.sorted.tsv
	my $stPEP = $PROTEINS_FOR_DOMAIN_IDENTIFICATION;
	my %ID;
	my $nbID="";
	my $count=0;
	my $stSeq="";

	open (FH4, ">$stPrefix.geneIncluding$stDomain[$i].fasta.temp") || die ("Can't open myfile");

	open(DATA, "$stData");
	while(my $stLine = <DATA>)
	{
			chomp($stLine);
			my @stList = split /[\t]/, $stLine;
			$ID{$stList[0]}=$ID{$stList[0]}."$stList[6]	$stList[7],";
	}
	close (DATA);

	open(DATA, "$stPEP");
	while(my $stLine2 = <DATA>)
	{
			chomp($stLine2);
			if($stLine2 =~ />([^\s]+)/)
			{
			if($stSeq ne "")
			{
					if($count == 1)
					{
							print FH4 ">$nbID\n$stSeq\n";
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
					print FH4 ">$nbID\n$stSeq\n";
		}
		$count=0;
		$stSeq="";
	}
	close (DATA);
close(FH4);

	my $stData="$stPrefix.geneIncluding$stDomain[$i].fasta.temp";

open (FH5, ">$stPrefix.geneIncluding$stDomain[$i].fasta") || die ("Can't open myfile");

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
		print FH5 "$stList{$stKey}\n$stKey\n";
	}
close(FH5);

	system("rm -rf $stPrefix.geneIncluding"."$stDomain[$i]".".fasta.temp $stPrefix."."$stDomain[$i]"."domain.$nPerc.fasta.temp");

###	ClustalW
	system("$CLUSTALW_PATH/clustalw2 -ITERATION=ALIGNMENT -ALIGN -type=protein -NUMITER=10000 -INFILE=$stPrefix."."$stDomain[$i]"."domain.$nPerc.fasta -OUTFILE=$stPrefix."."$stDomain[$i]"."domain.$nPerc.result");

	my $stData="$stPrefix.$stDomain[$i]domain.$nPerc.result";
	my $stOut=$stData;
	my $nIdx1=0;
	my %stList;
	my $IDlen=0;
	$stOut =~ s/\.result/\.msa/;

	open(DATA, "$stData");

	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		$nIdx1++;
		next if($nIdx1 == 1);
		if($stLine =~ /^[^\s]+/)
		{
			my @stInfo=split(/[\s]+/,$stLine);
			$stList{$stInfo[0]}=$stList{$stInfo[0]}.$stInfo[1];
			if($IDlen < length($stInfo[0]))
			{
				$IDlen = length($stInfo[0]);
			}
		}
	}
	close DATA;

open (FH6, ">$stOut") || die ("Can't open myfile");;
	print FH6 "# STOCKHOLM 1.0\n\n#=GF ID   NoName\n\n";
	foreach my $nkey (keys %stList)
	{
		my $add = $IDlen - length($nkey);
		my $empty = " " x ($add + 6);
		print FH6 "$nkey$empty$stList{$nkey}\n";
	}
	print FH6 "//\n";
close (FH6);

###	HMM build	###
	system("$HMMER_BIN_PATH/hmmbuild $stPrefix."."$stDomain[$i]"."domain.hmm $stPrefix."."$stDomain[$i]"."domain.$nPerc.msa");

###	HMM search	###
	system("$HMMER_BIN_PATH/hmmsearch $stPrefix."."$stDomain[$i]"."domain.hmm $stPrefix.SixFrameTranslation.PEP.fasta > $stPrefix."."$stDomain[$i]"."domain.search.out");
	$stMerge=$stMerge."$stPrefix."."$stDomain[$i]"."domain.search.out ";
}

###	Cat multiple hmm search result	###
sleep 3;
system("cat $stMerge> $stPrefix.Final."."$stRepreDom"."domain.search.out");

	my $stData = "$stPrefix.Final.$stRepreDom\domain.search.out";
	my $ID = "";
	my $num=0;

open FH7, ">", "$stPrefix.Final.$stRepreDom\domain.search.out.filter" || die ("Can't open myfile");
	open(DATA, "$stData");
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		if($stLine =~ /^>>/)
		{
			$ID = $stLine;
			$ID =~ s/>> //g;
			$ID =~ s/[\s]+$//g;	
		}
		elsif($stLine =~ /---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----/)
		{
			$num++;	
		}
		elsif(($num == 1)&&($stLine ne ""))
		{
			my @stList = split /[\s]+/, $stLine;
			my @stInfo = split /_/, $ID;
			
			if($stInfo[4] eq "+")
			{
				my $aminoStr = ($stList[10] - 1) * 3 + 1;
				my $aminoEnd = $stList[11] * 3;
				my $str = $stInfo[2] + $aminoStr - 1;
				my $end = $stInfo[2] + $aminoEnd - 1; 
				print FH7 "$stInfo[0]_$stInfo[1]_$stInfo[2]_$stInfo[3]_$stInfo[4]	Domain_position	$stInfo[0]	$str	$end\n";
			}
			elsif($stInfo[4] eq "-")
			{
				my $aminoStr = ($stList[10] - 1) * 3 + 1;
							my $aminoEnd = $stList[11] * 3;
							my $str = $stInfo[3] - $aminoEnd + 1;
							my $end = $stInfo[3] - $aminoStr + 1;          
							print FH7 "$stInfo[0]_$stInfo[1]_$stInfo[2]_$stInfo[3]_$stInfo[4]	Domain_position	$stInfo[0]	$str	$end\n";
			}
			else
			{
				print FH7 "error	$stLine\n"
			}
		}
		elsif($stLine eq "")
		{
			$num=0;
		}
	}
	close DATA;
close FH7;

	my $stData = "$stPrefix.Final.$stRepreDom\domain.search.out";
	my $num=0;
	my $sum=0;

	open(DATA, "$stData");
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		if($stLine =~ /------- ------ -----/)
		{
			$num++;
		}	
		elsif(($num == 1)&&($stLine !~ /------ inclusion threshold ------/))
		{
			if($stLine eq "")
			{
				last;
			}
			my @stList = split /[\s]+/, $stLine;
			$sum = $sum + $stList[8];
		}
		
	}
	close DATA;
	#	print "Validation_HMMsearch_out	$sum\n";

	my $stData = "$stPrefix.Final.$stRepreDom\domain.search.out.filter";
	my @NB=();
	my $num=0;

open FH8, ">", "$stPrefix.Final.$stRepreDom\domain.search.out.filter.sorted" || die ("Can't open myfile");
	open(DATA, "$stData");
	while(my $stLine = <DATA>)
	{
			chomp($stLine);
		push(@NB, "$stLine");   
	}
	close DATA;

	my @sort= sort { ($a =~ /	tempCh([0-9]+)/)[0] <=> ($b =~ /	tempCh([0-9]+)/)[0] || (split /[\t]+/, $a)[3] <=> (split /[\t]+/, $b)[3]}@NB;
	for(my $i=0; $i<@sort; $i++)
	{
			print FH8 "$sort[$i]\n";
	}
close FH8;

	my $stData = "$stPrefix.Final.$stRepreDom\domain.search.out.filter.sorted";
	my $str="";
	my $end="";

open FH9, ">", "$stPrefix.Final.$stRepreDom\domain.search.out.Extended" || die ("Can't open myfile");
	open(DATA, "$stData");
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		my @stList = split /[\t]+/, $stLine;
		$str = $stList[3] - $EXTENSION_LENGTH;
		$end = $stList[4] + $EXTENSION_LENGTH;
		if($str < 0)
		{
			$str = 1;
		}

		print FH9 "Extended	$stList[0]	$stList[1]	$stList[2]	$str	$end\n";
		
		if(($str < 0)||($end < 0))
		{
			print FH9 "error	$stLine\n";
			last;
		}
		
	}
	close DATA;
close FH9;

	my $stData = "$stPrefix.Final.$stRepreDom\domain.search.out.Extended";
	my $str="";
	my $end="";
	my $num=0;
	my $ID="";

open FH10, ">", "$stPrefix.Final.$stRepreDom\domain.search.out.Extended.RemoveRedundant" || die ("Can't open myfile");

	open(DATA, "$stData");
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		$num++;
		my @stList = split /[\t]/, $stLine;

		if($num != 1)
		{
			if($stList[3] eq "$ID")
			{
				if($stList[4] <= $end)
				{
					if($stList[5] > $end)
					{
						$end= $stList[5];
					}
				}
				else
				{	
					print FH10 "$ID	$str	$end\n";
					$str= $stList[4];
					$end= $stList[5];		
				}
			}
			else
			{
				print FH10 "$ID	$str	$end\n";
				$str= $stList[4];
				$end= $stList[5];
			}
			$ID = $stList[3];
		}
		else
		{
			$str= $stList[4];
			$end= $stList[5];
			$ID= $stList[3]; 
		}
	}
	close DATA;

	print FH10 "$ID	$str	$end\n";
close FH10;

###	Genome masking	###
	my $stData = "$stPrefix.Final.$stRepreDom\domain.search.out.Extended.RemoveRedundant";
	my $stGenome = "$stPrefix.tempID.fasta";
	my %Info;
	my $stID="";
	my $num=0;
	my $exEnd="";
	my $seq="";
	my $length="";


	open(DATA, "$stData");
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		my @stList = split /[\t]/, $stLine;
		$Info{$stList[0]} = $Info{$stList[0]}."$stList[1]	$stList[2],";	
	}
	close DATA;


open FH11, ">", "$stPrefix.tempID.Masked.Except"."$stRepreDom".".fasta" || die ("Can't open myfile");
	open(DATA, "$stGenome");
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		if($stLine =~ />/)
		{
			print FH11 "$stLine\n";
			$stLine =~ />([^\s]+)/; 
			$stID = $1;
		}
		else
		{
			$seq = $stLine;
			if($Info{$stID} ne "")
			{
				$Info{$stID} =~ s/,$//g;
				my @StrEnd = split /,/, $Info{$stID};
		
				for(my $i=0; $i<@StrEnd; $i++)
				{
					$num++;
					my @stList = split /[\t]/, $StrEnd[$i];
					if($num == 1)
					{
						$length = $stList[0]-1;
						substr ($seq, 0, $length) = "X" x ($length);
						$exEnd = $stList[1];
					}
					else
					{
						$length = $stList[0] - $exEnd +1 -2; #except both side nt
						substr ($seq, $exEnd, $length) = "X" x ($length); #$exEnd-1+1 (substr starts 0, repeat masking except NB-ARC region)
						$exEnd = $stList[1];
					}	
				}
				$num=0;

				if(length($seq) > $exEnd)
				{
					$length = length($seq) - $exEnd;
					substr ($seq, $exEnd, $length) = "X" x ($length);
					print FH11 "$seq\n";
				}
				else
				{
					print FH11 "$seq\n";
				}
			}
			else
			{
				$length = length($seq);
				substr ($seq, 0, $length) = "X" x ($length);
				print FH11 "$seq\n";	
			}
		}
	}
	close DATA;
close FH11;

open FH12, ">", "$stPrefix.tempID.Masked.Except$stRepreDom.fasta.validation" || die ("Can't open myfile");
	my $stData = "$stPrefix.tempID.Masked.Except$stRepreDom.fasta";
	my $stID="";
	my $num=0;
	my $num2=0;
	my $num3=0;

	open(DATA, "$stData");
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		if($stLine =~ />([^\s]+)/)
		{
			$stID = $1;
		}
		else
		{
			my @stList = split //, $stLine;
			for(my $i=0; $i<@stList; $i++)
			{
				$num++;
				if($stList[$i] eq "X")
				{
					if($num3 > 1)
					{
						my $end= $num-1;
						print FH12 "$end\n";
						$num3=0;
						$num2=0;
					}
					next;
				}
				else
				{
					$num2++;
					if($num2 == 1)
					{
						print FH12 "$stID	$num	";
					}
					else
					{
						$num3++;
						if($i eq $#stList)
						{
							print FH12 "$num\n"
						}
					}
					
				}
			}
		}
		$num=0;
		$num2=0;
		$num3=0;	

	}
	close DATA;
close FH12;

open (FH13, ">$stPrefix.tempID.Masked.Except"."$stRepreDom".".fasta.Main.fa") || die ("Can't open myfile");

	my $stData = "$stPrefix.tempID.Masked.Except"."$stRepreDom".".fasta";
	my $stID="";

	open(DATA, "$stData");
	while(my $stLine = <DATA>)
	{
		chomp($stLine);
		if($stLine =~ /^>/)
		{
			$stID = $stLine;		
		}
		else
		{
			if($stLine =~ /[^X]/)
			{
				print FH13 "$stID\n$stLine\n";
			}
		}
	}
	close DATA;
close(FH13);
system("$BLAST_BIN_PATH/makeblastdb -in $RESOURCE_PROTEIN -dbtype prot -out $BLAST_DB_NAME");
print "\n############ 1.Domain_Identification.pl is finished ############\n\n";
