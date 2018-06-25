use lib $ARGV[0];

BEGIN {
	require $ARGV[1];
		$ARGV[1]->import( qw(TGFAM_CONFIG) );
	require $ARGV[2];
		$ARGV[2]->import( qw(TGFAM_PATH) );
	}

use strict;

system("cp -rf $RUNNING_PATH/$PROTEIN_MAPPING_ANALYSIS_PATH/temp*/*RefPEP.Protein.sorted.PosiModi.gff3 $RUNNING_PATH/$PROTEIN_MAPPING_ANALYSIS_PATH/temp*/*RefPEP.Mapped.CDS.fa $RUNNING_PATH/$AUGUSTUS_ANALYSIS_PATH/*.filter $RUNNING_PATH/$PM_AUGUSTUS_ANALYSIS_PATH");

my $stExonerateGff = "$OUTPUT_PREFIX.RefPEP.Protein.sorted.PosiModi.gff3"; #RefPEP.Protein.sorted.PosiModi.gff3
my $stExonerateCDS = "$OUTPUT_PREFIX.RefPEP.Mapped.CDS.fa"; 
my $stAugustusGff = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.augustus.gff3.filter";
my $stAugustusCDS = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.augustus.CDS.fa.filter";

my $stTemp = "gene";
my $stID;
my %stSeq;
my $nCnt = 0;
my %stAugPM;
system("grep $stTemp $stAugustusGff > Temp.Info");
system("cat $stExonerateCDS $stAugustusCDS > Temp.CDS");

open (FH, ">$OUTPUT_PREFIX.PM_Augustus.list") || die ("Can't open myfile");
open(DATA, "Temp.CDS");
while(my $stLine = <DATA>)
{
	chomp($stLine);
	if($stLine =~ /^>([^\s]+)/) 
	{
		$stID = $1;
	}
	else
	{
		$stSeq{$stID} = $stLine;
	}
}
close(DATA);

open(DATA, "$stExonerateGff");
open(OUT, ">>Temp.Info");
while(my $stLine = <DATA>)
{
	chomp($stLine);
	my @stInfo = split /[\t]+/, $stLine;
	next if ($stInfo[2] ne "gene");
	$stInfo[$#stInfo] =~ s/-Exonerate-prediction-.+/_$nCnt/g;
	$nCnt++;
	$stLine = join("	",@stInfo);
	print OUT "$stLine\n";
}
close(DATA);
close(OUT);


open(DATA, "Temp.Info");
my @stTotal = <DATA>;
chomp(@stTotal);
close(DATA);

@stTotal = sort
{
	($a =~ /^([^\s]+)/)[0] cmp ($b =~ /^([^\s]+)/)[0]||
	($a =~ /^[^\s]+	[^\s]+	[^\s]+	([0-9]+)/)[0] <=> ($b =~ /^[^\s]+	[^\s]+	[^\s]+	([0-9]+)/)[0]||
	($b =~ /^[^\s]+	[^\s]+	[^\s]+	[0-9]+	([0-9]+)/)[0] <=> ($a =~ /^[^\s]+	[^\s]+	[^\s]+	[0-9]+	([0-9]+)/)[0]||
	($a =~ /^[^\s]+	([^\s]+)/)[0] cmp ($b =~ /^[^\s]+	([^\s]+)/)[0]
}@stTotal;

foreach(@stTotal)
{
#	print "$_\n";
}

for(my $i=0; $i<@stTotal; $i++)
{
	my @stInfo1 = split /[\t]+/, $stTotal[$i];
	next if ($stInfo1[1] ne "AUGUSTUS");

	for(my $j=$i-1; $j<@stTotal; $j++)
	{
		my @stInfo2 = split /[\t]+/, $stTotal[$j];
		next if($stInfo2[1] ne "REF");
		next if ($stInfo2[4]<$stInfo1[3]);
		if (($stInfo1[3]<=$stInfo2[3])&&($stInfo1[4]>=$stInfo2[4]))
		{
			my $stID1= $stInfo1[$#stInfo1];
			my $stID2= $stInfo2[$#stInfo2];
			$stID1 =~ s/ID=//g;
			$stID2 =~ s/ID=//g;
			if($stSeq{$stID1} =~ /$stSeq{$stID2}/)
			{
				$stAugPM{$stID1} = $stAugPM{$stID1}."$stID2,";
			}
			else
			{
	#			print "$stID1	$stID2\n";
			}
		}
		if($stInfo1[0] ne $stInfo2[0])
		{
			last;
		}
		if($stInfo1[4]<$stInfo2[3])
		{
			last;
		}
	}
}

foreach my $stKey (keys %stAugPM)
{
	next if ($stAugPM{$stKey} eq "");
	print FH "$stKey	$stAugPM{$stKey}\n";
}

close(FH);

my $stList = "$OUTPUT_PREFIX.PM_Augustus.list";
my $stPEP = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.augustus.PEP.fa.filter";
my $stCDS = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.augustus.CDS.fa.filter";
my $stGff = "$OUTPUT_PREFIX\_$REPRESENTATIVE_DOMAIN_NAME.augustus.gff3.filter";
my $stPrefix = "PM.Augustus";

my $stOut = $stPEP;
$stOut =~ s/.PEP.fa.filter//g;
my $stName;
my %stList = {};
my $stID = "";
my $nCnt1=0;
my $nCnt2=0;
my $nIdx =0;
my $nCnt = -1;

open(DATA, "$stList");
while(my $stLine = <DATA>)
{
	chomp($stLine);
	my @stInfo = split /[\t]+/, $stLine;
	$stList{$stInfo[0]}++;
}
close(DATA);

open(DATA, "$stPEP");
open(OUT, ">$stOut.PEP.fa.$stPrefix");
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
		if ($stList{$stName} ne "")
		{
                        print OUT ">$stName\n$stSeq\n";
					
		}
	}
}
close(DATA);


open(DATA, "$stCDS");
open(OUT, ">$stOut.CDS.fa.$stPrefix");
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
open(OUT, ">$stOut.gff3.$stPrefix");
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

system("cp *.PM.Augustus $RUNNING_PATH/$MERGING_ANALYSIS_PATH");
