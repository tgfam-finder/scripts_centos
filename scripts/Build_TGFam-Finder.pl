#!/usr/bin/perl
use lib $ARGV[0];

BEGIN {
        require $ARGV[1];
                $ARGV[1]->import( qw(add) );
        require $ARGV[2];
                $ARGV[2]->import( qw(add1) );
        }
use strict;
use warnings;

## CHECK FOR PROGRAM PATH ##

print "\n*********************************** PROGRAM PATH *********************************** \n";
print "* \n";

  if (-e "$TGFAM_SCRIPTS_PATH/0.SixFrameTranslation.pl") {
  } else {
    print "* Error: We could't find scripts file in TGFAM_SCRIPTS_PATH of PROGRAM_PATH.config \n";
	print "* \n";
	print "************************************************************************************ \n\n";
    exit 1;
  }

  if ($RUNNING_PATH ne '') {
  } else {
    print "* Error: Empty information of RUNNING_PATH in PROGRAM_PATH.config \n";
	print "* \n";
	print "************************************************************************************ \n\n";
    exit 1;
  }

  if (-e "$CLUSTALW_PATH/clustalw2") {
  } else {
    print "* Error: We could't find execution file in CLUSTALW_PATH of PROGRAM_PATH.config \n";
  }
  if (-e "$HMMER_BIN_PATH/hmmsearch") {
  } else {
    print "* Error: We could't find execution file in HMMER_BIN_PATH of PROGRAM_PATH.config \n";
  }
  if (-e "$BLAST_BIN_PATH/blastx") {
  } else {
    print "* Error: We could't find execution file in BLAST_BIN_PATH of PROGRAM_PATH.config \n";
  }
  if (-e "$IPRSCAN_PATH/interproscan.sh") {
  } else {
    print "* Error: We could't find execution file in IPRSCAN_PATH of PROGRAM_PATH.config \n";
  }
  if (-e "$BOWTIE_PATH/bowtie2-build") {
  } else {
    print "* Error: We could't find execution file in BOWTIE_PATH of PROGRAM_PATH.config \n";
  }
  if (-e "$TOPHAT_PATH/tophat") {
  } else {
    print "* Error: We could't find execution file in TOPHAT_PATH of PROGRAM_PATH.config \n";
  }
  if (-e "$EXONERATE_BIN_PATH/exonerate") {
  } else {
    print "* Error: We could't find execution file in EXONERATE_BIN_PATH of PROGRAM_PATH.config \n";
  }
  if (-e "$AUGUSTUS_PATH/scripts/new_species.pl") {
  } else {
    print "* Error: We could't find execution file for training in AUGUSTUS_PATH of PROGRAM_PATH.config \n";
  }
  if (-e "$AUGUSTUS_BIN_PATH/augustus") {
  } else {
    print "* Error: We could't find augustus execution file in AUGUSTUS_BIN_PATH of PROGRAM_PATH.config \n";
  }
  if (-d "$AUGUSTUS_CONFIG_PATH") {
  } else {
    print "* Error: We could't find config directory in AUGUSTUS_CONFIG_PATH of PROGRAM_PATH.config \n";
  }
  if (-e "$SCIPIO_PATH/yaml2gff.1.4.pl") {
  } else {
    print "* Error: We could't find execution file in SCIPIO_PATH of PROGRAM_PATH.config \n";
  }
  if (-e "$BLAT_PATH/blat") {
  } else {
    print "* Error: We could't find execution file in BLAT_PATH of PROGRAM_PATH.config \n";
  }
  if (-e "$CUFFLINKS_PATH/cufflinks") {
  } else {
    print "* Error: We could't find execution file in CUFFLINKS_PATH of PROGRAM_PATH.config \n";
  }

  if ((-e ("$TGFAM_SCRIPTS_PATH/0.SixFrameTranslation.pl") && -e ("$CLUSTALW_PATH/clustalw2") && -e ("$HMMER_BIN_PATH/hmmsearch") && -e ("$BLAST_BIN_PATH/blastx") && -e ("$IPRSCAN_PATH/interproscan.sh") && -e ("$BOWTIE_PATH/bowtie2-build") && -e ("$TOPHAT_PATH/tophat") && -e ("$EXONERATE_BIN_PATH/exonerate") && -e ("$AUGUSTUS_PATH/scripts/new_species.pl") && -e ("$AUGUSTUS_BIN_PATH/augustus") && -e ("$SCIPIO_PATH/yaml2gff.1.4.pl") && -e ("$BLAT_PATH/blat") && -e ("$CUFFLINKS_PATH/cufflinks") && ($RUNNING_PATH ne '') && -d ("$AUGUSTUS_CONFIG_PATH"))) {
	print "* The required information in 'PROGRAM_PATH' has been verified. \n";
  } else {
  }

print "* \n";
print "************************************************************************************ \n\n";


## CHECK FOR RESOURCE FILE ##

print "******************************** RESOURCE (Required) ******************************* \n";
print "* \n";

  if (-e "$TARGET_GENOME") {
  } else {
    print "* Error: We could't find TARGET_GENOME in RESOURCE.config \n";
  }
  if (-e "$PROTEINS_FOR_DOMAIN_IDENTIFICATION") {
  } else {
    print "* Error: We could't find PROTEINS_FOR_DOMAIN_IDENTIFICATION in RESOURCE.config \n";
  }
  if (-e "$TSV_FOR_DOMAIN_IDENTIFICATION") {
  } else {
    print "* Error: We could't find TSV_FOR_DOMAIN_IDENTIFICATION in RESOURCE.config \n";
  }
  if (-e "$RESOURCE_PROTEIN") {
  } else {
    print "* Error: We could't find RESOURCE_PROTEIN in RESOURCE.config \n";
  }
  if ($BLAST_DB_NAME ne '') {
  } else {
    print "* Error: You should input BLAST_DB_NAME in RESOURCE.config \n";
  }
  if ($OUTPUT_PREFIX ne '') {
  } else {
    print "* Error: You should input OUTPUT_PREFIX in RESOURCE.config \n";
  }
  if ($TARGET_DOMAIN_ID ne '') {
  } else {
    print "* Error: You should input TARGET_DOMAIN_ID in RESOURCE.config \n";
  }
  if ($TARGET_DOMAIN_NAME ne '') {
  } else {
    print "* Error: You should input TARGET_DOMAIN_NAME in RESOURCE.config \n";
  }
  if ($REPRESENTATIVE_DOMAIN_NAME ne '') {
  } else {
    print "* Error: You should input REPRESENTATIVE_DOMAIN_NAME in RESOURCE.config \n";
  }
  if ($EXTENSION_LENGTH ne '' && ($EXTENSION_LENGTH !~ /\D/)) {
  } else {
    print "* Error: You should input EXTENSION_LENGTH in RESOURCE.config \n";
  }
  if ($MAX_INTRON_LENGTH ne '' && ($MAX_INTRON_LENGTH !~ /\D/)) {
  } else {
    print "* Error: You should input MAX_INTRON_LENGTH in RESOURCE.config \n";
  }
  if ($THREADS ne '' && ($THREADS !~ /\D/)) {
  } else {
    print "* Error: You should input THREADS in RESOURCE.config \n";
  }

  if (-e ("$TARGET_GENOME" && "$PROTEINS_FOR_DOMAIN_IDENTIFICATION" && "$TSV_FOR_DOMAIN_IDENTIFICATION" && "$RESOURCE_PROTEIN") && ($BLAST_DB_NAME ne '') && ($OUTPUT_PREFIX ne '') && ($TARGET_DOMAIN_ID ne '') && ($TARGET_DOMAIN_NAME ne '') && ($REPRESENTATIVE_DOMAIN_NAME ne '') && ($EXTENSION_LENGTH ne '' && ($EXTENSION_LENGTH !~ /\D/)) && ($MAX_INTRON_LENGTH ne '' && ($MAX_INTRON_LENGTH !~ /\D/)) && ($THREADS ne '' && ($THREADS !~ /\D/)) ) {
	print "* The required values in 'RESOURCE' has been verified. \n";
  } else {
  }

print "* \n";
print "************************************************************************************ \n\n";

## CHECK FOR OPTIONAL INFORMATION ##

print "******************************* RESOURCE (Optional) ******************************** \n";
print "* \n";

  if (-e "$CDS_OF_TARGET_GENOME") {
  } else {
    print "* Warning: Incorrect or empty information of CDS_OF_TARGET_GENOME in RESOURCE.config\n";
  }
  if (-e "$GFF3_OF_TARGET_GENOME") {
  } else {
    print "* Warning: Incorrect or empty information of GFF3_OF_TARGET_GENOME in RESOURCE.config \n";
  }
  if (-e "$RNASEQ_FORWARD_PATH") {
  } else {
    print "* Warning: Incorrect or empty information of RNASEQ_FORWARD_PATH in RESOURCE.config \n";
  }
  if (-e "$RNASEQ_REVERSE_PATH") {
  } else {
    print "* Warning: Incorrect or empty information of RNASEQ_REVERSE_PATH in RESOURCE.config \n";
  }

  if ($EXCLUDED_DOMAIN_ID ne '') {
  } else {
    print "* Warning: Incorrect or empty information of EXCLUDED_DOMAIN_ID in RESOURCE.config \n";
  }

if (-e ("$CDS_OF_TARGET_GENOME" && "$GFF3_OF_TARGET_GENOME" && "$RNASEQ_FORWARD_PATH" && "$RNASEQ_REVERSE_PATH") && ($EXCLUDED_DOMAIN_ID ne '')) {
	print "* The optional values in 'RESOURCE' has been verified. \n";
} else {
  }

print "* \n";
print "************************************************************************************ \n\n";

## CHECK FOR PATH OF THE ANALYSIS RESULT FOLDER ##

print "*********************** PROGRAM PATH (Analysis result folder) ********************** \n";
print "* \n";

  if (-d "$RUNNING_PATH/$AUGUSTUS_ANALYSIS_PATH") {
    print "* Error: A folder having same name is located in AUGUSTUS_ANALYSIS_PATH \n";
  }
  if (-d "$RUNNING_PATH/$ISGAP_ANALYSIS_PATH") {
    print "* Error: A folder having same name is located in ISGAP_ANALYSIS_PATH \n";
  }
  if (-d "$RUNNING_PATH/$PROTEIN_MAPPING_ANALYSIS_PATH") {
    print "* Error: A folder having same name is located in PROTEIN_MAPPING_ANALYSIS_PATH \n";
  }
  if (-d "$RUNNING_PATH/$RNASEQ_AUGUSTUS_ANALYSIS_PATH") {
    print "* Error: A folder having same name is located in RNASEQ_AUGUSTUS_ANALYSIS_PATH \n";
  }
  if (-d "$RUNNING_PATH/$PM_AUGUSTUS_ANALYSIS_PATH") {
    print "* Error: A folder having same name is located in PM_AUGUSTUS_ANALYSIS_PATH \n";
  }
  if (-d "$RUNNING_PATH/$ISGAP_AUGUSTUS_ANALYSIS_PATH") {
    print "* Error: A folder having same name is located in ISGAP_AUGUSTUS_ANALYSIS_PATH \n";
  }
  if (-d "$RUNNING_PATH/$MERGING_ANALYSIS_PATH") {
    print "* Error: A folder having same name is located in MERGING_ANALYSIS_PATH \n";
  }
  if (-d "$RUNNING_PATH/$FINAL_ANALYSIS_PATH") {
    print "* Error: A folder having same name is located in FINAL_ANALYSIS_PATH \n";
  }

  if (-d ("$RUNNING_PATH/$AUGUSTUS_ANALYSIS_PATH" || "$RUNNING_PATH/$ISGAP_ANALYSIS_PATH" || "$RUNNING_PATH/$PROTEIN_MAPPING_ANALYSIS_PATH" || "$RUNNING_PATH/$RNASEQ_AUGUSTUS_ANALYSIS_PATH" || "$RUNNING_PATH/$PM_AUGUSTUS_ANALYSIS_PATH" || "$RUNNING_PATH/$ISGAP_AUGUSTUS_ANALYSIS_PATH" || "$RUNNING_PATH/$MERGING_ANALYSIS_PATH" || "$RUNNING_PATH/$FINAL_ANALYSIS_PATH")) {
	  
	print "* \n";
	print "************************************************************************************ \n\n";
    exit 1;

  }

  else {

	system("cp $TGFAM_SCRIPTS_PATH/setup_shell Run_TGFam-Finder.sh");
	system("perl -pi -e 's|CLUSTALW_PATH|$CLUSTALW_PATH|g' Run_TGFam-Finder.sh");
	system("perl -pi -e 's|HMMER_BIN_PATH|$HMMER_BIN_PATH|g' Run_TGFam-Finder.sh");
	system("perl -pi -e 's|BLAST_BIN_PATH|$BLAST_BIN_PATH|g' Run_TGFam-Finder.sh");
	system("perl -pi -e 's|IPRSCAN_PATH|$IPRSCAN_PATH|g' Run_TGFam-Finder.sh");
	system("perl -pi -e 's|BOWTIE_PATH|$BOWTIE_PATH|g' Run_TGFam-Finder.sh");
	system("perl -pi -e 's|TOPHAT_PATH|$TOPHAT_PATH|g' Run_TGFam-Finder.sh");
	system("perl -pi -e 's|EXONERATE_BIN_PATH|$EXONERATE_BIN_PATH|g' Run_TGFam-Finder.sh");
	system("perl -pi -e 's|AUGUSTUS_PATH|$AUGUSTUS_PATH|g' Run_TGFam-Finder.sh");
	system("perl -pi -e 's|AUGUSTUS_BIN_PATH|$AUGUSTUS_BIN_PATH|g' Run_TGFam-Finder.sh");
	system("perl -pi -e 's|AUGUSTUS_CONFIG_PATH|AUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH|g' Run_TGFam-Finder.sh");
	system("perl -pi -e 's|SCIPIO_PATH|$SCIPIO_PATH|g' Run_TGFam-Finder.sh");
	system("perl -pi -e 's|BLAT_PATH|$BLAT_PATH|g' Run_TGFam-Finder.sh");
	system("perl -pi -e 's|CUFFLINKS_PATH|$CUFFLINKS_PATH|g' Run_TGFam-Finder.sh");
	system("perl -pi -e 's|ISGAP_ANALYSIS_PATH|$ISGAP_ANALYSIS_PATH|g' Run_TGFam-Finder.sh");
	system("perl -pi -e 's|MERGING_ANALYSIS_PATH|$MERGING_ANALYSIS_PATH|g' Run_TGFam-Finder.sh");
	system("perl -pi -e 's|PROTEIN_MAPPING_ANALYSIS_PATH|$PROTEIN_MAPPING_ANALYSIS_PATH|g' Run_TGFam-Finder.sh");
	system("perl -pi -e 's|RNASEQ_AUGUSTUS_ANALYSIS_PATH|$RNASEQ_AUGUSTUS_ANALYSIS_PATH|g' Run_TGFam-Finder.sh");
	system("perl -pi -e 's|PM_AUGUSTUS_ANALYSIS_PATH|$PM_AUGUSTUS_ANALYSIS_PATH|g' Run_TGFam-Finder.sh");
	system("perl -pi -e 's|ISGAP_AUGUSTUS_ANALYSIS_PATH|$ISGAP_AUGUSTUS_ANALYSIS_PATH|g' Run_TGFam-Finder.sh");
	system("perl -pi -e 's|TGFAM_SCRIPTS_PATH|$TGFAM_SCRIPTS_PATH|g' Run_TGFam-Finder.sh");
	system("perl -pi -e 's|RUNNING_PATH|$RUNNING_PATH|g' Run_TGFam-Finder.sh");
	system("perl -pi -e 's|BLAST_DB_NAME|$BLAST_DB_NAME|g' Run_TGFam-Finder.sh");
	system("perl -pi -e 's|RESOURCE_PROTEIN|$RESOURCE_PROTEIN|g' Run_TGFam-Finder.sh");
	system("perl -pi -e 's|DOMAIN_NAME|$REPRESENTATIVE_DOMAIN_NAME|g' Run_TGFam-Finder.sh");
	system("perl -pi -e 's|OUTPUT_PREFIX|$OUTPUT_PREFIX|g' Run_TGFam-Finder.sh");
	system("perl -pi -e 's|AUGUSTUS_ANALYSIS_PATH|$AUGUSTUS_ANALYSIS_PATH|g' Run_TGFam-Finder.sh");
	system("perl -pi -e 's|CONF1|$ARGV[0]|g' Run_TGFam-Finder.sh");
	system("perl -pi -e 's|CONF2|$ARGV[1]|g' Run_TGFam-Finder.sh");
	system("perl -pi -e 's|CONF3|$ARGV[2]|g' Run_TGFam-Finder.sh");

	system("mkdir -p $RUNNING_PATH/$AUGUSTUS_ANALYSIS_PATH");
	system("mkdir -p $RUNNING_PATH/$ISGAP_ANALYSIS_PATH");
	system("mkdir -p $RUNNING_PATH/$PROTEIN_MAPPING_ANALYSIS_PATH");
	system("mkdir -p $RUNNING_PATH/$RNASEQ_AUGUSTUS_ANALYSIS_PATH");
	system("mkdir -p $RUNNING_PATH/$PM_AUGUSTUS_ANALYSIS_PATH");
	system("mkdir -p $RUNNING_PATH/$ISGAP_AUGUSTUS_ANALYSIS_PATH");
	system("mkdir -p $RUNNING_PATH/$MERGING_ANALYSIS_PATH");
	system("mkdir -p $RUNNING_PATH/$FINAL_ANALYSIS_PATH");
 
	print "* The result folders of annotation in 'PROGRAM_PATH' has been created. \n";
	print "* \n";
	print "************************************************************************************ \n\n";

    print STDERR "Type '. Run_TGFam-Finder.sh' for run TGFam-Finder.\n";
  }
