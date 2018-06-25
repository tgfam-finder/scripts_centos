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

  if (-d "$RUNNING_PATH/$AUGUSTUS_ANALYSIS_PATH") {
    print "You should be check AUGUSTUS_ANALYSIS_PATH \n";
    exit 1;
  }
  if (-d "$RUNNING_PATH/$ISGAP_ANALYSIS_PATH") {
    print "You should be check ISGAP_ANALYSIS_PATH \n";
    exit 1;
  }
  if (-d "$RUNNING_PATH/$PROTEIN_MAPPING_ANALYSIS_PATH") {
    print "You should be check PROTEIN_MAPPING_ANALYSIS_PATH \n";
    exit 1;
  }
  if (-d "$RUNNING_PATH/$RNASEQ_AUGUSTUS_ANALYSIS_PATH") {
    print "You should be check RNASEQ_AUGUSTUS_ANALYSIS_PATH \n";
    exit 1;
  }
  if (-d "$RUNNING_PATH/$PM_AUGUSTUS_ANALYSIS_PATH") {
    print "You should be check PM_AUGUSTUS_ANALYSIS_PATH \n";
    exit 1;
  }
  if (-d "$RUNNING_PATH/$ISGAP_AUGUSTUS_ANALYSIS_PATH") {
    print "You should be check ISGAP_AUGUSTUS_ANALYSIS_PATH \n";
    exit 1;
  }
  if (-d "$RUNNING_PATH/$MERGING_ANALYSIS_PATH") {
    print "You should be check MERGING_ANALYSIS_PATH \n";
    exit 1;
  }
  if (-d "$RUNNING_PATH/$FINAL_ANALYSIS_PATH") {
    print "You should be check FINAL_ANALYSIS_PATH \n";
    exit 1;
  }

  else {
	system("mkdir -p $RUNNING_PATH/$AUGUSTUS_ANALYSIS_PATH");
	system("mkdir -p $RUNNING_PATH/$ISGAP_ANALYSIS_PATH");
	system("mkdir -p $RUNNING_PATH/$PROTEIN_MAPPING_ANALYSIS_PATH");
	system("mkdir -p $RUNNING_PATH/$RNASEQ_AUGUSTUS_ANALYSIS_PATH");
	system("mkdir -p $RUNNING_PATH/$PM_AUGUSTUS_ANALYSIS_PATH");
	system("mkdir -p $RUNNING_PATH/$ISGAP_AUGUSTUS_ANALYSIS_PATH");
	system("mkdir -p $RUNNING_PATH/$MERGING_ANALYSIS_PATH");
	system("mkdir -p $RUNNING_PATH/$FINAL_ANALYSIS_PATH");

    print STDERR "You can Enjoy TGFam!\n";
  }
