#!/usr/bin/perl
use strict;

if (@ARGV !=1 || $ARGV[0] eq "-h") {
	printUsage();
	}

	print "Do you want full-installation (Y:full-installation, N:Excluding InterProScan)? (Y/N)\n";
	my $rc = <STDIN>;
	$rc =~ s/\s+//g;
        if ($rc =~ /N/i){
		system("cp -rf auto_install_tgfam_noiprscan tgfam_install_centos.sh");
	} else {
		system("cp -rf auto_install_tgfam_all tgfam_install_centos.sh");
	}

my $stTGFamSoftLoc = $ARGV[0]; ### TGFam Program Installation Path

system("perl -pi -e 's|setTGFamSoftLoc|$stTGFamSoftLoc|g' tgfam_install_centos.sh");

system("cp -rf ./test_sample/PROGRAM_PATH_sample.config ./test_sample/PROGRAM_PATH_test.config");
system("perl -pi -e 's|setTGFamSoftLoc|$stTGFamSoftLoc|g' ./test_sample/PROGRAM_PATH_test.config");

system("cp -rf ./test_sample/RESOURCE_sample.config ./test_sample/RESOURCE_test.config");
system("perl -pi -e 's|setTGFamSoftLoc|$stTGFamSoftLoc|g' ./test_sample/RESOURCE_test.config");

system("mkdir $stTGFamSoftLoc");

system("cp -rf ./test_sample/PROGRAM_PATH_test.config $stTGFamSoftLoc/PROGRAM_PATH.config");
system("cp -rf ./RESOURCE.config $stTGFamSoftLoc/RESOURCE.config");
system("cp -rf ./scripts $stTGFamSoftLoc");
system("cp -rf ./test_sample $stTGFamSoftLoc");
system("rm -rf $stTGFamSoftLoc/test_sample/PROGRAM_PATH_sample.config");
system("rm -rf $stTGFamSoftLoc/test_sample/RESOURCE_sample.config");
system("rm -rf ./test_sample/PROGRAM_PATH_test.config");
system("rm -rf ./test_sample/RESOURCE_test.config");
system("mv tgfam_install_centos.sh $stTGFamSoftLoc/");
system(". $stTGFamSoftLoc/tgfam_install_centos.sh");

sub printUsage {
	print "TGFam-Finder Installation for CentOS User\n";
	print "Help : perl tgfam_install.pl -h\n";
	print "Usage : perl tgfam_install.pl [Absolute Path for TGFam-Finder_Install_Location]\n";
	print "    ex) perl tgfam_install.pl /home/TGFam\n";
	exit(0);
}
