# Scripts_CentOS
TGFam-Finder scripts (CentOS)

## 1. Introduction ##
Target Gene Family Finder (TGFam-Finder) is an annotation tool that works in the Linux OS environment for structural annotation of protein-coding genes containing target domains of interest. TGFam-Finder includes scripts that allow starting from automatic installation to generation of gene models for target-gene family of interest. Through auto-installation process, users can automatically setup prerequisite tools for further annotation using TGFam-Finder. Users also can use TGFam-Finder scripts for annotation without auto-installation process after self-installation of the prerequisite tools. 

Prior to operation of TGFam-Finder, users should input location information of genomic resources in a ‘RESOURCE.config’ file and the prerequisite tools in a ‘PROGRAM_PATH.config’ file. The TGFam-Finder package automatically generates ‘RESOURCE_PATH.config’ for provided sample data and ‘PROGRAM_PATH.config’ for the installed tools. This manual provides instructions for installation processes, preparation of configuration information, annotation command of TGFam-Finder, and detailed information of the prerequisite tools.

## 2. TGFam-Finder installation ##
Prior to installation of TGFam-Finder, users should complete installation of following programs, **Perl 5.6.1 or higher version, Bioperl, and YAML (perl module)**. Then, an auto-installation script in TGFam-Finder performs installation of the other prerequisite tools, including Bowtie2-2.3.1, HMMER-3.1b2, BLAST 2.6.0+, InterproScan-5.22-61.0, Exonerate-2.2.0, Blat v35, Tophat-2.1.1 and Cufflinks-2.2.1, Augustus-3.2.3, Scipio-1.4, and ClustalW-2.1, with yum (apt-get) update to perform annotation using TGFam-Finder. If users already updated **yum (apt-get)**, the auto-installation process can be performed using general user account as well as root authority.

> 2.1 System requirements for installation of TGFam-Finder

- Hard drive capacity: over **20 Gb** storage
- OS: Fedora **CentOS version 6** or higher (**Ubuntu version 16.04** or higher)
- Prerequisite softwares: 1) Perl 5.6.1 or higher version, Bioperl, and YAML (perl module), 2) Yum (Ubuntu: apt-get), and 3) Bowtie2-2.3.1, HMMER-3.1b2, BLAST 2.6.0+, InterproScan-5.22-61.0, Exonerate-2.2.0, Blat v35, Tophat-2.1.1 and Cufflinks-2.2.1, Augustus-3.2.3, Scipio-1.4, and ClustalW-2.1
- Prerequisite softwares in the first category should be installed before installation of the tools in the second and the third categories. TGFam-Finder provides a script that allows automatic installation of tools in the second and the third categories. If users cannot access root account, Yum (or apt-get) should be installed or updated as the latest version before the auto-installation process. Then, users can perform the auto-installation using any account

> 2.2 Download TGFam-Finder script

- You can download TGFam-Finder script in TGFam-Finder web site (http://tgfam-finder.snu.ac.kr) or Github (https://github.com/tgfam-finder).

	**Usage:** 

		$wget http://tgfam-finder.snu.ac.kr/download/TGFam-Finder_'latest'.tar.gz
		ex) $wget http://tgfam-finder.snu.ac.kr/download/TGFam-Finder_CentOS_v1.01.tar.gz

	Github download. 

		CentOS (Fedora, Redhat series)
		$git clone https://github.com/tgfam-finder/scripts_centos.git
		
		Ubuntu Series (Debian series)
		$git clone https://github.com/tgfam-finder/scripts_ubuntu.git

	
> 2.3 Extraction of the downloaded file

- **Usage:**

		$tar zxvf file name
		ex) tar zxvf TGFam-Finder_CentOS_v1.01.tar.gz


> 2.4 Installation



- Basically, the installation script performs full auto-installation of the prerequisite tools that we mentioned in the section **“System requirements for installation of TGFam-Finder”**. Considering huge size of InterproScan causing download delay, users can do partial installation except for installation of InterproScan v5. We recommend that if users are hard to download InterproScan through our package, users should use the install script after direct download from the website (https://www.ebi.ac.uk/interpro/download.html) and self-installation of InterproScan. For users who want to manually setup the prerequisite tools, users can only use the annotation scripts. If users consider manual installation of all prerequisite tools, see the section 5 **“Information of the used prerequisite tools”**.

	**Usage:** 

		$perl install_TGFam-Finder.pl absolute_installation_path
		(ex. $perl install_TGFam-Finder.pl /home/TGFam-Finder/TGFam)

	If you want to confirm completion of installation, you need to check list of programs that installed in the installation path.


## 3. Registration for location of programs and resources
To run TGFam-Finder, users should prepare RESOURCE.config including full location of genomic resources, and PROGRAM_PATH.config containing absolute location of pre-installed programs and name of output directories. Through the full auto-installation, RESOURCE.config for sample data and PROGRAM_PATH.config are automatically generated. Users who do not want to use our installation script should manually input location of the programs in PROGRAM_PATH.config as well as the resources in RESOURCE.config. If users perform partial-installation except for InterproScan using our script, the location of InterproScan should be inserted in PROGRAM_PATH.config.

> 3.1. PROGRAM_PATH.config

    ## DEFAULT PATH ##
    1. $TGFAM_SCRIPTS_PATH="";  ## Directory path containing annotation scripts
    2. $RUNNING_PATH="";  ## Location of directory running TGFam-Finder 

    ## PROGRAM PATH ##
    1. $CLUSTALW_PATH = "";  ## Directory path containing ClustalW execution file
    2. $HMMER_BIN_PATH = "";  ## Directory path containing Hmmer execution file
    3. $BLAST_BIN_PATH = "";  ## Directory path containing BLAST execution file
    4. $IPRSCAN_PATH = "";  ## Directory path including InterproScan execution file 
    5. $BOWTIE_PATH = "";  ## Directory path including Bowtie2 execution file
    6. $TOPHAT_PATH = "";  ## Directory path including Tophat execution file
    7. $EXONERATE_BIN_PATH = "";  ## Directory path containing exonerate execution file
    8. $AUGUSTUS_PATH = "";  ## Directory path of augustus containing ‘bin’ directory
    9. $AUGUSTUS_BIN_PATH = "";  ## Directory path including augustus execution file
    10. $AUGUSTUS_CONFIG_PATH = "";  ## Location of config directory of augustus
    11. $SCIPIO_PATH = "";  ## Directory path containing scipio execution file 
    12. $BLAT_PATH = "";  ## Directory path including blat execution file
    13. $CUFFLINKS_PATH = "";  ## Directory path including cufflinks execution file
    
	## OUTPUT PATH OF ANALYSIS RESULTS ##
	1. $AUGUSTUS_ANALYSIS_PATH="";  ## Augustus annotation output directory name that will be created in RUNNING_PATH ex) Augustus
	2. $ISGAP_ANALYSIS_PATH="";  ## ISGAP output directory name
	3. $MERGING_ANALYSIS_PATH="";  ## Output directory name for merging initial gene models
	4. $PROTEIN_MAPPING_ANALYSIS_PATH="";  ## Protein mapping output directory name
	5. $RNASEQ_AUGUSTUS_ANALYSIS_PATH="";  ## Output directory name of augustus annotation for assembled transcripts
	6. $PM_AUGUSTUS_ANALYSIS_PATH="";  ## Output directory name for comparison between protein mapping and augustus
	7. $ISGAP_ AUGUSTUS_ANALYSIS_PATH="";  ## Output directory for comparison of results between ISGAP and augustus annotations
	8. $FINAL_ANALYSIS_PATH="";  ## Final result output directory name

> 3.2. RESOURCE.config

    ## REQUIRED INFORMATION ##
    1. $TARGET_GENOME = "";  ## Full path of assembled genome (ex, path/filename), fasta format
    2. $PROTEINS_FOR_DOMAIN_IDENTIFICATION = "";  ## Full path of peptide sequences of target or allied species, fasta format
    3. $TSV_FOR_DOMAIN_IDENTIFICATION = "";  ## Full path of InterPro result of the peptide sequences, tsv format
    4. $RESOURCE_PROTEIN = "";  ## Full path of target peptide sequences in multiple species that users prepare as resources for each annotation step, fasta format
    5. $BLAST_DB_NAME = "";  ## Full path of BLAST DB that will be automatically created
    6. $OUTPUT_PREFIX = "";  ## Output prefix generated results
    7. $TARGET_DOMAIN_ID = "";  ## Target domain ID saved in fifth column of tsv files such as PF00000, SSF00000, and SM00000, not IPR number
    8. $TARGET_DOMAIN_NAME = "";  ## target domain name(s)
    9. $REPRESENTATIVE_DOMAIN_NAME = "";  ## Target gene-family name
    10. $EXTENSION_LENGTH = "";  ## Extension length for flanking regions of target domain
    11. $MAX_INTRON_LENGTH = "";  ## Max intron length
    12. $THREADS = "";  ## Number of threads

    ## OPTIONAL INFORMATION ##
    1. $CDS_OF_TARGET_GENOME = "";  ## Full path of coding DNA sequences of the target genome, fasta format
    2. $GFF3_OF_TARGET_GENOME = "";  ## Full path of GFF3 of the target genome, gff3 
    3. $RNASEQ_FORWARD_PATH = "";  ## Full path of RNA-seq (forward), fastaq format
    4. $RNASEQ_REVERSE_PATH = "";  ## Full path of RNA-seq (reverse), fastaq format
    5. $EXCLUDED_DOMAIN_ID = "";  ## Excluded Target domain ID
    

**Note.** If you have problems for preparation of config files, see the section 6.2 and 6.3.

## 4. Operation command of TGFam-Finder ##

After the completion of configuration for PROGRAM_PATH.config and RESOURCE.config, users can annotate target-gene family of interest in an assembled genome. 

> 4.1. Preparation for running TGFam-Finder

- 4.1.1. Moving to running directory where users want to save annotation results (RUNNING_PATH in PROGRAM_PATH.config)

	**Usage:**

	    $cd folder_name_saving_output_of_annotation.    
	    (ex. $cd /home/TGFam-Finder/TGFam/test_sample)

- 4.1.2. Construction of output directories in the running folder

	Build_TGFam-Finder.pl makes output directories in the running folder and check whether a directory has the same name of output folder in PROGRAM_PATH . If there is a directory with the same name of output folder (it gives an error message). In this case delete or rename the directory. It also creates Run_TGFam-Finder.sh, a script that runs TGFam-Finder.

	**Usage:**

	    $perl Build_TGFam-Finder.pl absolute_path_containing_config_files RESOURCE.config PROGRAM_PATH.config
   		(ex. $perl /home/TGFam-Finder/TGFam/scripts/Build_TGFam-Finder.pl /home/TGFam-Finder/TGFam/test_sample PROGRAM_PATH_test.config RESOURCE_test.config)

> 4.2. Run TGFam-Finder

- Basically, TGFam-Finder performs full annotation of target-gene family of interest. Note that you must use the .(dot) when using the shell command.

	**Usage:**

		$. Run_TGFam-Finder.sh
		(ex. $. Run_TGFam-Finder.sh)

> 4.3. Final outputs of TGFam-Finder

- Final gene models including gff3, peptides, coding DNA sequences and tsv files are generated

	**Usage:**

		$cd Final_directory && ls
		(ex. $cd /home/TGFam-Finder/TGFam/test_sample/Final && ls)

	**Usage:**
	
		$cat Prefix.TGFam.CDS.fa

	**Usage:**

		$cat Prefix.TGFam.gff3

	**Usage:**

		$cat Prefix.TGFam.PEP.fa

	**Usage:**

		$cat Prefix.TGFam.PEP.fa.nostop.tsv

## 5. Information of the used prerequisite tools

TGFam-Finder provides an automatic install script, but if you want to use TGFam-Finder through already installed programs, you should install the following programs.

> **List of programs used in TGFam-Finder**


- If you do not want Auto Installation, or if you are using the self-installer, you need to install the program below to use TGFam-Finder.

	**Bioperl**

	Homepage: https://bioperl.org/index.html

	Manual: https://bioperl.org/INSTALL.html

	
	**Bowtie2**

	Homepage: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

	Manual: http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml

	TGFam-Finder Current Version: 2.3.1-x86_64

	
	**Samtools**

	Homepage: http://www.htslib.org/ 

	Manual: http://www.htslib.org/doc/#manual-pages

	TGFam-Finder Current Version: 0.1.19-x86_64

	**Hmmer**

	Homepage: http://hmmer.org/ 

	Manual: http://hmmer.org/documentation.html

	TGFam-Finder Current Version: 3.1b2-x86_64

	
	**BLAST+** 

	Homepage: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download 

	Manual: https://www.ncbi.nlm.nih.gov/books/NBK279690/
	TGFam-Finder Current Version: 2.6.0+-x86_64
	
	**Exonerate** 

	Homepage: https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate
 
	Manual: https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate-user-guide
	TGFam-Finder Current Version: 2.2.0-x86_64
	
	**Tophat** 

	Homepage: https://ccb.jhu.edu/software/tophat/index.shtml 

	Manual: https://ccb.jhu.edu/software/tophat/manual.shtml
	TGFam-Finder Current Version: 2.1.1-x86_64
	
	**Cufflinks**

	Homepage: http://cole-trapnell-lab.github.io/cufflinks/ 

	Manual: https://github.com/cole-trapnell-lab/cufflinks

	TGFam-Finder Current Version: 2.2.1-x86_64

	
	**Scipio** 

	Homepage: http://www.webscipio.org/ 

	Manual: https://www.webscipio.org/webscipio/download_scipio

	TGFam-Finder Current Version: scipio-1.4

	
	**Augustus**

	Homepage: http://bioinf.uni-greifswald.de/augustus/ 

	Manual: http://bioinf.uni-greifswald.de/augustus/

	TGFam-Finder Current Version: augustus-3.2.3

	
	**Blat** 

	Homepage: http://genome.ucsc.edu/cgi-bin/hgBlat?command=start 

	Manual: http://genome.ucsc.edu/goldenPath/help/blatSpec.html

	TGFam-Finder Current Version: 35

	
	**ClustalW** 

	Homepage: http://www.clustal.org/clustal2/ 

	Manual: http://www.clustal.org/clustal2/#Documentation

	TGFam-Finder Current Version: 2.1
	
	**InterProScan**

	Homepage: https://www.ebi.ac.uk/interpro/

	Manual: https://www.ebi.ac.uk/interpro/documentation.html

	TGFam-Finder Current Version: 5.22-61-64bit


## 6. Points to note ##
> 6.1 If TGFam-Finder installation is failed

- TGFam-Finder is a program based on known bioinformatics programs. Auto installation script is provided to make it easier for the user to install, since the specifications and status of all servers or computers are not the same, sometimes the installation is not perfect and errors may occur.
 
	At this time, check whether the following libraries and modules are properly installed, and then execute the automatic installation again.


	**Check list**

	Yum (Ubuntu: apt-get): Use yum to update the following packages
	- Packages for rpm or yum-based Linux distributions (RedHat / Fedora / CentOS) are:
	zlib-devel, bzip2-devel, xz-devel
	- Packages for dpkg-based Linux distributions (Debian / Ubuntu) are:
	zlib1g-dev, libbz2-dev, liblzma-dev
	
	Perl 5.6.1 or higher Version 5.8 or higher is highly recommended. Modules are tested against version 5.8 and above on 5.8 and above.

	BioPerl, perl module – yaml, Glib version 3.20 or higher, KentLib, Libpng, Ncurses, 
	GNU make & C compiler (e.g. gcc or clang), Boost C++ libraries (version 1.47 or higher)


> 6.2 Configuration 

- Users can input multiple ‘TARGET_DOMAIN_ID’, ‘TARGET_DOMAIN_NAME’ using comma delimiter but target domain names should be distinct.

		ex) $TARGET_DOMAIN_ID = “PF00319,cd00265”
	        $TARGET_DOMAIN_NAME = “SRF,MEF”

- If users want to obtain final gene model including existing gene model of the target genome, users need to input CDS_OF_TARGET_GENOME and GFF3_OF_TARGET_GENOME. 
- IDs in peptides (not RESOURCE_PROTEIN), tsv, coding DNA sequences, and gff3 should be matched. If users use merged peptides and tsv of genomes including target and allied species, IDs of the peptides and tsv should contain IDs of CDS and gff3. 
- If users don’t insert RNASeq information, ISGAP and related analyses will not be performed and final gene model is generated combining protein mapping and augustus annotation except for ISGAP and related analyses.

> 6.3 Cautions for preparation of genomic resources

- Gff3 file should be simplified. Like the capture below, information of coding genes in the third column should contain gene, mRNA, and CDS. The last column only includes ID name like ’ID=XP00000’ (see below).

- RESOURCE_PROTEIN in RESOURCE.config means combined target protein sequences in multiple genomes and users should prepare them. For example, we extracted and used protein sequences as RESOURCE PROTEIN containing FAR1 (PF3101) NBARC (PF00931) in plant genomes as well as C2H2 zinc finger(PF00096) and homeobox(PF00046) in animal genomes.
The IDs in protein and tsv files of target or allied species should be less than 25 characters.
