#!/usr/bin/env perl

use FindBin;
use lib $FindBin::RealBin;

use Plgd::Utils;
use Plgd::Script;
use Plgd::Project;

use strict;

sub defaultConfig() {
    return (
        PROJECT=>"",
        RAWREADS=>"",
        GENOME_SIZE=>"",
        THREADS=>4,
        MIN_READ_LENGTH=>2000,
        CNS_OVLP_OPTIONS=>"-kmer_size 13",
        CNS_PCAN_OPTIONS=>"-p 100000 -k 100",
        CNS_OPTIONS=>"",
        CNS_OUTPUT_COVERAGE=>30,
        TRIM_OVLP_OPTIONS=>"-skip_overhang",
        TRIM_PM4_OPTIONS=>"-p 100000 -k 100",
        TRIM_LCR_OPTIONS=>"",
        TRIM_SR_OPTIONS=>"",
        ASM_OVLP_OPTIONS=>"",
        CLEANUP=>0,
        FSA_OL_FILTER_OPTIONS=>"--max_overhang=-1 --min_identity=-1",
        FSA_ASSEMBLE_OPTIONS=>"",
    );
}

sub loadMecatConfig($) {
    my ($fname) = @_;
    my %cfg = defaultConfig();

    loadConfig($fname, \%cfg);

    my @required = ("PROJECT", "GENOME_SIZE", "RAWREADS");
    foreach my $r (@required) {
        if (not exists($cfg{$r}) or $cfg{$r} eq "")  {
            plgdError("Not set config $r");
        }
    }


    return %cfg;
}

sub loadMecatEnv($) {
    my ($cfg) = @_;

    my %env = loadEnv($cfg);
    $env{"BinPath"} = $FindBin::RealBin;
    return %env;    
}

sub initializeMecatProject($) {
    my ($cfg) = @_;
    initializeProject($cfg);
}

sub runCorrectRawreads($$) {
    my ($env, $cfg) = @_;

    my $prjDir = %$env{"WorkPath"} ."/". %$cfg{"PROJECT"};
    my $workDir = "$prjDir/1-consensus";
    mkdir $workDir;
    my $rawreads = %$cfg{"RAWREADS"};
    
    my $thread = %$cfg{"THREADS"};
    my $genomeSize = %$cfg{"GENOME_SIZE"};
    my $minReadSize = %$cfg{"MIN_READ_LENGTH"};
    my $coverage = %$cfg{"CNS_OUTPUT_COVERAGE"};
    my $binPath = %$env{"BinPath"};
    my $cnsOvlpOptions = %$cfg{'CNS_OVLP_OPTIONS'};
    my $cnsPcanOptions = %$cfg{'CNS_PCAN_OPTIONS'};
    my $cnsOptions = %$cfg{'CNS_OPTIONS'};
    
    my $jobPw = Job->new(
        name => "cns_pw",
        ifiles => [$rawreads],
        ofiles => ["$workDir/cns_pm.seqidx"],
        gfiles => ["$workDir/cns_pm*"],
        mfiles => [],
        cmds => ["$binPath/mecat2map $cnsOvlpOptions -task pm -outfmt seqidx -num_threads $thread -db_dir $workDir/cns_pm_dir -keep_db -min_query_size $minReadSize -out $workDir/cns_pm.seqidx $rawreads $rawreads"],
        #cmds => ["$binPath/mecat2pw -j 0 -d $rawreads -o $workDir/cns_pm.can -w $workDir/cns_pm_dir -t $thread $cnsOvlpOptions"],
        msg => "correcting rawreads step 1 mecat2map",
    );

    my $jobPcan = Job->new(
        name => "cns_pcan",
        ifiles => ["$workDir/cns_pm.seqidx"],
        ofiles => [],
        gfiles => [],
        mfiles => [],
        cmds   => ["$binPath/mecat2pcan $cnsPcanOptions -t $thread $workDir/cns_pm_dir $workDir/cns_cns_dir $workDir/cns_pm.seqidx"],
        msg    => "partition correction candidates step 2 mecat2pcan",
    );

    my $jobCns = Job->new(
        name => "cns_cns",
        ifiles => [],
        ofiles => [],
        gfiles => [],
        mfiles => [],
        cmds   => ["$binPath/mecat2cns $cnsOptions -num_threads $thread $workDir/cns_pm_dir $workDir/cns_cns_dir"],
        #cmds => ["$binPath/mecat2cns -i 0 -t $thread $cnsOptions $workDir/cns_pm.can $rawreads $workDir/cns_reads.fasta"],
        msg => "correcting rawreads step 3 mecat2cns",
    );

    my $jobMakeCnsReadList = Job->new(
        name   => "cns_make_list",
        ifiles => [],
        ofiles => [],
        gfiles => [],
        mfiles => [],
        cmds   => ["ls $workDir/cns_cns_dir/p*.cns.fasta > $workDir/cns_reads_list.txt"],
        msg    => "correcting rawreads step 4 make consensus reads list",
    );

    my $jobExtr = Job->new(
        name => "cns_extract",
        ifiles => [],
        ofiles => [],
        gfiles => [],
        mfiles => [],
        #cmds => ["$binPath/extract_sequences $workDir/cns_reads.fasta $workDir/cns_final $genomeSize $coverage"],
        #cmds => ["$binPath/mecat2elr $workDir/cns_reads.fasta $genomeSize $coverage $workDir/cns_final.fasta"],
        cmds   => ["$binPath/mecat2extseqs $genomeSize $coverage $workDir/cns_reads_list.txt > $workDir/cns_final.fasta"],
        msg => "correcting rawreads step 5 extract_sequences",
    );

    my $job = Job->new (
        name => "cns_job",
        ifiles => [$rawreads],
        ofiles => ["$workDir/cns_final.fasta"],
        mfiles => [],
        jobs => [$jobPw, $jobPcan, $jobCns, $jobMakeCnsReadList, $jobExtr],        
        msg => "correcting rawreads",
    );
    
    serialRunJobs($env, $cfg, $job);
}


sub runTrimReads($$) {
    my ($env, $cfg) = @_;

    my $prjDir = %$env{"WorkPath"} . "/" .%$cfg{"PROJECT"};
    my $workDir = "$prjDir/2-trim_bases";
    mkdir $workDir;
    my $volDir = "$workDir/trim_pm_dir";
    mkdir $volDir;
    my $pm4Dir = "$workDir/trim_pm4_dir";
    mkdir -p $pm4Dir;

    my $cnsReads = "$prjDir/1-consensus/cns_final.fasta";
    my $trimReads = "$prjDir/2-trim_bases/trimReads.fasta"; 
    my $trimPm = "$volDir/trim_pm.m4x";
    my $binPath = %$env{"BinPath"};
    my $pmOptions = %$cfg{"TRIM_OVLP_OPTIONS"};
    my $pm4Options = %$cfg{"TRIM_PM4_OPTIONS"};
    my $lcrOptions = %$cfg{"TRIM_LCR_OPTIONS"};
    my $srOptions = %$cfg{"TRIM_SR_OPTIONS"};
    my $thread = %$cfg{"THREADS"};

    my $jobPm = Job->new (
        name   => "tr_pm",
        ifiles => [],
        ofiles => [],
        gfiles => [],
        mfiles => [],
        cmds   => ["$binPath/mecat2map $pmOptions -num_threads $thread -db_dir $volDir -keep_db -task pm -outfmt m4x -out $trimPm $cnsReads $cnsReads"],
        msg    => "pairwise mapping for trimming",
    );

    my $lcrResult = "$workDir/lcr.txt";
    my $srResult = "$workDir/sr.txt";

    my $jobPm4 = Job->new (
        name   => "tr_pm4",
        ifiles => [],
        ofiles => [],
        gfiles => [],
        mfiles => [],
        cmds   => ["$binPath/mecat2pm4 $pm4Options -t $thread $volDir $pm4Dir $trimPm"],
        msg    => "partition m4 records for trimming",        
    );

    my $jobLcr = Job->new (
        name   => "tr_lcr",
        ifiles => [],
        ofiles => [],
        gfiles => [],
        mfiles => [],
        cmds   => ["$binPath/mecat2lcr $lcrOptions -num_threads $thread -out $lcrResult $volDir $pm4Dir"],
        msg    => "find largest cover range for trimming",
    );

    my $jobSr = Job->new (
        name   => "tr_sr",
        ifiles => [],
        ofiles => [],
        gfiles => [],
        mfiles => [],
        cmds   => ["$binPath/mecat2splitreads $srOptions -num_threads $thread -out $srResult $volDir $pm4Dir $lcrResult"],
        msg    => "find largest clear range for trimming",
    );

    my $jobTb = Job->new (
        name   => "tr_tb",
        ifiles => [],
        ofiles => [],
        gfiles => [],
        mfiles => [],
        cmds   => ["$binPath/mecat2trimbases $volDir $srResult 1 $trimReads"],
        msg    => "trim sequences to their largest clear ranges",
    );

    my $job = Job->new (
        name => "tr_job",
        ifiles => [$cnsReads],
        ofiles => [$trimReads],
        mfiles => [],
        jobs => [$jobPm, $jobPm4, $jobLcr, $jobSr, $jobTb],
        msg => "trimming corrected reads",
    );
    
    serialRunJobs($env, $cfg, $job);
}

sub runAlignTReads($$) {
    my ($env, $cfg) = @_;

    my $prjDir = %$env{"WorkPath"} ."/". %$cfg{"PROJECT"};
    my $workDir = "$prjDir/3-assembly";
    mkdir $workDir;

    my $volDir = "$workDir/asm_pm_dir";
    mkdir $volDir;

    my $binPath = %$env{"BinPath"};
    my $trimReads = "$prjDir/2-trim_bases/trimReads.fasta";

    my $asmPm = "$workDir/asm_pm.m4";
    my $options = %$cfg{"ASM_OVLP_OPTIONS"};
    my $thread = %$cfg{"THREADS"};

    my $jobPm = Job->new (
        name   => "altr_pm",
        ifiles => [$trimReads],
        ofiles => [],
        gfiles => [],
        mfiles => [],
        cmds => ["$binPath/mecat2map $options -task pm -num_threads $thread -db_dir $volDir -out $asmPm $trimReads $trimReads"],
        msg => "pairwise mapping of trimmed reads",
    );
    
    my $job = Job->new (
        name => "altr_job",
        ifiles => [$trimReads],
        ofiles => [$asmPm],
        mfiles => [$volDir],
        jobs => [$jobPm],
        msg => "aligning trimmed reads for assembling",
    );
    
    serialRunJobs($env, $cfg, $job);
}

sub runAssemble($$) {
    my ($env, $cfg) = @_;

    my $prjDir = %$env{"WorkPath"} . "/" .%$cfg{"PROJECT"};
    my $workDir = "$prjDir/4-fsa";
    mkdir $workDir;

    my $script = "$prjDir/scripts/assemble.sh";
    my $overlaps = "$prjDir/3-assembly/asm_pm.m4";
    my $reads = "$prjDir/2-trim_bases/trimReads.fasta";
    my $contigs = "$workDir/contigs.fasta";
    my $filtered_overlaps = "$workDir/filter.m4";

    my $binPath = %$env{"BinPath"}; 
    my $thread = %$cfg{"THREADS"};
    my $filterOptions = %$cfg{"FSA_OL_FILTER_OPTIONS"};
    if (%$cfg{"GENOME_SIZE"}) {
        $filterOptions = $filterOptions . " --genome_size=" . %$cfg{"GENOME_SIZE"};
    }
    my $assembleOptions = %$cfg{"FSA_ASSEMBLE_OPTIONS"};

    my $job = Job->new(
        name => "ass_job",
        ifiles => [$overlaps, $reads],
        ofiles => [$filtered_overlaps, $contigs],
        gfiles => [$filtered_overlaps, $contigs],
        mfiles => [],
        cmds => ["$binPath/fsa_ol_filter $overlaps $filtered_overlaps --thread_size=$thread --output_directory=$workDir $filterOptions", 
                 "$binPath/fsa_assemble $filtered_overlaps --read_file=$reads --thread_size=$thread --output_directory=$workDir $assembleOptions"],
        msg => "assembling",
    );

    serialRunJobs($env, $cfg, $job);
}


sub statCorrectedReads($$) {
    my ($env, $cfg) = @_;
    my $prjDir = %$env{"WorkPath"} . "/" .%$cfg{"PROJECT"};
 
    plgdInfo("Information of corrected reads $prjDir/1-consensus/cns_final.fasta");
    my $cmd = %$env{"BinPath"} . "/mecat2viewdb $prjDir/1-consensus/cns_final.fasta";
    system($cmd);
}

sub statTrimmedReads($$) {
    my ($env, $cfg) = @_;
    my $prjDir = %$env{"WorkPath"} . "/" .%$cfg{"PROJECT"};
 
    plgdInfo("Information of trimmed reads $prjDir/2-trim_bases/trimReads.fasta");
    my $cmd = %$env{"BinPath"} . "/mecat2viewdb $prjDir/2-trim_bases/trimReads.fasta";
    system($cmd);
}

sub statContigs($$) {
    my ($env, $cfg) = @_;
    my $prjDir = %$env{"WorkPath"} . "/" .%$cfg{"PROJECT"};
 
    my $cmd = %$env{"BinPath"} . "/mecat2viewdb $prjDir/4-fsa/contigs.fasta";
    plgdInfo("N50 of contigs: $prjDir/4-fsa/contigs.fasta");
    system($cmd);
}

my %cfg = ();
my %env = ();

sub cmdCorrect($) {
    my ($fname) = @_;

    %cfg = loadMecatConfig($fname);
    %env = loadMecatEnv(\%cfg);
    initializeMecatProject(\%cfg);

    runCorrectRawreads(\%env, \%cfg);
    statCorrectedReads(\%env, \%cfg);
}

sub cmdTrim($) {
    
    my ($fname) = @_;

    %cfg = loadMecatConfig($fname);
    %env = loadMecatEnv(\%cfg);
    initializeMecatProject(\%cfg);

    runCorrectRawreads(\%env, \%cfg);
    runTrimReads(\%env, \%cfg); 
    statTrimmedReads(\%env, \%cfg);
}

sub cmdAssemble($) {
    
    my ($fname) = @_;

    %cfg = loadMecatConfig($fname);
    %env = loadMecatEnv(\%cfg);
    initializeMecatProject(\%cfg);

    runCorrectRawreads(\%env, \%cfg);
    runTrimReads(\%env, \%cfg);
    runAlignTReads(\%env, \%cfg);
    runAssemble(\%env, \%cfg); 
    statContigs(\%env, \%cfg);
}

sub cmdConfig($) {
    my ($fname) = @_;

    my %cfg = defaultConfig();

    my @items = ("PROJECT", "RAWREADS", "GENOME_SIZE", "THREADS", "MIN_READ_LENGTH", 
                 "CNS_OVLP_OPTIONS", "CNS_PCAN_OPTIONS", "CNS_OPTIONS", "CNS_OUTPUT_COVERAGE", 
                 "TRIM_OVLP_OPTIONS", "TRIM_PM4_OPTIONS", "TRIM_LCR_OPTIONS", "TRIM_SR_OPTIONS",
                 "ASM_OVLP_OPTIONS", "FSA_OL_FILTER_OPTIONS", "FSA_ASSEMBLE_OPTIONS");

    open(F, "> $fname") or die; 
    foreach my $k (@items) {
        if ($k =~ /OPTIONS/) {
            print F "$k=\"$cfg{$k}\"\n"
        } else {
            print F "$k=$cfg{$k}\n";
        }
    }
    
    foreach my $k (keys %cfg) {
        if (not grep /^$k$/, @items) {
            print F "$k=$cfg{$k}\n";
        }
    }
    close(F);

}



sub usage() {
    print "Usage: mecat.pl correct|assemble|config cfg_fname\n".
          "    correct:     correct rawreads\n" .
          "    trim:        trim corrected reads\n" .
          "    assemble:    generate contigs\n" .
          "    config:      generate default config file\n"
}

sub main() {
    if (scalar @ARGV >= 2) {
        my $cmd = @ARGV[0];
        my $cfgfname = @ARGV[1];

        if ($cmd eq "correct") {
            cmdCorrect($cfgfname);
        } elsif ($cmd eq "trim") {
            cmdTrim($cfgfname);
        } elsif ($cmd eq "assemble") {
            cmdAssemble($cfgfname);
        } elsif ($cmd eq "config") {
            cmdConfig($cfgfname);
        } else {
            usage();
        }
    } else {
        usage();
    }
}


$SIG{TERM}=$SIG{INT}=\& catchException;
sub catchException { 
    plgdInfo("Catch an Exception, and do cleanup");
    stopRunningScripts(\%env, \%cfg);
    exit -1; 
} 

#eval {
    main();
#};

if ($@) {
    catchException();
}

END {
    stopRunningScripts(\%env, \%cfg);
}
