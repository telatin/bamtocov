#!/usr/bin/env perl
# ABSTRACT - Test multiple tools against multiple datasets 

use 5.012;
use warnings;
use File::Basename;
use File::Spec;
use Getopt::Long;
use JSON::PP;
use Data::Dumper;

my $out = "./benchmark";
my $minruns = 3;
my $warmup = 1;
# Absolute path of $0
my $abs_script_path = File::Spec->rel2abs(__FILE__);
my $script_dir = dirname($abs_script_path);

GetOptions(
    # Flags
    'help' => \my $help,
    'verbose' => \my $verbose,
    'debug' => \my $debug,
    'r|run'     => \my $run,
    'export' => \my $export,

    # Parameters
    'bam=s' => \my $bam,
    'o|out=s' => \$out,
    'bed=s' => \my $bed,
    'minruns=i' => \$minruns,
    'warmup=i'  => \$warmup,
    'filter=s'    => \my $filter,
    'd|datasets=s' => \my $input_datasets,
    't|tools=s'    => \my $input_tools,
) or die "Error in command line arguments\n";
 

my $datasets = {
    'mini' => {
        'bam' => '{BASE}/../../bamtocov/input/mini.bam',
        'bed' => '{BASE}/../../bamtocov/input/mini.bed',
        'runs' => 10,
    },
    'panel' => {
        'bam' => '{BASE}/../../covtobed/example_data/panel_01.bam',
        'bed' => '{BASE}/../../covtobed/example_data/lstarget.bed',
        'runs' => 5,
    },
    'exome-sub-chr2' => {
        'bam' => '{BASE}/../../exome/HG00258.bam.04-chr2.bam',
        'bed' => '{BASE}/../../exome/target.bed',
        'runs' => 3,
    },
    'exome-sub' => {
        'bam' => '{BASE}/../../exome/HG00258.bam.04.bam',
        'bed' => '{BASE}/../../exome/target.bed',
        'runs' => 3,
    },
    'exome' => {
        'bam' => '{BASE}/../../exome/HG00258.bam',
        'bed' => '{BASE}/../../exome/target.bed',
        'runs' => 2,
    },
    
};
my $tools = {
    "covtobed-1.2.0" => {
        "path" => "/local/miniconda3/bin/covtobed",
        "args" => "{BAM}",
        "check" => "--version",
    },
    "2.0.2-conda"  => {
        "path" => "/local/miniconda3/bin/bamtocov",
        "args" => "{BAM}",
        "check" => "--version" 
    },
    "2.0.4-conda" => {
        "path" => "/local/miniconda3/envs/test-cov/bin/bamtocov",
        "args" => "{BAM}",
        "check" => "--version" 
    },
    "2.1.0-conda" => {
        "path" => "/local/miniconda3/envs/bamtocov210/bin/bamtocov",
        "args" => "{BAM}",
        "check" => "--version" 
    },
    "2.1.0-nimble" => {
        "path" => "{BASE}/../../bamtocov/bin/bamtocov",
        "args" => "{BAM}",
        "check" => "--version" 
    },
    "2.1.0-nimc"   => {
        "path" => "{BASE}/../../bamtocov//test/bamtocov",
        "args" => "{BAM}",
        "check" => "--version" 
    },
    "mosdepth-0.3.1-conda" => {
        "path" => "/local/miniconda3/bin/mosdepth",
        "args" => "/tmp/prefix {BAM}",
        "check" => "--version",
    },  
    "megadepth-1.1.2-conda" => {
        "path" => "/local/miniconda3/envs/mega/bin/megadepth",
        "args" => "--coverage {BAM}",
        "check" => "--version",
    }
};


### Export

if ($export) {
    if (not defined $out) {
        die "Please specify an output prefix with --out PREFIX\n";
    }
    my $json = JSON::PP->new->pretty->canonical;
    open my $datasetsFH, '>', $out . '_dataset.json' or die "Cannot open $out.json: $!\n";
    print $datasetsFH $json->encode($datasets);
    close $datasetsFH;
    open my $toolsFH, '>', $out . '_tools.json' or die "Cannot open $out.json: $!\n";
    print $toolsFH $json->encode($tools);
    close $toolsFH;
    exit;
}


### Import

if (defined $input_tools) { 

    say "Importing datasets from $input_tools";
    my $tools = load_json_file($input_tools);
    my @keys = ('path', 'args', 'check');
    for my $t (keys %$tools) {
        for my $key (@keys) {
            if (not defined $tools->{$t}{$key}) {
                die "ERROR: Malformed JSON *tools* ($input_tools),  $t is missing key $key [@keys]\n";
            }
        }
    }
    
    say "Imported tools.";
  
}

if (defined $input_datasets) { 
    say "Importing datasets from $input_datasets";
    my $datasets = load_json_file($input_datasets);
    my @keys = ('bam', 'runs');
    for my $t (keys %$datasets) {
        for my $key (@keys) {
            if (not defined $datasets->{$t}{$key}) {
                die "ERROR: Malformed JSON *datasets* ($input_datasets),  $t is missing key $key [@keys]\n";
            }
        }
    }

    say "Imported datasets.";
    
}


## Fix paths in datasets and tools
for my $t (keys %$tools) {
    
    $tools->{$t}->{'path'} =~ s/^~/$ENV{HOME}/;
    $tools->{$t}->{'path'} =~ s/\{BASE\}/$script_dir/;
    if ( ! -e $tools->{$t}->{'path'}) {
        die "ERROR: Tool $t does not exist at $tools->{$t}->{'path'}\n";
    }
}
for my $t (keys %$datasets) {
    $datasets->{$t}->{'bam'} =~ s/\{BASE\}/$script_dir/;
    $datasets->{$t}->{'bed'} =~ s/\{BASE\}/$script_dir/;
    if (not defined $datasets->{$t}->{'runs'}) {
        say STDERR Dumper $datasets->{$t}->{'runs'};
        die "ERROR: Dataset $t is missing runs\n";    
    }
    if ( ! -e $datasets->{$t}->{'bam'}) {
        die "ERROR: Dataset $t does not exist at $datasets->{$t}->{'path'}\n";
    }
}

# Print all dataset bams
if ($verbose) {
    say STDERR "DATASETS";
    for my $t (keys %$datasets) {
        say STDERR "Dataset $t: $datasets->{$t}->{'bam'}";
    }
    say STDERR "TOOLS";
    for my $t (keys %$tools) {
        say STDERR "Tool $t: $tools->{$t}->{'path'}";
    }
}
my $hyper = "hyperfine --min-runs {RUNS} --warmup $warmup ";


# Sort $dataset keys by "runs" ascending
my @datasets = sort { $datasets->{$b}->{runs} <=> $datasets->{$a}->{runs} } keys %$datasets;
for my $dset (@datasets) {
        if ( -e "$out-$dset.md" ) {
                say STDERR " Markdown found, skipping $dset";
                next;
            }
        # Hyperfine requires minimum two runs
        my $runs = $datasets->{$dset}->{runs} > 1 ? $datasets->{$dset}->{runs} : 2;
    my @commands;
        if (defined $filter and $dset!~m/$filter/) {
                say "Skipping $dset" if $verbose;
                next;
            }  
        my $file = $datasets->{$dset}->{bam};
        
        say " +++ Processing $dset" if $verbose;
        for my $pkg (sort keys %{ $tools }) {
                my $path = $tools->{$pkg}->{path};
                my $cmd  = $tools->{$pkg}->{check};
                my $args = $tools->{$pkg}->{args};
                if (check_bin($path . " " . $cmd)) {
                            say "Testing $pkg on $dset";
                                
                                my $command = $path . " " . $args;
                                $command =~ s/{BAM}/"$file"/g;
                                push(@commands, "'". $command."'");
                                
                                
                        }
            }
    my $cmd = join(" ", @commands);
my $benchmark_cmd = $hyper . "  --export-markdown $out-$dset.md " . $cmd;
$benchmark_cmd =~ s/{RUNS}/$runs/g;
if ($run) {
    say "\tRunning...";
    system($benchmark_cmd);
    fix_markdown("$out-$dset.md", $tools, $datasets);
}

}


sub fix_markdown {
    my $file = shift;
    my $tools = shift;
    my $datasets = shift;
    my $text = read_file($file);
    for my $t (keys %{ $tools} ) {
        my $path = $tools->{$t}->{path};
        $text =~ s/$path/$t/g;
    }
    for my $d (keys %{ $datasets } ) {
        my $bam = $datasets->{$d}->{bam};
        my $bed = $datasets->{$d}->{bed};
        $text =~ s/$bam/$d.bam/g;
        $text =~ s/$bed/$d.bed/g;
    }
    write_file($file, $text);
}

sub read_file {
    my $file = shift;
    open(my $fh, '<:encoding(UTF-8)', $file)
        or die "Could not open file '$file' $!";
    return <$fh>;
}

sub write_file {
    my $file = shift;
    my $text = shift;
    open(my $fh, '>:encoding(UTF-8)', $file)
        or die "Could not open file '$file' $!";
    print $fh $text;
}

sub countBam {
    my $bam = shift;
    my $cmd = "samtools view -c $bam";
    my $count = `$cmd`;
    chomp($count);
    return $count;
}

sub check_bin {

    my $cmd = shift;
    my $out = `$cmd 2>/dev/null`;
    return 1 if $? == 0;
    return 0;
}

sub load_json_file {
    my $file = shift;
    if ( ! -e $file ) {
        die "ERROR: Cannot find $file\n";
    }
    my $json = JSON::PP->new->pretty->canonical;
    open my $FH, '<', $file or die "Cannot open $file: $!\n";
    my $data = $json->decode(join '', <$FH>);
    close $FH;
    return $data;
}
__DATA__

{BASE}/../../covtobed/example_data/panel_01.bam: 599967
../exome/HG00258.bam: 223754931
input/mini.bam: 25