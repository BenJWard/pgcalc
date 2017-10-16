# pgcalc

*A calculator of population genetic, and evolution statistics, from sequence
files, built on the package ecosystem of BioJulia.*

## Installation

### Compiling from the source code

1. First you must make sure you have an up to date installation of `julia`
   installed on your machine.

2. Download this repository e.g. `git clone https://github.com/Ward9250/pgcalc.git`.

3. `cd` in the downloaded repository, and then run the `make` command.
   `make` will download the BioJulia packages this tool depends on,
   and will compile the pgcalc.jl script into a binary for your machine.

4. When `make` is finished, there will be a folder called `builddir` with two
   files inside, one is a dynamic library file (called `libpgcalc`) needed by
   the executable, and the other file is the executable itself, called `pgcalc`
   (or `pgcalc.exe` on Windows).

## Useage

Using `pgcalc` with the `--help` flag will show you which options are available.

Currently pgcalc only supports sequence (FASTA) files.
As the BioJulia package ecosystem develops, we intend other formats to be
supported in the future.

```sh
N82106:pgcalc bward$ pgcalc --help
usage: <PROGRAM> [-o OUTDIR] [-I SINPUT] [-h] input
                 {sample|pairwise|dndscodon}

commands:
  sample               Compute statistics that apply to the whole
                       sample.
  pairwise             Compute pairwise statistics for your sequence
                       sample.
  dndscodon            Compute pairwise dNdS according to the Nei &
                       Gojoborei method, on a per codon basis.

positional arguments:
  input                A FASTA input file containing sequence
                       alignment.

optional arguments:
  -o, --outdir OUTDIR  Specify an output directory for table files.
                       (default: "./")
  -I, --sinput SINPUT  The FASTA input file containing similated
                       sequence alignments. (default: "")
  -h, --help           show this help message and exit
```

### Computing sample statistics for your sequences

You can compute sample statistics for your sequence file with the `sample` command.

The currently supported statistics are:

* The number of segregating sites. (Use flag `-s` or `--seg`)
* The nucleotide diversity. (Use flag `-n` or `--nucdiv`)
* Tajima's D (Use flag `-d` or `--tajima`)

To calculate any or all of these statistics, use the sample command of pgcalc, with the
flags of the statistics you want to calculate.

```sh
./pgcalc myfile.fas sample -s -n -d
```

or

```sh
./pgcalc myfile.fas sample -snd
```

In the above examples, an output file called myfile_sample.csv will be generated.

### Computing pairwise statistics for your sequences

You can compute pairwise statistics for your sequences file with the `pairwise` command.

Currently supported commands are:

* The number of polymorphisms between two sequences. (Use flag `-m` or `--mut`)
* dNdS statistics, computed according to the method of Nei and Gojoborei. (Use flag `-d` or `--dnds`)
* dNdS statistics, according to Nei and Gojoborei, with a correction for low
  levels of polymorphism. (Use flag `-D` or `--dndscor`)

Example:

```sh
./pgcalc myfile.fas pairwise -dm
```

or with dNdS correction:

```sh
./pgcalc myfile.fas pairwise -dDm
```

In the above examples, an output file called myfile_pairwise.csv will be generated.

### Computing dNdS statistics on a per codon level

You can compute dN and dS statistics as with the `pairwise` command, but output
numbers for every single codon position.

You do this with the `dndscodon` command.

```sh
./pgcalc myfile.fas dndscodon
```

In the above examples, an output file called myfile_dndscodon.csv will be generated.
