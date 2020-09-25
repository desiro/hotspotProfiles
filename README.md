# flexibleRegions

```flexibleRegions``` reduces the interaction strength between a given area of the first sequence in regard to all other sequences. The tool is written in ```Python 3.7.1``` and the calculations are performed with the ```RNAcofold``` python site-package of the ```ViennaRNA Package 2.4.13```. For command-line options, please refer to the [manual](https://github.com/desiro/flexibleRegions/blob/master/manual.md)

## Mandatory Prerequisites

* [python 3.7.1](https://www.python.org/downloads/release/python-385/)
* [viennaRNA 2.4.13](https://www.tbi.univie.ac.at/RNA/documentation.html#install)

## Optional Prerequisites

* [Miniconda3](https://docs.conda.io/en/latest/miniconda)

## Installation with Miniconda

Please download the [Miniconda3](https://docs.conda.io/en/latest/miniconda.html) application for your system. The following will demonstrate the installation and set up of miniconda on Linux, which should be similar on other platforms.

```
bash Miniconda3-latest-Linux-x86_64.sh -p ~/miniconda3
conda create --name flexibleRegions
conda activate flexibleRegions
conda install -c bioconda viennarna
```

## Examples

The basic input for ```flexibleRegions.py``` includes the following parameters:
* ```-pfx``` - name of the output folder
* ```-fsa``` - input fasta file

### Fasta Formation

The tool always takes the first sequence in the fasta file as the target mutation sequence. If there is no specified range, the whole sequence will be taken and compared against all other sequences present in the fasta file. 

### Basic Example

```
python3 flexibleRegions.py -pfx example -fsa example.fa -pss 32 -pse 96 -thr 4 -ovr 
```

## Authors

* [Daniel Desirò](https://github.com/desiro)

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## Reference

Please cite [flexibleRegions](https://doi.org/10.1101/424002) if you find our tool useful.

## Workflow overview

![workflow](https://github.com/desiro/flexibleRegions/blob/master/workflow.png "(a) creates all k-mers of the query and target sequences (b) predicts structures with RNAcofold between all query and all target k-mers (c) iteratively takes the query k-mer that has the most stable structure with a target k-mer as the regress query k-mer (d) creates all three nucleotide mutations of the middle nucleotide in the regress query and predicts structures with RNAcofold between all three mutants and all target k-mers (e.1) penalizes the current most stable query k-mer if there is no regress query k-mer mutant less stable than the original (e.2) or otherwise takes the regress query k-mer mutant with the least stable structure, replaces it with the original and recalculates all query k-mers that share the mutated nucleotide (f) terminates and returns the mutant query sequence if there is no query k-mer that is more stable than the predefined stopping energy threshold")
