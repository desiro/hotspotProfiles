# [<samp>hotspotProfiles</samp>](https://github.com/desiro/hotspotProfiles)
[![License: GPL v3](https://img.shields.io/badge/License-GPL_v3-bd0000.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python v3.9.7](https://img.shields.io/badge/Language-Python_v3-75a8d3.svg)](https://www.python.org/)
[![Conda v4.11.0](https://img.shields.io/badge/Uses-Conda-43b02a.svg)](https://docs.conda.io/en/latest/miniconda.html)

***

## Description

This tool evaluates vRNAsite like interaction tables and creates hotspot profiles for each individual sequence, showing regions with high interaction counts.

### Mandatory Prerequisites

* [![Python v3.9.7](https://img.shields.io/badge/Python_v3.9.7-75a8d3.svg)](https://www.python.org/downloads/release/python-397/)
* [![NumPy v1.22.2](https://img.shields.io/badge/NumPy_v1.22.2-013243.svg)](http://www.numpy.org/)
* [![Matplotlib v3.5.1](https://img.shields.io/badge/Matplotlib_v3.5.1-11557c.svg)](https://matplotlib.org/)
* [![Scipy v1.8.0](https://img.shields.io/badge/Scipy_v1.8.0-013243.svg)](https://scipy.org/)

### Optional Prerequisites

* [![Conda v4.11.0](https://img.shields.io/badge/Conda_v4.11.0-43b02a.svg)](https://docs.conda.io/en/latest/miniconda.html)

***

## Installation

To run <samp>hotspotProfiles</samp>, I recommend using Miniconda and following the steps below. If this is the first time using conda, you should probably restart your shell after the installation of Miniconda. The following will demonstrate the installation and set up of Miniconda on Linux, which should be similar on other platforms. For Windows 10 users, I advise using the [Ubuntu 20.04 LTS](https://www.microsoft.com/en-us/p/ubuntu-2004-lts/9n6svws3rx71?cid=msft_web_chart) subsystem. More information can be found on the [Miniconda](https://docs.conda.io/en/latest/miniconda.html) and [Bioconda](https://bioconda.github.io/user/install.html) pages.

### Conda Installation

Installing Miniconda:
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Updating Miniconda and setting channels:
```
conda update conda
conda update python
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

Installing Conda packages:
```
conda create --name hotspotProfiles python=3.9.7
conda activate hotspotProfiles
conda install -c conda-forge numpy=1.22.2
conda install -c conda-forge matplotlib=3.5.1
conda install -c conda-forge scipy=1.8.0
git clone https://github.com/desiro/hotspotProfiles.git
cd hotspotProfiles
```

***

## Examples

The basic input for ```hotspotProfiles.py``` includes the following parameters:
* ```-pfx``` - name of the output folder
* ```-fsa``` - input fasta file

### Fasta Formation

The tool always takes the first sequence in the fasta file as the target mutation sequence. If there is no specified range, the whole sequence will be taken and compared against all other sequences present in the fasta file. 

### Basic Example

```
python hotspotProfiles.py -pfx example -vst example.tsv -ovr
```

### Options

For more command line options, see the [manual](https://github.com/desiro/hotspotProfiles/blob/master/manual.md).

***

## Authors

* [Daniel Desirò](https://github.com/desiro)

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## Reference

Please cite <samp>HotspotProfiles</samp> if you find our tool useful.

```
D. Desirò, A. Borodavka and M. Marz.
"HotspotProfiles: Prediction of hotspot regions capable of mediating the selective but flexible genome packaging of influenza A viruss."
In Preparation, 2022.
```
