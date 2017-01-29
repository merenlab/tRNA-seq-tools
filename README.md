# tRNA-seq Tools

This is a library and a set of scripts designed to process tRNA-seq data. Currently it is quite an early development stage, but feel free to e-mail us if you have questions.

## Installation

The codebase requires Python 3.5 or later. Follow these steps to get it working:

```bash
# setup some directories
mkdir ~/github/
mkdir ~/virtual-envs/

# get the code
cd ~/github/
git clone https://github.com/merenlab/tRNA-seq-tools.git

# start a new virtual environment
virtualenv-3.5 ~/virtual-envs/tRNA-seq-tools

# activate the env temporarily to install requirements:
source ~/virtual-envs/tRNA-seq-tools/bin/activate
pip install -r tRNA-seq-tools/requirements.txt
deactivate

# update the activation script
echo 'export PYTHONPATH="$HOME/github/tRNA-seq-tools"' >> ~/virtual-envs/tRNA-seq-tools/bin/activate
echo 'export PATH="$PATH:$HOME/github/tRNA-seq-tools/bin"' >> ~/virtual-envs/tRNA-seq-tools/bin/activate

# add an alias to your profile for easy activation:
echo 'alias trna-seq-tools-activate="source ~/virtual-envs/tRNA-seq-tools/bin/activate"' >> ~/.bash_profile

```

Now you can activate your virtual environment and deactivate it the following way:

``` bash
$ trna-seq-tools-activate
(tRNA-seq-tools) $ trna-profile -h
usage: trna-profile [-h] [-n SAMPLE_NAME] [-o OUTPUT_PATH] readfile

Sort tRNAs

positional arguments:
  readfile              name of read file

optional arguments:
  -h, --help            show this help message and exit
  -n SAMPLE_NAME, --sample_name SAMPLE_NAME
                        sample name (to be used for naming output files
  -o OUTPUT_PATH, --output_path OUTPUT_PATH
                        output path (path where output will be redirected
(tRNA-seq-tools) $ deactivate
```

You want to run a small example?

``` bash
$ trna-seq-tools-activate
$ cd ~/github/tRNA-seq-tools/testing
$ ./01_mini_run.sh
Sample name ..................................: test_sample
Input FASTA ..................................: raw_tRNA_sequences.fa
Total raw seqs processed .....................: 1,833
Total tRNA seqs recovered ....................: 516
Total full length tRNA seqs ..................: 48
Output DB path ...............................: test-output/test_tRNA_profile.db
```