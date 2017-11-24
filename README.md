# SigSeekr

SigSeekr is now new and improved - it uses kmers to find sequences in a set of inclusion genomes that are not present in an exclusion set.

### Installation

Just clone this repository - you'll need the SigSeekr.py script.

### External Dependencies

To run SigSeekr, you will need to have the following external programs installed and present on your $PATH:
- BBTools >= 37.23 (https://jgi.doe.gov/data-and-tools/bbtools/)
- KMC >= 3.0 (http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=kmc&subpage=download)
- Python >= 3.5
- bedtools >= 2.26.0 (https://github.com/arq5x/bedtools2/releases/download/v2.26.0/bedtools-2.26.0.tar.gz)
 
### Python Package Dependencies

Included in requirements.txt - to install, use pip: `pip install -r requirements.txt`

### Usage

SigSeekr needs a folder containing genomes you want signature sequences (`--inclusion`) for and a folder containing genomes you do not want the signature sequences
to match to (`--exclusion`). The genomes in these folders must be in FASTA format, end with .fasta, and be uncompressed.

The output folder you specify will end up containing your signature sequences - they will be in a file called `unique_sequence.fasta`.

```
usage: SigSeekr.py [-h] -i INCLUSION -e EXCLUSION -o OUTPUT_FOLDER

optional arguments:
  -h, --help            show this help message and exit
  -i INCLUSION, --inclusion INCLUSION
                        Path to folder containing genome(s) you want signature
                        sequences for. Genomes must be in FASTA format
                        (.fasta), and should not be compressed.
  -e EXCLUSION, --exclusion EXCLUSION
                        Path to folder containing exclusion genome(s) - those
                        you do not want signature sequences for. Genomes must
                        be in FASTA format (.fasta) and should not be
                        compressed.
  -o OUTPUT_FOLDER, --output_folder OUTPUT_FOLDER
                        Path to folder where you want to store output files.
                        Folder will be created if it does not exist.
```

