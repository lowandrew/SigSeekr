# USS-Ppip
Unique Strain-Specific Python Probe Identifying Program

## Requirements
- rMLST data from [PubMLST](http://pubmlst.org/rmlst/) (requires an account) in *.fasta format
- Genomes folder in your working directory containing *.fa files
- Target Genome(s)
- BLASTn and Biopython


## Purpose
Allows the user to find probes that will differentiate between closely related genomes. rMLST is required for the software to run typing methods. A mass downloader is recommended for downloading rMLST genes

## Use
Core.py is the main handler for typing and unique sequence searching methods
* Output folder where a csv typing file and Unique sequences folder will be placed `(-o, --output)`
* Input folder containing *Genomes* folder `(-i, --input)`
* Markers folder containing all 53 rMLST markers in fasta format `(-m, --markers)`

 * Other markers may be used as well, anything to differentiate genomes or proteins

#### Optional Parameters
* Evalue start (default 1e-40) `(-e, --evalue)`
* Evalue stop (default 1e-90) `(-s, --estop)`


#### Example
`Core.py -o USS_Ppip/ -i USS_Ppip/ -m /rMLST/ -t Ecoli_EDL933.fa`
