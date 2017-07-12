# USS-Ppip
Unique Strain-Specific Python Probe Identifying Program

## Requirements
- mash (https://github.com/marbl/mash) installed and present in your $PATH.
- Genomes folder in your working directory containing *.fa files
- Target Genome(s)
- BLASTn and Biopython


## Purpose
Allows the user to find probes that will differentiate between closely related genomes.

## Use
Core.py is the main handler for typing and unique sequence searching methods
* Output folder where a csv typing file and Unique sequences folder will be placed `(-o, --output)`
* Input folder containing *Genomes* folder `(-i, --input)`
* Target genome that you want to make probes for `(-t, --target)`

#### Optional Parameters
* Evalue start (default 1e-40) `(-e, --evalue)`
* Evalue stop (default 1e-90) `(-s, --estop)`


#### Example
`Core.py -o USS_Ppip/ -i USS_Ppip/ -t Ecoli_EDL933.fa`
