__author__ = 'mikeknowles'
from os.path import isfile
from glob import glob
from json import load
from shutil import copy
wgspath = "/nas/WGS_Spades/AssemblyData/"
sam = "/nas/Genomes/SalmonellaSequenced/"
count = 0

for run in glob(wgspath + "JsonReports/*.json"):
    report = load(open(run))
    ref = report["1.General"]["referenceGenome"]
    quality = report["1.General"]["averageDepthofCov"]
    contlen = report["2.Assembly"]["TotalLength"]
    numcont = report["2.Assembly"]["NumContigs"]
    if "bon" in ref:
        print ref
    if "Salmonella" in ref and float(quality) > 20.0 and 5e6 > float(contlen) > 3.8e6 and int(numcont) < 200:
        # print run.replace("metadataReport.json", "filteredAssembled.fasta").replace("JsonReports", "Assemblies")
        # print quality, numcont
        # print numcont

        count += 1
        # copy(run.replace("metadataReport.json", "filteredAssembled.fasta").replace("JsonReports", "Assemblies"), sam)
print count
