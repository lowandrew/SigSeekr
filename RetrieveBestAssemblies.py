__author__ = 'mikeknowles'
""" Retrieve Adam's Best Assemblies
"""
import glob, os, shutil


def retriever(genomes, output):
    for folders in glob.glob(genomes):
        if os.path.exists(folders + "Best_Assemblies"):
            for fasta in glob.glob(folders + "Best_Assemblies"):
                shutil.copy(fasta, output)