import os
from os import walk
from os import path
import argparse


# set parameters
parser = argparse.ArgumentParser()
parser.add_argument("-inputPath", help="Path to all the input fastq files.")
parser.add_argument("-templatejson", help="template json file, can be generated from wdltools inputs command.")
parser.add_argument("-outputPath", default = "./", help="output path for json files.")
parser.add_argument("-mappingIndex", default = "/hpcfs/users/a1763076/genome_index/Kallisto/pig/pig_kallisto_ensembl_index", help="index path for alignment.")
parser.add_argument("-AfterQCSoftwareDir", default = "/hpcfs/users/a1763076/local/virtualenvs/python_package/lib/python3.8/AfterQC-master/after.py", help="annotation file for featurecount.")

args = parser.parse_args()

# get input fastq files with path
mypath = args.inputPath

if not mypath.endswith("/"):
    mypath = mypath + "/"

files = []
for (dirpath, dirnames, filenames) in walk(mypath):
    for file in filenames:
        if file.endswith("fastq.gz"):
            files.extend(filenames)


def getinput(file):
    name = file.split("/")[-1].rstrip(".fastq.gz")
    outputfile = args.outputPath+'/'+name+".json"
    if path.exists(outputfile):
        os.remove(outputfile)
    with open(args.templatejson, "r") as template, open(outputfile, "a") as outjson:
        for line in template:
            if line.startswith("{"):
                print("{", file=outjson)
            else:
                data = line.split(": ")
                if "rawdata" in data[0]:
                    data[1]="\""+file+"\""+","
                if "thread" in data[0]:
                    data[1]="\""+"16"+"\""+","
                if "indexdir" in data[0]:
                    data[1]="\""+args.mappingIndex+"\""+","
                if "name" in data[0]:
                    data[1]="\""+name+"\""+","
                if "AfterQCSoftwareDir" in data[0]:
                    data[1]="\""+args.AfterQCSoftwareDir+"\""
                print(": ".join(i for i in data), file=outjson)

for file in files:
    fullpathfile = mypath+file
    getinput(fullpathfile)
