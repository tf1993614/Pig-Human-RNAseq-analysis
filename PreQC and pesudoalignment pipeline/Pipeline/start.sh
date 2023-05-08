#!/bin/bash


ml Python/3.6.1-foss-2016b

# 1st para is input path
# 2nd para is template
# 3rd para is output path
input="/hpcfs/users/a1763076/test/fastq_files/"
output="/hpcfs/users/a1763076/test/json_output"
template="/hpcfs/users/a1763076/test/json_input/template.json"
index="/hpcfs/users/a1763076/genome_index/Kallisto/pig/pig_kallisto_ensembl_index"

echo "start to gerenate input json files."
python getInput.py -inputPath ${input} -templatejson ${template} -outputPath ${output} -mappingIndex ${index}
echo "finish generating json files."

for i in ${output}/*.json; do sbatch slurm_RNAseqWDL.sh ${i}; done
