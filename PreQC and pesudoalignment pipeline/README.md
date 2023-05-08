## This is a pipiline to parallely run QC steps and pseudoalignment

- **WDL**: java scripts to run WDL frame work can be downloaded via this [link](https://github.com/broadinstitute/cromwell/releases/tag/85). More details can be found [here](https://github.com/openwdl/wdl)
- **Pipeline**: Python and shell scripts to chain up biominformatcis analysis steps under WDL framework
  * `start.sh`: shell script to load necessary modules and initiate the whole pipeline
  * `getInput.py`: automically generate the json files required for WDL frame work
  * `template.json`: a template json file for running `getInput.py`
  * `slurm_RNAseqWDL.sh`: shell script to initiate WDL frame work
  * `RNAseq.wdl`: WDL script to chain up anlysis steps
  



In this pipeline, those following steps were done:

- **Fastqc**: check fastq sequencing files quality
- **AfterQC**: remove adapters and low-quality reads from sequencing files
- **Kallisto**: do a pseudoalignment at transcript level
