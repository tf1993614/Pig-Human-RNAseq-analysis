workflow RNApipeline{
  File rawdata
  String thread
  String name
  String indexdir
  String AfterQCSoftwareDir

  call fastqc{
    input:
      fastq=rawdata,
      thread=thread
  }

  call trimming{
    input:
      fastq=rawdata,
      thread=thread,
      name=name,
      AfterQC=AfterQCSoftwareDir
  }

  call alignment{
    input:
      trimFile=trimming.trimfile,
      thread=thread,
      name=name,
      index=indexdir
  }
}

task fastqc{
  File fastq
  String thread
  command {
    fastqc ${fastq} --threads ${thread}
  }
}

task trimming{
  File fastq
  String thread
  String name
  String AfterQC
  command {
    python ${AfterQC} -1 ${fastq} --read1_flag ${name}
  }
  output {
    File trimfile = "${name}_R1.good.fq.gz"
  }
}

task alignment{
  File trimFile
  String thread
  String name
  String index
  command {
    kallisto quant -i ${index} --single -l 200 -s 20 ${trimFile}
  }
}
