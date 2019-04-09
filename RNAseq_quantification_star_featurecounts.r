library(stringr)

make.bsubfile <- function(n.core, job.name, cue.name, project.name, n.hrs, ptile, mem, data_dir, star_index, outpath, gtf)
{

  str.vec <- c("#BSUB -L /bin/bash",
  paste("#BSUB -n",n.core),
  paste("#BSUB -R span[ptile=",ptile,"]",sep = ""),
  paste("#BSUB -R rusage[mem=",mem,"]",sep = ""),
  paste("#BSUB -J",basename(job.name)),
  paste("#BSUB -o ",basename(job.name),".out",sep = ""),
  paste("#BSUB -e ",basename(job.name),".err",sep = ""),
  paste("#BSUB -q",cue.name),
  paste("#BSUB -P",project.name),
  paste0("#BSUB -W ",n.hrs,":00", "\n\n"),
  paste("module load star/2.7.0c fastqc/0.11.8 subread/1.6.3 stringtie/1.3.3b", 
  
  paste(paste0("fastqc ", data.dir, job.name, "_R1_001.fastq"), "-t", n.core, "-o", outpath, sep = " "), 
  paste(paste0("fastqc ", data.dir, job.name, "_R2_001.fastq"), "-t", n.core*2, "-o", outpath, sep = " "), 
  paste("STAR --genomeDir", star_index, "--runThreadN", n.core, '--readFilesIn', paste0(data.dir, job.name, "_R1_001.fastq"), paste0(data.dir, job.name, "_R2_001.fastq"), "--outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard", "--outFileNamePrefix", paste0(outpath, basename(job.name)), sep = " "),
  paste0("featureCounts -T ", n.core, " -t exon -g gene_id -s 1 -a ", gtf, " -o ", outpath, basename(job.name), "_featureCounts_exon.txt", " ", outpath, basename(job.name), "Aligned.sortedByCoord.out.bam"),
  paste("stringtie", paste0(outpath, basename(job.name), "Aligned.sortedByCoord.out.bam"), "-G", gtf, "-A", paste0(outpath, basename(job.name),"_stringtie_exon.txt"), sep = " "),
  sep = "\n\n"))

  output <- paste(str.vec,collapse = "\n")

  fileConn<-file(paste0(basename(job.name), ".lsf"))
  writeLines(output, fileConn)
  close(fileConn)
}

#set lsf parameters
ncore = 12
mem = 2500
cue.name <- "premium"
pj.name <- "acc_zhangb03a"
n.hrs <- "4"

#data and work directory
data.dir <- "/sc/orga/projects/zhangb03a/lei_guo/Dongming_Drug_RNAseq/raw/"
star.index <- "/sc/orga/projects/zhangb03a/lei_guo/GRCm38_release95/GRCm38_release95_Index_Star/"
out.path <- "/sc/orga/projects/zhangb03a/lei_guo/Dongming_Drug_RNAseq/rnaseq_quant/"
gtf <- "/sc/orga/projects/zhangb03a/lei_guo/GRCm38_release95/Mus_musculus.GRCm38.95.chr.gtf"

#fastq file names
data.files <- list.files(path = data.dir,full.names = F,pattern = "\\.fastq$", recursive = T)

file_name_split <- str_split(data.files, "_R")
file_name_unique <- unique(sapply(file_name_split, "[[", 1))

#write lsf files
sapply(file_name_unique, make.bsubfile, n.core = ncore, cue.name = cue.name, project.name = pj.name, n.hrs = n.hrs, ptile = ncore, mem = mem, data_dir = data.dir, star_index = star.index, outpath = out.path, gtf = gtf)


#make bsub.sh
bsub.cmd <- paste0("bsub", " ", "< ", paste(list.files(pattern = "\\.lsf$"), sep = "\n"))
filecon<-file(paste0("bsub_cmd", ".sh"))
writeLines(bsub.cmd, filecon)
close(filecon)

system("bash ./bsub_cmd.sh")


