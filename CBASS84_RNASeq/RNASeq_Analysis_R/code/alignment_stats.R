library(tidyverse)
library(ggplot2)
library(reshape2)
library(fs)
library(rjson)
stats.file <- read_delim("data/kallisto_alignment.tsv",delim = "\t") %>%
  separate(Category,sep = "_", into = c("a","b","c","Sample","d")) %>%
  select(Sample,Pseudoaligned,`Not aligned`) %>%
  arrange(desc(Pseudoaligned)) %>%
  mutate(`% aligned` = (Pseudoaligned/(Pseudoaligned + `Not aligned`))*100)

stats.file2 <- melt(stats.file,id.vars = "Sample", measure.vars = c("% aligned"), variable.name = "Type", value.name = "% of reads aligned")

pdf(file="data/kallisto_alignment_summary.pdf")
p <- ggplot(stats.file2, aes(Sample,`% of reads aligned` , fill=Type))
p <- p + geom_bar(stat="identity")
p <- p + theme(axis.text.x = element_text(angle = 90, size=6))
p <- p + coord_flip()
p

dev.off()


#### Get list of all snp callings from bcftools
snp.calls <- dir_ls(path = "/Volumes/omics4tb2/Collaborations/Vulcan/data/181203_K00235_0152_BHWNJ7BBXX_CBASS84_RNASeq/", recurse = T,regexp = "*samtools_final.vcf.gz")
snp.table <- as_tibble(snp.calls) %>%
  separate(value, sep = "/", into = c("a","b","c","d","e","f","g","folder","i","file")) %>%
  select(folder, file)

snp.folder <- "~/Downloads/coral_snp_calls"
dir_create(path = snp.folder)

for(snp.call in snp.calls){
  file_copy(snp.call, new_path = snp.folder)

}


#### Get list of all kallisto abundance measures
counts.files <- dir_ls(path = "/Volumes/omics4tb2/Collaborations/Vulcan/data/CBASS84_RNASeq_Final/", recurse = T,regexp = "*abundance.tsv")
counts.table <- as_tibble(counts.files) %>%
  separate(value, sep = "/", remove = F, into = c("a","b","c","d","e","f","g","folders","i","file")) %>%
  select(value,folders, file)

alignment.stats <- tibble()
for(folder in counts.table$folders){
  sample = folder
  kallisto.file <- counts.table %>%
    filter(folders == folder) %>%
    pull(value)

  print(sample)

  log.file <- sub("abundance.tsv","run_info.json",paste(kallisto.file, sep = ""))

  read.counts <- read_delim(kallisto.file, delim = "\t")

  total.reads.processed <- fromJSON(file=paste(log.file, sep = ""))["n_processed"][[1]]


  spis.reads <- filter(read.counts, grepl("Spis", target_id))
  spis.aligned.total <- round(sum(spis.reads$est_counts), digits = 0)
  spis.alignet.pct <- round((spis.aligned.total/total.reads.processed)*100, digits = 2)

  smic.reads <- filter(read.counts, grepl("Smic", target_id))
  smic.aligned.total <- round(sum(smic.reads$est_counts), digits = 0)
  smic.alignet.pct <- round((smic.aligned.total/total.reads.processed)*100, digits = 2)

  alignment.stats <- bind_rows(alignment.stats,
                               bind_cols(Sample = sample,
                                         Total.reads = total.reads.processed,
                                         Spis.Aligned.Total= spis.aligned.total,
                                         Spis.Aligned.Percent = spis.alignet.pct,
                                         Smic.Aligned.Total = smic.aligned.total,
                                         Smic.Aligned.Percent = smic.alignet.pct
                                         ))

}

write_delim(alignment.stats, file="/Volumes/omics4tb2/Collaborations/Vulcan/data/CBASS84_RNASeq_Final/kallisto_alignment_stats.txt", delim = "\t")





snp.folder <- "~/Downloads/coral_snp_calls"
dir_create(path = snp.folder)

for(snp.call in snp.calls){
  file_copy(snp.call, new_path = snp.folder)

}





