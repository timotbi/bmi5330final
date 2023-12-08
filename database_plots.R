library(dplyr)

npm1_vars <- read.table("databases/dbSNP/vcf/NPM1_snvs.tsv")

npm1_snvs <- npm1_vars[grep(pattern = "VC=SNV", x = npm1_vars[, 8]), ]
npm1_mnvs <- npm1_vars[grep(pattern = "VC=MNV", x = npm1_vars[, 8]), ]
npm1_indels <- rbind(npm1_vars[grep(pattern = "VC=INDEL", x = npm1_vars[, 8]), ],
                     npm1_vars[grep(pattern = "VC=INS", x = npm1_vars[, 8]), ],
                     npm1_vars[grep(pattern = "VC=DEL", x = npm1_vars[, 8]), ])

# NPM1 SNVs Coding vs. Non-coding
ncoding_snvs_npm1 <- sum(length(grep(pattern = "NSM", npm1_snvs[,8])),
                  length(grep(pattern = "NSN", npm1_snvs[,8])),
                  length(grep(pattern = "SYN", npm1_snvs[,8])))

# NPM1 Indel Coding vs. Non-coding
ncoding_indels_npm1 <- sum(length(grep(pattern = "NSF", npm1_indels[,8])),
                           length(grep(pattern = "NSN", npm1_indels[,8])))

npm1_table <- matrix(data = c(3, 0, 
                              c(nrow(npm1_indels) - ncoding_indels_flt3), ncoding_indels_flt3,
                              c(nrow(npm1_snvs) - ncoding_snvs_npm1), ncoding_snvs_npm1), nrow = 2, ncol = 3)
rownames(npm1_table) <- c("noncoding", "coding")
colnames(npm1_table) <- c("MNVs", "Indels", "SNVs")

npm1_plot <- barplot(prop.table(npm1_table), col = c("grey", "red"), ylim = c(0, 1),
                     legend.text = c("Noncoding", "Coding"), main = "dbVar NPM1", args.legend = c(x = "topleft"))
text(x = npm1_plot, y = prop.table(npm1_table)[1, ]/2, labels = c("", paste0(npm1_table[1,2:3], " NC")))
text(x = npm1_plot, y = colSums(prop.table(npm1_table)) + 0.025, labels = c("3 NC", paste0(npm1_table[2, 2:3], " C")))

# NPM1 SNV coding breakdown
pie(x = c(length(grep(pattern = "NSM", npm1_snvs[,8])),
          length(grep(pattern = "NSN", npm1_snvs[,8])),
          length(grep(pattern = "SYN", npm1_snvs[,8]))),
    labels = c(paste0("Missense\n", length(grep(pattern = "NSM", npm1_snvs[,8]))), 
              paste0("Nonsense\n", length(grep(pattern = "NSN", npm1_snvs[,8]))),
              paste0("Synonymous\n", length(grep(pattern = "SYN", npm1_snvs[,8])))),
    #cex = 1.8,)
)

# NPM1 Indel coding breakdown
pie(x = c(length(grep(pattern = "NSF", npm1_indels[,8])),
          length(grep(pattern = "NSN", npm1_indels[,8]))),
    label = c(paste0("Non-Syn\nFrameshift\n", length(grep(pattern = "NSF", npm1_indels[,8]))),
              paste0("Nonsense\n", length(grep(pattern = "NSN", npm1_indels[,8])))),
    cex = 1.6)


## FLT3 Plots
flt3_vars <- read.table("databases/dbSNP/vcf/FLT3_snvs.tsv")

flt3_snvs <- flt3_vars[grep(pattern = "VC=SNV", x = flt3_vars[, 8]), ]
flt3_mnvs <- flt3_vars[grep(pattern = "VC=MNV", x = flt3_vars[, 8]), ]
flt3_indels <- rbind(flt3_vars[grep(pattern = "VC=INDEL", x = flt3_vars[, 8]), ],
                     flt3_vars[grep(pattern = "VC=INS", x = flt3_vars[, 8]), ],
                     flt3_vars[grep(pattern = "VC=DEL", x = flt3_vars[, 8]), ])

# FLT3 SNVs Coding vs. Non-coding
ncoding_snvs_flt3 <- sum(length(grep(pattern = "NSM", flt3_snvs[,8])),
                         length(grep(pattern = "NSN", flt3_snvs[,8])),
                         length(grep(pattern = "SYN", flt3_snvs[,8])))

ncoding_indels_flt3 <- sum(length(grep(pattern = "NSF", flt3_indels[,8])),
                           length(grep(pattern = "NSN", flt3_indels[,8])))

ncoding_mnvs_flt3 <- length(grep(pattern = "NSM", flt3_mnvs[,8]))

flt3_table <- matrix(data = c((nrow(flt3_mnvs) - ncoding_mnvs_flt3), ncoding_mnvs_flt3,
                              (nrow(flt3_indels) - ncoding_indels_flt3), ncoding_indels_flt3,
                              (nrow(flt3_snvs) - ncoding_snvs_flt3), ncoding_snvs_flt3), ncol = 3, nrow = 2)
rownames(flt3_table) <- c("noncoding", "coding")
colnames(flt3_table) <- c("MNVs", "Indels", "SNVs")

flt3_plot <- barplot(prop.table(flt3_table), col = c("grey", "red"), ylim = c(0, 1),
                     legend.text = c("Noncoding", "Coding"), main = "dbVar FLT3", args.legend = c(x = "topleft"))
text(x = flt3_plot, y = prop.table(flt3_table)[1, ]/2, labels = c("", paste0(flt3_table[1,2:3], " NC")))
text(x = flt3_plot, y = colSums(prop.table(flt3_table)) + 0.025, labels = c("10 NC, 1C", paste0(flt3_table[2, 2:3], " C")))

# FLT3 SNV coding breakdown
pie(x = c(length(grep(pattern = "NSM", flt3_snvs[,8])),
          length(grep(pattern = "NSN", flt3_snvs[,8])),
          length(grep(pattern = "SYN", flt3_snvs[,8]))),
    labels = c(paste0("Missense\n", length(grep(pattern = "NSM", flt3_snvs[,8]))), 
               paste0("Nonsense\n", length(grep(pattern = "NSN", flt3_snvs[,8]))),
              paste0("Synonymous\n", length(grep(pattern = "SYN", flt3_snvs[,8])))), 
    cex = 1)

# FLT3 indel coding breakdown
pie(x = c(length(grep(pattern = "NSF", flt3_indels[,8])),
          length(grep(pattern = "NSN", flt3_indels[,8]))),
    label = c(paste0("Non-Syn\nFrameshift\n", length(grep(pattern = "NSF", flt3_indels[,8]))),
              paste0("Nonsense\n", length(grep(pattern = "NSN", flt3_indels[,8])))),
    cex = 1.6)


#clinVar stuff
clinVar <- read.table("databases/clinVar/vcf/clinvar.vcf.gz")
npm1_vars_clin <- npm1_vars[grep(pattern = "CLN", x = npm1_vars[, 8]), ]
clinvar_npm1 <- clinVar[grep(pattern = "NPM1", x = clinVar$V8), ]

npm1_vars_clin_mut_class <- numeric(length = 23)
for(i in 0:22){
  npm1_vars_clin_mut_class[i + 1] <- length(grep(paste0("CLNSIG=[[:punct:]]*", i), npm1_vars_clin$V8))
}
npm1_vars_clin[grep(paste0("CLNSIG=[[:punct:]]*", 5), npm1_vars_clin$V8), 8]

flt3_vars_clin <- flt3_vars[grep(pattern = "CLN", x = flt3_vars[, 8]), ]
clinvar_flt3 <- clinVar[grep(pattern = "FLT3:", x = clinVar$V8), ]

flt3_vars_clin_mut_class <- numeric(length = 23)
for(i in 0:22){
  flt3_vars_clin_mut_class[i + 1] <- length(grep(paste0("CLNSIG=[[:punct:]]*", i), flt3_vars_clin$V8))
}
flt3_path <- rbind(flt3_vars_clin[grep(paste0("CLNSIG=[[:punct:]]*", 4), flt3_vars_clin$V8), ],
                   flt3_vars_clin[grep(paste0("CLNSIG=[[:punct:]]*", 16), flt3_vars_clin$V8), ],
                   flt3_vars_clin[grep(paste0("CLNSIG=[[:punct:]]*", 5), flt3_vars_clin$V8), ])



# COSMIC stuff
cosmic_npm1 <- do.call(rbind, 
                       strsplit(x = readLines("databases/COSMIC/NPM1_consensus_mutations.tsv"), split = "\t"))
cosv_npm1 <- unique(cosmic_npm1[, 7])
cosv_npm1 <- cosv_npm1[-which(cosv_npm1 == "")]

n_per_cosv <- numeric()
for(id in cosv_npm1){
  n_per_cosv <- c(n_per_cosv, length(grep(id, cosmic_npm1[, 7])))
}
n_per_cosv[order(n_per_cosv, decreasing = T)[1:20]]

cosv_phenotype <- sapply(X = cosv_npm1[order(n_per_cosv, decreasing = T)][1:20], 
                         function(x){cosmic_npm1[cosmic_npm1[, 7] == x, ][1, 11]})

cosv_npm1_pie <- pie(x = n_per_cosv[order(n_per_cosv, decreasing = T)[1:20]], 
                     labels = cosv_phenotype[1:4], cex = 1.5)

plot.new()
text(x = 0.5, y = seq(0.94, 0, -0.06), labels = paste0(seq(5, 20, 1), ". ", cosv_phenotype[5:20]), pos = 4)


cosmic_flt3 <- do.call(rbind, 
                       strsplit(x = readLines("databases/COSMIC/FLT3_consensus_mutations.tsv"), split = "\t"))
cosv_flt3 <- unique(cosmic_flt3[, 7])
cosv_flt3 <- cosv_flt3[-which(cosv_flt3 == "")]

n_per_cosv <- numeric()
for(id in cosv_flt3){
  n_per_cosv <- c(n_per_cosv, length(grep(id, cosmic_flt3[, 7])))
}
cosv_phenotype <- sapply(X = cosv_flt3[order(n_per_cosv, decreasing = T)][1:20], 
                         function(x){cosmic_flt3[cosmic_flt3[, 7] == x, ][1, 11]})
pie(n_per_cosv[order(n_per_cosv, decreasing = T)[1:20]], labels = cosv_phenotype[1:6], cex = 1.5)

plot.new()
text(x = 0.5, y = seq(0.80, 0, -0.06), labels = paste0(seq(7, 20, 1), ". ", cosv_phenotype[7:20]), pos = 4)



for(i in cosv_flt3[order(n_per_cosv, decreasing = T)][1:20]){
  print(cosmic_flt3[cosmic_flt3[, 7] == i, ][1, 11])
}

cosmic_flt3