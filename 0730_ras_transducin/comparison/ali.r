## name: ali.r
## date: 08/12/2015
## date: 08/19/2015

load("/u2/hyangl/project/ras/results/2015/0730_ras_transducin/transducin/gt_networks.RData")

align <- read.fasta("/u2/hyangl/project/ras/results/2015/0730_ras_transducin/alignment/ras_transducin.fa")
align$ali[c("1TND_A","1QRA_A"),]
# we can see 1TND and 1TAD are well aligned.

# only select the useful regions
ali <- align$ali[c("1TND_A","1QRA_A"), c(30:117,119:366)]

# row 3: resno in Gt(1TND)
# row 4: resno in Ras(1QRA, the same as 5P21/4Q21)
# row 5: membership of the 146 aligned residues
# row 6: membership of Gt
# row 7: membership of Ras

ali <- rbind(ali, matrix("0", nrow=5, ncol=dim(ali)[2]))
rownames(ali)[3:7] <- c("resno_gt", "resno_ras", "membership_146",
                        "membership_gt", "membership_ras")

# in 1TND, resno begins at 27
ali[3,(ali[1,]!="-")] <- 27:349
ali[4,(ali[2,]!="-")] <- 1:166

membership_ras_146 <- rep(0,166)
membership_ras_146[c(5:9,18:26,29:30,38:46,50:56)] <- 1
membership_ras_146[10:17] <- 2
membership_ras_146[c(31:37)] <- 3
membership_ras_146[c(57:75)] <- 4
membership_ras_146[c(76:85,110:121,140:153)] <- 5
membership_ras_146[c(86:103,108:109)] <- 6
membership_ras_146[c(126:139)] <- 7
membership_ras_146[c(154:163)] <- 8
inds <- which(membership_ras_146 != 0)
membership_146 <- membership_ras_146[membership_ras_146!=0]
for (i in 1:length(membership_146)) {
  ali["membership_146", which(ali["resno_ras",]==inds[i])] <- membership_146[i]
}

membership_gt <- rep(0,339)
membership_gt[c(31:35,44:53,180:195)] <- 1
membership_gt[36:43] <- 2
membership_gt[c(54:72,148:172)] <- 10
membership_gt[c(73:107,110:147)] <- 11
membership_gt[173:179] <- 3
membership_gt[196:214] <- 4
membership_gt[c(215:226,259:274,316:329)] <- 5
membership_gt[227:258] <- 6
membership_gt[275:299] <- 9
membership_gt[c(300:307,309:312,314:315)] <- 7
membership_gt[330:339] <- 8
ali["membership_gt",which(ali["resno_gt",]!="0")[5:313]] <- membership_gt[31:339]

membership_ras <- rep(0,166)
membership_ras[c(1:9,18:24,38:56)] <- 1
membership_ras[10:17] <- 2
membership_ras[25:37] <- 3
membership_ras[57:75] <- 4
membership_ras[c(76:85,110:121,140:153)] <- 5
membership_ras[86:109] <- 6
membership_ras[122:130] <- 9
membership_ras[131:139] <- 7
membership_ras[154:166] <- 8
ali["membership_ras",ali["resno_ras",]!="0"] <- membership_ras

save(ali,
     file="ali.RData")

