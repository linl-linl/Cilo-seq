# 导入所有reads的count值，header = T保证后续计算不会将列名算入而造成错误
# 参数：read.table(file="")
n <- read.table(file="ERCCexon_length_gtf_fadj.txt", header = T, row.names = 1)
count1 <- read.table(file="1012ZQQ_L4_count.txt", row.names = 1)
count2 <- read.table(file="1013ZQQ_L4_count.txt", row.names = 1)
count3 <- read.table(file="1016ZQQ_L3_count.txt", row.names = 1)

countToFpkm <- function(counts, effLen)
{
    N <- sum(counts)
    exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}

count5 <- mapply(countToFpkm, count1, n)
count5 <- data.frame(count5)
rn5 <- rownames(count1)
rownames(count5) <- rn5
count6 <- mapply(countToFpkm, count2, n)
count6 <- data.frame(count6)
rn6 <- rownames(count2)
rownames(count6) <- rn6
count7 <- mapply(countToFpkm, count3, n)
count7 <- data.frame(count7)
rn7 <- rownames(count3)
rownames(count7) <- rn7

# 生成表格
# 参数：filename
write.table(count5,"1012ZQQ_L4_FPKM.txt", row.names = TRUE, col.names = FALSE, sep = '\t')
write.table(count6,"1013ZQQ_L4_FPKM.txt", row.names = TRUE, col.names = FALSE, sep = '\t')
write.table(count7,"1016ZQQ_L3_FPKM.txt", row.names = TRUE, col.names = FALSE, sep = '\t')

