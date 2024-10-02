# 获取当前目录下所有的 .txt 文件
files <- list.files(pattern = "\\.txt$")

# 循环处理每个 .txt 文件
for (file in files) {
  # 读取文件
  df <- read.delim(file, header = TRUE, row.names = 1)
  
  # 截取需要的部分
  df_cut <- df[, 1:2]
  
  # 计算RPK
  df_cut$RPK <- df_cut$mapped_read * 1000 / df_cut$length
  
  # 移除最后一行
  n <- nrow(df_cut)
  result <- df_cut[-n, ]
  
  # 计算TotalRPK和TPM
  TotalRPK <- sum(result$RPK)
  result$TPM <- (result$RPK * 10e6) / TotalRPK
  
  # 生成输出文件名
  output_file <- sub("\\.txt$", "_TPM.csv", file)
  
  # 写入CSV文件
  write.csv(result, file = output_file)
}
