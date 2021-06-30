
library(edgeR)  #load edger package

count_data = '/data_244/RNA/mapped/count/table_data/count_table.csv'

count_df = read.csv(count_data)

head(count_df)

y = DGEList(count=count_df)

y

