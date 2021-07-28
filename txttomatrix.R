v <- readr::read_tsv('/data/Brown_lab/UKBB_LD_Alkes/1_rs670318_LD.txt.gz')
v <- as.vector(v)
is.vector(v)

df[lower.tri(df, diag = T)] <- v

print(df[1:10,1:10])
print(dim(df))

