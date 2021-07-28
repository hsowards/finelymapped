library(reticulate)

args <- commandArgs(trailingOnly=T)
np <- import("numpy")

npz <- np$load(args[1])
print("1")
lst <- npz$files
print("2")

print(lst)
print("3")

print(npz$f[["data"]][1:10,])
print(dim(npz$f[["row"]]))
print(npz$f[["col"]][1:10,])
