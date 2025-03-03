# ----
library(JGLtry)
library(data.table)

f1 = "../iddn_data/pbmc_subset/pbmc_rna_ct1_0p25_0p25.csv"
f2 = "../iddn_data/pbmc_subset/pbmc_rna_ct2_0p25_0p25.csv"

df1 = fread(f1, header=TRUE, drop=c("V1"))
df2 = fread(f2, header=TRUE, drop=c("V1"))

dat1 = scale(df1)[,1:1000]
dat2 = scale(df2)[,1:1000]
dat_lst = list(dat1, dat2)

# ----
start_time = Sys.time()
out = JGL(
  dat_lst,
  lambda1 = 0.05,
  lambda2 = 0.005,
  return.whole.theta = TRUE
)
end_time <- Sys.time()
print(end_time - start_time)
