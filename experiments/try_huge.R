# ----
library(huge)

# ----
net <- huge.generator(
  n = 20,
  d = 5000,
  graph = "random",
  v = NULL,
  u = NULL,
  g = NULL,
  prob = 1/500,
  vis = FALSE,
  verbose = TRUE
)

# ----

sum(net$sigma>0)
