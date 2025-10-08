
# devtools::source_url("https://raw.githubusercontent.com/ebrenc/rstats/main/var_normalize.R")

# df = df %>% mutate(vot = VOT_ms %>% var_normalize(get_coefs = F))
# norm = df2$ori %>% var_normalize()
# denorm = var_denormalize(x = norm$norm, coefs = norm$coef, poly = norm$poly)

var_normalize = function(x, poly_n = 3, get_coefs = FALSE) {
  library(tidyverse)
  if (!require(bestNormalize, quietly = TRUE)) { install.packages("bestNormalize", dependencies = TRUE) }
  library(bestNormalize)
  ordernorm = x %>% bestNormalize::orderNorm(warn = FALSE) %>% pluck("x.t")
  fit = lm(ordernorm ~ poly(x, poly_n, raw = TRUE), NULL)
  poly = poly(x, poly_n, raw = TRUE) %>% as_tibble() #%>% pluck("1")
  normalized = predict(fit, as.data.frame(x))
  if (get_coefs == TRUE) {result = list(
    norm = normalized, 
    coef = enframe(coef(fit), name = NULL) %>% rename(coef = value), 
    poly = poly)
  } else {result = normalized}
  return(result)
}

var_denormalize = function(x, coefs, poly) {
  library(tidyverse)
  if (!require(bestNormalize, quietly = TRUE)) { install.packages("bestNormalize", dependencies = TRUE) }
  library(bestNormalize)
  coefs = as.vector(coefs$coef)
  poly_function = function(x, coefs) { sum(coefs * x^(0:(length(coefs) - 1))) }
  denormalized = numeric(length(x))
  for (i in 1:length(x)) {
    target = x[i]
    func = function(x) poly_function(x, coefs) - target
    denormalized[i] = uniroot(func, lower = -1e6, upper = 1e6)$root
  }
  return(denormalized)
}
