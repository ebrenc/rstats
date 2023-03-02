var_normalize = function(x, poly_n = 3, get_coefs = F) {
  require(tidyverse)
  x = x + abs(min(x, na.rm=T)) + 1
  ordernorm = x %>% bestNormalize::orderNorm(warn = F) %>% pluck("x.t")
  fit = lm(ordernorm ~ poly(log10(x), poly_n, raw = T), NULL)
  poly = poly(log10(x), poly_n, raw = T) %>% as_tibble() %>% pluck("1")
  normalized = predict(fit, as.data.frame(x))
  if (get_coefs == T) {result = list(norm = normalized, coef = enframe(coef(fit), name = NULL) %>% rename(coef = value), poly = poly)} else {result = normalized}
  return(result)
}
