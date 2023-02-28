var_normalize = function(x, poly_n = 3, print_coefs = F) {
  require(tidyverse)
  x = x + abs(min(x)) + 1
  ordernorm = x %>% bestNormalize::orderNorm(warn = F) %>% pluck("x.t")
  fit = lm(ordernorm ~ poly(log10(x), poly_n, raw = T), NULL)
  if (print_coefs == T) {print(coef(fit))}
  normalized = predict(fit, as.data.frame(x))
  return(normalized)
}

# set.seed(123)
# df = bind_rows(rgamma( 1500, 3)       %>% enframe() %>% mutate(group = "A", .before = name), 
#                (rt(    1500, 3)+4)    %>% enframe() %>% mutate(group = "B", .before = name), 
#                (rgamma(1500, 3)^2-40) %>% enframe() %>% mutate(group = "C", .before = name))
# df %>% ggplot(aes(value)) + geom_histogram(bins = 50) + facet_grid(.~group, scales = "free")
# 
# df = df %>% group_by(group) %>% mutate(value_norm = value %>% var_normalize())
# 
# df %>% ggplot(aes(value_norm)) +        geom_histogram() + facet_grid(.~group, scales = "free")
# df %>% ggplot(aes(value, value_norm)) + geom_point() + facet_grid(.~group, scales = "free")
# df %>% ggplot(aes(value_norm)) +        geom_histogram()
