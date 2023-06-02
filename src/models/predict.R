
STEC <-
  dat %>%
  filter(caseDef == "STEC") %>%
  rename(y = cases) %>%
  mutate(monthInYear = as.integer(format(Date, "%m"))) %>%
  select(-caseDef)

tmp <- aeddo(data = STEC,
      theta = rep(1,7),
      formula = y ~ -1 + ageGroup,
      method = "L-BFGS-B",
      lower = c(rep(1e-6,6), -6),
      upper = rep(1e2, 7),
      model = "PoissonNormal",
      k = 36,
      excludePastOutbreaks = TRUE, CI = FALSE)

# library(readr)
# write_rds(x = tmp, "toy_CI_example.RDS")

# tmp %>% select(ref.date, par) %>% unnest(cols = c(par)) %>%
#   ggplot(mapping = aes(x = ref.date, group = Parameter)) +
#   geom_line(mapping = aes(y = theta, colour = Parameter)) +
#   geom_line(mapping = aes(y = CI.lwr)) +
#   geom_line(mapping = aes(y = CI.upr)) +
#   facet_wrap(facets = vars(Parameter), scales = "free_y")

tmp %>% ggplot(mapping = aes(x = ref.date, y = LogS)) +
  geom_line()
  




predict.aeddo <- function(object){
  
}