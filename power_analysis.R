library("pwr")

cohen.ES(test = "anov", size = "large")
# f = 0.4 is a large effect size

# k = number of groups, f = effect size
pwr.anova.test(k = 4, f = 0.4, sig.level = 0.05, power = 0.8)
