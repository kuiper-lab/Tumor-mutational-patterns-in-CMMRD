##STATISTICS

wilcox.test(<Variable1> ~ <Variable2>, data = <DATASET>)

kruskal.test(<Variable1> ~ <Variable2>, data = <DATASET>)

fisher.test(<CONTINGENCYTABLE>)

library(FSA)
dunnTest(<Variable1> ~ <Variable2>, data = <DATASET>, method = "bonferroni")

library(contingencytables)
FisherFreemanHalton_asymptotic_test_rxc(<CONTINGENCYTABLE>)
