library(magrittr)
set.seed()
imp <- mice(boys, m = 1)
data <- complete(imp)

indic <- sample(1:748, 375)
d1 <- data[indic,]
d2 <- data[-indic,]

data %$% lm(age ~ wgt + tv) %>% summary
A <- d1 %$% lm(age ~ wgt + tv) %>% summary
B <- d2 %$% lm(age ~ wgt + tv) %>% summary


list(A = d1 %$% lm(age ~ wgt + tv),
     B = d2 %$% lm(age ~ wgt + tv)) %>% pool

(coef(A)*375 + coef(B)*373)/748
(coef(A) + coef(B))/2

