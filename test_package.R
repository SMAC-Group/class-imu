library(simts)
library(gmwm)
library(wv)
library(classimu)
library(progress)



n1 = 100000
n2 = 10000
n3 = 10000

model1 =  AR1(.85,sigma2 = 1e-4 ) + WN(.005) + RW (1e-7)
model = 3*AR1() + WN() + RW ()


Xt =  gen_gts(n1, model1)
Yt =  gen_gts(n2, model1)
Zt =  gen_gts(n3, model1)
Wt =  gen_gts(n3, model1)

mimu = make_wvar_mimu_obj(Wt ,Xt, Yt, Zt, freq = 100, unit = "s",
         sensor.name = "MTiG - Gyro. X", exp.name = c("today", "yesterday", "a few days ago"))


test.optim = mgmwm(model, mimu, stationarity_test = FALSE, B = 30, fast = F, alpha_near_test = 0.05)

# Resutlat du test

test_model_selection = model_selection(mimu,model,s_test = 2, test_pval = TRUE)

class(test_model_selection)

plot(test_model_selection, process.decomp = TRUE)



