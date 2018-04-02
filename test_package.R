library(simts)
library(gmwm)
library(wv)
library(classimu)
library(progress)
library(iterpc)

n1 = 1000000
n2 = 1000000
n3 = 1000000

model1 =  AR1(.995,sigma2 = 1e-6) + WN(.005) + RW (1e-7)
model = AR1() + WN() + RW ()

tamer = c(inv_transform_phi(.095), log(1e-6) ,log(.005), log(1e-7))

Xt =  gen_gts(n1, model1)
Yt =  gen_gts(n2, model1)
Zt =  gen_gts(n3, model1)
Wt =  gen_gts(n3, model1)

mimu = make_wvar_mimu_obj(Wt ,Xt, Yt, Zt, freq = 100, unit = "s",
         sensor.name = "MTiG - Gyro. X", exp.name = c("today", "yesterday", "a few days ago"))


test.optim = mgmwm(model, mimu, stationarity_test = FALSE, B = 30, alpha_near_test = 0.05)

plot(test.optim)

# Resutlat du test

test_model_selection = model_selection(mimu,model,s_test = 2, test_pval = TRUE)


plot(test_model_selection, process.decomp = TRUE)



