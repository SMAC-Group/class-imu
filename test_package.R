library(simts)
library(gmwm)
library(wv)
library(classimu)



n1 = 10000
n2 = 10000
n3 = 1000000

model1 =  AR1(.85,sigma2 = 1 ) + WN(.5) + RW (1e-4)
model = 3*AR1() + WN() + RW ()


Xt =  gen_gts(n1, model1)
Yt =  gen_gts(n2, model1)
Zt =  gen_gts(n3, model1)
Wt =  gen_gts(n3, model1)

mimu = make_wvar_mimu_obj(Wt ,Xt, Yt, Zt, freq = 100, unit = "s",
         sensor.name = "MTiG - Gyro. X", exp.name = c("today", "yesterday", "a few days ago"))


test.optim = mgmwm(model, mimu, stationarity_test = FALSE, B = 30)


plot(test.optim, process.decomp = TRUE)

test_model_selection = model_selection(mimu,model,s_test = 2)





