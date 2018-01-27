library(simts)
library(gmwm)
library(wv)
library(classimu)



n1 = 10000
n2 = 10000
n3 = 1000000

model =  AR1(.85,sigma2 = 1 ) + WN(.5) + RW (1e-4)
model1 = AR1() + WN() + RW ()


Xt =  gen_gts(n1, model)
Yt =  gen_gts(n2, model)
Zt =  gen_gts(n3, model)

mimu = make_wvar_mimu_obj(Xt, Yt, Zt, freq = 100, unit = "s",
         sensor.name = "MTiG - Gyro. X", exp.name = c("today", "yesterday", "a few days ago"))


test.optim = mgmwm(model1, mimu)

obj_list = test.optim

obj = mimu







