library(simts)
library(gmwm)
library(wv)
library(classimu)

n1 = 100000
n2 = 100000
n3 = 100000

model =  AR1(.85,sigma2 = 1 ) + WN(1.5e-10) + RW (1e-7)
model1 = AR1() + WN() + RW ()


Xt =  gen_gts(n1, model)
Yt =  gen_gts(n2, model)
Zt =  gen_gts(n3, model)

obj = make_wvar_mimu_obj(Xt, Yt, Zt, freq = 100, unit = "s",
         sensor.name = "MTiG - Gyro. X", exp.name = c("today", "yesterday", "a few days ago"))
tau = 2^(1:12)

wv.theo = wv_theo(model, tau)

mimu[[2]]$variance

theta = c(.85,1,1.5e-10,1e-7)
mgmwm_obj_function (theta, model, mimu)

test.optim = matrix(NA,30,4)

for(i in 1:30){

  test.optim[i,] = mgmwm(model1, mimu)$par
}





