using WienerARD

m=Model(2,3,4,3,.7,PeriodicMaintenancePolicy(.7,15))
ΔX(m, 1.0)
#@run df=simulate(m,100,3)
df=simulate(m,100,1)