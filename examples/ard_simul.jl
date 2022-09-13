using WienerARD

m=Model(2,3,4,5,.7,PeriodicMaintenancePolicy(.7,3))
typeof(m)
ΔX(m, 1.0)
df=simulate(m,100,1)
