using WienerARD

m=Model(2,3,4,5,.7,PeriodicMaintenancePolicy(.7,3.0))
ΔX(m, 1)
df=simulate(m,100,1.0)
