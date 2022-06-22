using WienerARD

m=Model(2,3,4,5,.2,PeriodicMaintenancePolicy(.2,12))

@run df=simulate(m,100,1)
