function simulate(model::Model,n::Int64,p::Float64)::DataFrame
    tau_periode=model.mp.τ
    nb_maint=floor(Int64,(n - 1) / (tau_periode + 1))
    n=n-nb_maint
    t=0:p:n*p-1
   
     taus=Float64[0;t[tau_periode+1:tau_periode:length(t)-1];t[length(t)]]
     
     tausm=t[range(tau_periode+1,length=nb_maint,step=tau_periode)]
     Random.seed!(ale)
     inc=[0 0]
     cov=model.r*sqrt(model.sigma21*model.sigma22)
     for i in 1:length(diff(t))
         inc=[inc;reshape(rand(MvNormal([model.mu1*diff(t)[i],mu2*diff(t)[i]],[model.sigma21*diff(t)[i] cov*diff(t)[i] ; cov*diff(t)[i] sigma22*diff(t)[i]])),1,2)]
     end
     Xt=cumsum(inc,dims=1)
     Xt0=Xt[findall(t.<=taus[2]),:]
     Xts=Float64[] ; X1=Float64[] ; X2=Float64[]
     for i in 1:(length(taus)-1)
         Xt1=Xt[intersect(findall(t.<=taus[i+1]),findall(t.>=taus[i])),1]
         Xt2=Xt[intersect(findall(t.<=taus[i+1]),findall(t.>=taus[i])),2]
         Xtf=Xt1.-model.mp.ρ*Xt2[1]
         Xts=append!(Xts,Xtf)
         X1=append!(X1,Xt1)
         X2=append!(X2,Xt2)
     end
     if length(taus)==(length(tausm)+1)
         Xtp=X1[length(X1)].-model.mp.ρ *X2[length(X2)] 
         Xts=append!(Xts,Xtp)
      end
     ti=sort([tausm;t])
     dfxt=DataFrame(Temps=ti, Degradation=Xts)
     return dfxt, taus, tausm, t, Xt , X1, X2
end


#= function ARD1(al,ale,mu1,mu2,sigma21,sigma22,r,rho,n::Int64,p::Float64,tau_periode::Int64,perio::Bool)
    nb_maint=floor(Int64,(n - 1) / (tau_periode + 1))
    n=n-nb_maint
     if perio
         t=0:p:n*p-1
     else
         Random.seed!(al)
         t=sort(sample(1:(n*5),n-1,replace=false))
         t=append!(zeros(1),t)
         Random.seed!()
     end
     taus=Float64[0;t[tau_periode+1:tau_periode:length(t)-1];t[length(t)]]
     
     tausm=t[range(tau_periode+1,length=nb_maint,step=tau_periode)]
     Random.seed!(ale)
     inc=[0 0]
     cov=r*sqrt(sigma21*sigma22)
     for i in 1:length(diff(t))
         inc=[inc;reshape(rand(MvNormal([mu1*diff(t)[i],mu2*diff(t)[i]],[sigma21*diff(t)[i] cov*diff(t)[i] ; cov*diff(t)[i] sigma22*diff(t)[i]])),1,2)]
     end
     Xt=cumsum(inc,dims=1)
     Xt0=Xt[findall(t.<=taus[2]),:]
     Xts=Float64[] ; X1=Float64[] ; X2=Float64[]
     for i in 1:(length(taus)-1)
         Xt1=Xt[intersect(findall(t.<=taus[i+1]),findall(t.>=taus[i])),1]
         Xt2=Xt[intersect(findall(t.<=taus[i+1]),findall(t.>=taus[i])),2]
         Xtf=Xt1.-rho*Xt2[1]
         Xts=append!(Xts,Xtf)
         X1=append!(X1,Xt1)
         X2=append!(X2,Xt2)
     end
     if length(taus)==(length(tausm)+1)
         Xtp=X1[length(X1)].-rho*X2[length(X2)] 
         Xts=append!(Xts,Xtp)
      end
     ti=sort([tausm;t])
     dfxt=DataFrame(Temps=ti, Degradation=Xts)
     return dfxt, taus, tausm, t, Xt , X1, X2
 end
 

 =#