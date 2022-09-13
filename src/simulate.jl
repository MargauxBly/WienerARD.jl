function simulate(model::Model,n::Int64,period::Float64;δ=1.0)::DataFrame
            τ=model.mp.period*δ

            Δt=period*δ
            nb_maint=floor((n - 1) / (model.mp.period + 1))
            n_real=n-Int(nb_maint)
            time=range(start=0,length=n_real,step=p)
            Xt=zeros(n_real,2)
            Yt=zeros(n)
            j , t = 1, time[1]  
            Xt1, Xt2 = 0.0 , 0.0      
            for i in 2:n_real
                t=time[i]
                ΔX12 = ΔX(model,Δt)
                Xt1 += ΔX12[1]
                Xt2 += ΔX12[2] 
                Xt[i,1]=Xt1
                Xt[i,2]=Xt2
                Yt[j]=Xt1
                if  τ == t #maintenance
                    τ += τ
                    j+=1
                    Yt[j]= Xt1-model.mp.ρ*(Xt2-Xt[i-(model.mp.period+1),2])
                end
                j+=1
                
            end
        
        
        return DataFrame(Time=time,Degradation=Yt)
    end
    







# function simulate3(model::Model,n::Int64,p::Int64;δ=1.0)::DataFrame
#      τ_mp=model.mp.τ
#      Δt=p*δ

#      Xt=zeros(n-nb_maint,2) ; Yt=zeros(n,1)
    
#      Xt1, Xt2=0.0, 0.0
#      t=0.0
#      time=zeros(n)
#      nb_maint=floor(((n-1)*Δt)/τ_mp)
#      ind_τ=[0;(ceil(τ_mp/Δt)+1):(ceil(τ_mp/Δt)):(n-nb_maint)] #indices des temps de maintenance sur Xt
#      N=n-nb_maint
#      i=2
#      for l in 2:N     
#         Xt[l,1]=Xt1
#         Xt[l,2]=Xt2

#         if  τ_mp < t #min temps
#             Δt1=Δt-(t-τ_mp)
#             Δt2=t-τ_mp
#             time[i]=τ_mp       
#             Δy = ΔX(model,Δt1)
#             Yt[i] = Xt1      
#             t = τ_mp
            
            
#             i += 1
#             if τ_mp == t
#                 time[i]=τ_mp
#                 t += Δt2
#                 Yt[i] = Xt1 - model.mp.ρ*(Xt2-Xt[ind_τ[k],2])
#                 Δy = ΔX(model,Δt2)  
#                 i += 1
#             end 
#             τ_mp += τ_mp   
             
#         else
#             if τ_mp == t
#                 time[i]=τ_mp
#                 Yt[i] = Xt1 - model.mp.ρ*(Xt2-Xt[ind_τ[k],2]) 
#                 i += 1
#                 τ_mp += τ_mp  
#             end
#             time[i]=t
#             t += Δt     
#             Δy = ΔX(model,Δt)   
#             Xt1 += Δy[1]
#             Xt2 += Δy[2] 
#             Yt[i]=Xt1       
#             i += 1    
#         end
#         end

#         Xt1 += Δy[1]
#         Xt2 += Δy[2] 
        
        
        
#      end
     
#     return DataFrame(Time=time,Degradation=Yt)
# end





# function simulate2(model::Model,n::Int64,p::Int64)::DataFrame
#     tau_periode=model.mp.τ
#     nb_maint=floor(Int64,(n - 1) / (tau_periode + 1))
#     n=n-nb_maint
#     t=0:p:n*p-1
   
#      taus=Float64[0;t[tau_periode+1:tau_periode:length(t)-1];t[length(t)]]
     
#      tausm=t[range(tau_periode+1,length=nb_maint,step=tau_periode)]
#      Random.seed!(ale)
#      inc=[0 0]
#      cov=model.r*sqrt(model.sigma21*model.sigma22)
#      for i in 1:length(diff(t))
#          inc=[inc;reshape(rand(MvNormal([model.mu1*diff(t)[i],mu2*diff(t)[i]],[model.sigma21*diff(t)[i] cov*diff(t)[i] ; cov*diff(t)[i] sigma22*diff(t)[i]])),1,2)]
#      end
#      Xt=cumsum(inc,dims=1)
#      Xt0=Xt[findall(t.<=taus[2]),:]
#      Xts=Float64[] ; X1=Float64[] ; X2=Float64[]
#      for i in 1:(length(taus)-1)
#          Xt1=Xt[intersect(findall(t.<=taus[i+1]),findall(t.>=taus[i])),1]
#          Xt2=Xt[intersect(findall(t.<=taus[i+1]),findall(t.>=taus[i])),2]
#          Xtf=Xt1.-model.mp.ρ*Xt2[1]
#          Xts=append!(Xts,Xtf)
#          X1=append!(X1,Xt1)
#          X2=append!(X2,Xt2)
#      end
#      if length(taus)==(length(tausm)+1)
#          Xtp=X1[length(X1)].-model.mp.ρ *X2[length(X2)] 
#          Xts=append!(Xts,Xtp)
#       end
#      ti=sort([tausm;t])
#      dfxt=DataFrame(Temps=ti, Degradation=Xts)
#      return dfxt, taus, tausm, t, Xt , X1, X2
# end


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

