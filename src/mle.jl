function setΔyj(df, nb_maint, τ_periode)
    i, im = 1, 1
    yAv = 0.0
    Δyj = Float64[]
    while im <= nb_maint
        i += τ_periode
        yAp = df[!,:Degradation][i]
        push!(Δyj, yAp - yAv)
        i += 1
        yAv = df[!,:Degradation][i]
        im += 1
    end 
    Δyj
end

function manip_df(df::DataFrame)
  Δt = diff(df[!,:Time]) #Delta t_{j,i} et temps des sauts
  Δy = diff(df[!,:Degradation]) #Delta y_{j,i} et sauts
  δt = Δt[1]    # At least 1 observation between maintenances
  maints = Δt .== 0 # BitVector
  nb_maint=sum(maints) #si perio, ne marche qu avec un operateur comme + comprends pas pq...
  τ_periode= findfirst(maints) - 1
  zj=Δy[maints]  #juste les sauts
  Δτj=repeat([δt * τ_periode],nb_maint) 
  Δyj=setΔyj(df, nb_maint, τ_periode)
  Δtji=Δt[.!maints]
  Δyji=Δy[.!maints]
  (zj=zj, Δτj=Δτj, Δyj=Δyj, Δtji=Δtji, Δyji=Δyji, δt=δt)
end


# manip_df<-function(df){
#   Dt=diff(df[,"Temps"]) #Delta t_{j,i} et temps des sauts
#   Dy=diff(df[,"Degradation"]) #Delta y_{j,i} et sauts
#   df2=data.frame(DTemps=Dt,DDegradation=Dy)
#   nb_maint=sum(Dt==0)
#   tau_periode=(which(Dt==0)[1])-1
#   df2_NOj=df2[df2[,1]!=0,] #sans les sauts : Delta t_{j,i} et  Delta y_{j,i}
#   djumps=df2[df2[,1]==0,2]  #juste les sauts
#   Mrt=matrix(df2_NOj[1:(tau_periode*nb_maint),1],tau_periode,nb_maint)
#   Dtauj=apply(Mrt,2,sum) # Delta tau_j + on convertit en vecteur
#   Mrt2=matrix(df2_NOj[1:(tau_periode*nb_maint),2],tau_periode,nb_maint)
#   Dyj=apply(Mrt2,2,sum) # Delta y_j
#   dtji=df2_NOj[,1] ; dyji=df2_NOj[,2]
#   list(Dt=Dt, Dy=Dy, df2=df2, nb_maint=nb_maint, tau_periode=tau_periode, 
#        df2_NOj=df2_NOj, djumps=djumps, Dtauj=Dtauj, Dyj=Dyj, dtji=dtji, dyji=dyji)
# }

function mle(x,df)
    mdf=manip_df(df)
    μ1, μ2, σ21, σ22, ρ, cov = x
    cov *= sqrt(σ21 * σ22)
    ## zj, Δτj, Δyj, Δtji, Δyji, δt = mdf.zj, mdf.Δτj, mdf.Δyj, mdf.Δtji, mdf.Δyji, mdf.δt
    
    logvrais = 0.5 * (
        sum(log.(2 .* π .* σ21 .* mdf.Δtji) .+ ((mdf.Δyji .- μ1 .* mdf.Δtji).^ 2) ./ (σ21 .* mdf.Δtji))
        + sum( log.(2 .* π .* ρ .^ 2 .* mdf.Δτj .* (σ22 .- (cov .^ 2 ./ σ21)))
                .+ ((mdf.zj .+ ρ .* (μ2 .* mdf.Δτj .+ (cov ./ σ21) .* (mdf.Δyj .- μ1 .* mdf.Δτj))) .^ 2) ./ (ρ .^ 2 .* mdf.Δτj .* (σ22 .- (cov .^ 2 ./ σ21)))) 
        )
    logvrais
end

# logvrais=0.5*(sum(log(2*pi *sigma21 *dtji )+ ((dyji-mu1*dtji)^2)/(sigma21*dtji))
#              +sum(log(2*pi *rho^2 *Dtauj*(sigma22-(cov^2/sigma21)))
#                   +((djumps+rho*(mu2*Dtauj+(cov/sigma21)*(Dyj-mu1*Dtauj)))^2)/(rho^2*Dtauj*(sigma22-(cov^2/sigma21)))))
  

# manip_df4<-function(df){
#   Dt=diff(df[,"Temps"]) #Delta t_{j,i} et temps des sauts
#   Dy=diff(df[,"Degradation"]) #Delta y_{j,i} et sauts
#   p<-Dt[1]
#   df2=data.frame(DTemps=Dt,DDegradation=Dy)
#   nb_maint=sum((Dt)==(max(Dt))) #si perio, ne marche qu avec un operateur comme + comprends pas pq...
#   tau_periode=which((Dt)==(max(Dt)))[1]
#   df2_NOj=df2[(df2[,1])!=(max(Dt)),] #sans les sauts : Delta t_{j,i} et  Delta y_{j,i}
#   zj=df2[(df2[,1])==(max(Dt)),2]  #juste les sauts
#   Dtauj=rep(p*tau_periode,nb_maint) 
#   Mrt1=matrix(df2_NOj[1:(tau_periode-1),2],tau_periode-1,1)# \Delta y_{0,i}
#   Mrt2=matrix(df2_NOj[seq(tau_periode,length.out=(nb_maint-1)*(tau_periode-2)),2],tau_periode-2,nb_maint-1)
#   Δyj=c(sum(Mrt1),apply(Mrt2,2,sum)) # Delta y_j
#   Δtji=df2_NOj[,1] ; Δyji=df2_NOj[,2]
#   list(Dt=Dt, Dy=Dy, df2=df2, nb_maint=nb_maint, tau_periode=tau_periode, 
#        df2_NOj=df2_NOj, zj=zj, Dtauj=Dtauj, Δyj=Δyj, Δtji=Δtji, Δyji=Δyji, p=p)
# }



# LV4<-function(x,df){
#   mdf=manip_df4(df)
#   μ1=x[[1]] ; μ2=x[[2]] ; σ21=x[[3]] ; σ22=x[[4]] ; ρ=x[[5]] ; cov=x[[6]]*sqrt(σ21*σ22)
#   zj=mdf[[7]]; Dtauj=mdf[[8]]; Δyj=mdf[[9]]; Δtji=mdf[[10]]; Δyji=mdf[[11]] ; p=mdf[[12]]
  
#   nb_maint<-length(zj)
#   pvec<-rep(p,nb_maint)
#   Dt_obs<-c(Dtauj[1]-p,Dtauj[-1]-2*p)
  
#   μ_zj<-μ1*2*pvec-ρ*μ2*Dtauj #length=nb_maint
#   σ2_zj<-σ21*2*pvec+ρ^2*σ22*Dtauj-2*ρ*cov*pvec
  

#   logvrais=(1/2)*(sum(log(2*pi *σ21 *Δtji )+((Δyji-μ1*Δtji)^2)/(σ21*Δtji))
#              +sum(log(2*pi *(σ2_zj-ρ^2 *cov^2*((Dt_obs/σ21)+(c(0,(pvec[-1])^2)/c(1,σ2_zj[-(length(σ2_zj))]))))) #1 au pif
#                   +(zj-μ_zj+((ρ*cov)/σ21)*(Δyj-μ1*Dt_obs)
#                  +(ρ*cov*c(0,pvec[-1])/c(1,σ2_zj[-(length(σ2_zj))]))*(c(0,zj[-length(zj)])-c(0,μ_zj[-length(μ_zj)])))^2/(σ2_zj-ρ^2 *cov^2*((Dt_obs/σ21)+((c(0,pvec[-1])^2)/c(1,σ2_zj[-length(σ2_zj)]))))))
#   return(logvrais)
# } 