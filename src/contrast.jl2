function manip_df(df)
    Dt=diff(df[:,"Temps"]) #Δt_{j,i} et temps des sauts
    Dy=diff(df[:,"Degradation"]) #Δy_{j,i} et sauts
    df2=DataFrame(DTemps=Dt,DDegradation=Dy)
    nb_maint=sum(Dt.==0)
    tau_periode=(findall(Dt.==0)[1])-1
    df2_NOj=df2[df2[:,1].!=0,:] #sans les sauts : Δt_{j,i} et  Δy_{j,i}
    djumps=df2[df2[:,1].==0,2]  #juste les sauts
    #last_df2=last(df2_NOj,(size(Dt)[1]-nb_maint)%nb_maint) #valeurs après dernières maintenance
    Mrt=reshape(df2_NOj[1:tau_periode*nb_maint,1],(nb_maint,tau_periode))
    Dtauj=vec(sum(Mrt,dims=2)) # Δτ_j + on convertit en vecteur
    Mrt2=reshape(df2_NOj[1:tau_periode*nb_maint,2],(tau_periode,nb_maint))
    Dyj=vec(sum(Mrt2,dims=1)) # Δy_j
    dtji=df2_NOj[:,1] ; dyji=df2_NOj[:,2]
    return Dt, Dy, df2, nb_maint, tau_periode, df2_NOj, djumps, Mrt, Dtauj, Mrt2, Dyj, dtji, dyji
end

function LV1bis(x::Vector)  
    mu1=x[1] ; mu2=x[2] ; sigma21=x[3] ; sigma22=x[4] ; rho=x[5] ; r=x[6]
    Dt, Dy, df2, nb_maint, tau_periode, df2_NOj, djumps, Mrt, Dtauj, Mrt2, Dyj, dtji, dyji = manip_df(df)
    cov=r*sqrt(sigma21*sigma22)
    logvrais=-(sum(-log.(2π *sigma21 *dtji )/2- (1/2)*((dyji-mu1*dtji).^2)./(sigma21*dtji))
                -sum(log.(2π *rho^2 *Dtauj*(sigma22-(cov^2/sigma21)))/2
                +(1/2)*((djumps+rho*(mu2*Dtauj+cov*(sigma22/sigma21)*(Dyj-mu1*Dtauj))).^2)./(rho^2*Dtauj*(sigma22-(cov^2/sigma21)))))
                
            return logvrais
end


function mle()
    nb_traj=100
    MLE=Array{Float64}(undef, 101, 6)
    j=1
    for k in 200:200
        df=ARD1(24,k,2,3,4,5,0.8,0.5,20,1,3,true)[1]  
        lower = [-Inf,-Inf,0.01,0.01,0.01,-0.99]
        upper = [Inf, Inf,Inf,Inf,0.99,0.99]
        initial_x = [2,3,4,5,0.5,0.7]
        #inner_optimizer =  LBFGS()
        inner_optimizer = NelderMead() 
        results = optimize(LV1bis,lower, upper, initial_x, Fminbox(inner_optimizer))
        res2=Optim.minimizer(results)
        MLE[j,:]=res2
        j=j+1
    end
    return MLE
end

