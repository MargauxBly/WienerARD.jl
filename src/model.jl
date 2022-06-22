mutable struct Model 
    μ1::Float64
    μ2::Float64
    σ21::Float64
    σ22::Float64
    r::Float64
    mp::PeriodicMaintenancePolicy
end

function ΔX(model::Model, Δt::Float64)::Vector{Float64}
    cov=model.r*sqrt(model.σ21*model.σ22)
    return rand(MvNormal([model.μ1*Δt , model.μ2*Δt],[model.σ21*Δt cov*Δt ; cov*Δt model.σ22*Δt]))
end
