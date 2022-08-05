abstract type AbstractMaintenancePolicy end

mutable struct PeriodicMaintenancePolicy <: AbstractMaintenancePolicy
ρ::Float64
τ_periode::Float64
end