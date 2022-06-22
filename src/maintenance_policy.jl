abstract type AbstractMaintenancePolicy end

mutable struct PeriodicMaintenancePolicy <: AbstractMaintenancePolicy
ρ::Float64
τ::Float64
end