abstract type AbstractMaintenancePolicy end

mutable struct PeriodicMaintenancePolicy <: AbstractMaintenancePolicy
ρ::Float64
period::Int64
end