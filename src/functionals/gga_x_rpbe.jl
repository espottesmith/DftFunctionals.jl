struct RpbeExchange{CA} <: Functional{:gga,:x} where {CA<:ComponentArray{<:Number}}
    parameters::CA
    identifier::Symbol
end
function RpbeExchange(parameters::ComponentArray)
    RpbeExchange(parameters, :gga_x_rpbe_custom)
end

identifier(rpbe::RpbeExchange) = rpbe.identifier
parameters(rpbe::RpbeExchange) = rpbe.parameters
function change_parameters(rpbe::RpbeExchange, parameters::ComponentArray;
                           keep_identifier=false)
    if keep_identifier
        RpbeExchange(parameters, rpbe.identifier)
    else
        RpbeExchange(parameters)
    end
end

function energy(rpbe::RpbeExchange, ρ::T, σ::U) where {T<:Number,U<:Number}
    TT = arithmetic_type(rpbe, T, U)

    # TODO This function is quite sensitive to the floating-point type ...
    #      so for now we don't bother doing this in TT, but rather convert before return
    
    κ = rpbe.parameters.κ
    μ = rpbe.parameters.μ

    rpbe_x_f(s²) = 1 + κ * (1 - exp(-μ * s²^2 / κ))   # (RPBE, Eq. 15)
    # rₛ = cbrt(3 / (4π  * ρ))                 # page 2, left column, top
    # kF = cbrt(3π^2 * ρ)                      # page 2, left column, top
    # s  = sqrt(σ) / (2kF * ρ)                 # below (9)
    s² = σ / (ρ^(4 / 3) * 2cbrt(3π^2))^2

    res = energy(LdaExchange(), ρ) * rpbe_x_f(s²)     # (10)
    TT(res)
end

#
# Concrete functionals
#

"""
RPBE exchange.
Hammer, Hansen, Nørskov 1998 (DOI 10.1103/PhysRevB.59.7413)
"""
function DftFunctional(::Val{:gga_x_rpbe})
    RpbeExchange(ComponentArray(κ=0.8040, μ=pbe_μ_from_β(0.06672455060314922)), :gga_x_rpbe)
end