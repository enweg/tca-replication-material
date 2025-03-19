"""
    internal_instrument_SVAR(data, p, horizons; kwargs)

Estimates an SVAR using an internal instrument. Assumes the first structural 
shock is being instrumented. 

## Arguments
- `data::Matrix`: Data matrix for the estimation of the SVAR. Assumes that the
  instrument is in the first column and the normalising variable is in the 
  second column. 
- `p::Int`: Lag-length of the SVAR. 
- `horizons::Int`: Number of horizons for which the IRFs should be computed. 
   Contemporaneous horizon equals 0. 

## Keyword Arguments
- `kwargs`: Passed on to `create_XY`.
"""
function internal_instrument_SVAR(data, p, horizons; kwargs...)
    X, Y = create_XY(data, p; kwargs...)
    B_plus = OLS_Bplus(X, Y)
    errors = Y - X*B_plus
    Σ = cov(errors)
    Σ = 0.5 * (Σ' + Σ)
    A0, Aplus = reduced_to_structural(B_plus, Σ, I)
    irfs = structural_irf(A0, Aplus, p, vcat(0, horizons))
    irfs = irfs ./ irfs[2, 1, 1]
    irfs = irfs[2:end, 1:1, 2:end]
    return irfs
end
