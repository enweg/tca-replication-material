"""

Assumes instrument is in first column and to instrumenting variable is in second column. Assumes first structural 
shock is being instrumented. 
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