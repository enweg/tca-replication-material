

"""

Assumes that the first `num_instruments` columns in `data` are the instruments and that the `num_instruments+1`
column is the to instrumenting variable. Also assumes the structural shock being instrumented is the first.
"""
function external_instrument_SVAR(data, instruments, p, horizons; kwargs...)
    
    X, Y = create_XY(data, p; kwargs...)
    Bplus = OLS_Bplus(X, Y)
    errors = Y - X*Bplus
    Σ = cov(errors)
    Σ = Hermitian(0.5 * (Σ' + Σ))
    A0, Aplus = reduced_to_structural(Bplus, Σ, I)
    irfs = structural_irf(A0, Aplus, p, horizons)
    reduced_irfs = mapslices(x -> x*A0', irfs; dims = [1, 2])
        
    Z = instruments
    missing_rows = [any(ismissing.(r)) for r in eachrow(Z)]
    not_missing = map(x -> !x, missing_rows)
    @info "Using $(sum(not_missing)) observations for 2SLS. \n$(sum(missing_rows)) missing values were removed."
    X, Y = create_XY(data[not_missing, :], p; kwargs...)
    Z = Z[not_missing, :]
    X2sls = hcat(Z[(p+1):end, :], X)
    first_stage = (X2sls'*X2sls) \ (X2sls' * Y[:, 1])
    first_stage_fitted = X2sls*first_stage
    
    X2sls = hcat(first_stage_fitted, X)
    second_stage = (X2sls' * X2sls) \ (X2sls' * Y)
    theta_0i1 = second_stage[1, :]
    # return theta_0i1
        
    relative_irfs = mapslices(x -> x*theta_0i1, reduced_irfs; dims = [1,2])
    
    return relative_irfs

end