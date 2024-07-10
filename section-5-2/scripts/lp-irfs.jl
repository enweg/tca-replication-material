using LinearAlgebra

function make_lag_matrix!(X::AbstractMatrix{T}, Y::AbstractMatrix{T}) where {T<:Real}
  view(X, :, :) .= T(NaN)
  _, k = size(Y)
  nlag = floor(Int, size(X, 2) / k)
  for l=1:nlag
    @views X[(l+1):end, ((l-1)*k+1):(l*k)] .= Y[1:(end-l), :]
  end
  return X
end
function make_lag_matrix(Y::AbstractMatrix{T}, nlag::Int) where {T<:Real}
  r, c = size(Y)
  X = zeros(T, r, c*nlag)
  return make_lag_matrix!(X, Y)
end

function _lp_irf(
    data::DataFrame, 
    shock::Symbol, 
    contemp_controls::Vector{Symbol}, 
    p::Int; 
    include_constant::Bool=true
)

    Y = Matrix(data)
    X = make_lag_matrix(Y, p)    
    X = X[(p+1):end, :]
    Y = Y[(p+1):end, :]
    ix_contemp_controls = map(x -> findfirst(==(string(x)), names(data)), contemp_controls)
    ix_shock = findfirst(==(string(shock)), names(data))
    
    # adding contemporaneous controls
    X = hcat(Y[:, ix_contemp_controls], X)
    # adding constant if wanted
    if include_constant
        X = hcat(ones(size(X, 1)), X)
    end
    # finally include contemp part of shock variable
    X = hcat(Y[:, ix_shock], X)

    k = size(Y, 2)
    irfs = zeros(k, 1, length(horizons))
    for (i, h) in enumerate(horizons)
        Y_h = Y[(h+1):end, :]
        X_h = X[1:(end-h), :]
        beta_h = X_h' * X_h \ X_h'*Y_h
        irfs[:, 1, i] .= beta_h[1, :]
    end

    return irfs, names(data)
end
