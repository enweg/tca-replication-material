function lag_matrix(mat::AbstractMatrix, lag::Int; pad_missing = true)
    missings = fill(missing, lag, size(mat, 2))
    lag_mat = mat[1:(end-lag), :]
    !pad_missing && return lag_mat
    return vcat(missings, mat[1:(end-lag), :])
end

function lag_matrix(mat::AbstractMatrix, lags::AbstractVector{Int}; pad_missing = true)
    lag_mats = [lag_matrix(mat, l; pad_missing = pad_missing) for l in lags]
    return lag_mats
end

function create_XY(data_mat::AbstractMatrix{T}, p; include_constant = true, include_linear_trend = false) where {T}
    lag_mats = reduce(hcat, lag_matrix(data_mat, 1:p; pad_missing = true))
    constant = ones(size(lag_mats, 1))
    linear_trend = 1:size(lag_mats, 1)
    X = lag_mats
    if include_linear_trend 
        X = hcat(linear_trend, X)
    end
    if include_constant
        X = hcat(constant, X)
    end
    X = X[(p+1):end, :]
    Y = data_mat[(p+1):end, :]
    return T.(X), Y
end

function OLS_Bplus(X::AbstractMatrix{T}, Y::AbstractMatrix{T}) where {T}
    Bplus = (X'*X) \ (X'*Y)
    return Bplus
end

@doc raw"""
# Move from reduced form VAR coefficients to structural coeffcients

We assume that the SVAR is described by 

```math
y_t'A_0 = c_t + y_{t-1}'A_1 + ... + y_{t-p}'A_p + \varepsilon_t'. 
```
Thus, the SVAR can be written as 
```math
y_t'A_0 = x_t'A^+ + \varepsilon_t', 
```
where $x_t' = [c_t, y_{t-1}', ..., y_{t-p}']$ and $A^+ = [1'; A_1; A_2; ...; A_p]$, and where $Var(\varepsilon_t)=I$. 
Using the QL-decomposition of $A_0'$, we have $A_0 = RQ'$ where $R$ is an upper-triangular matrix and $Q$
is an orthogonal matrix. The reduced form VAR is given by 
```math
y_t' = x_t'B^+ + u_t'
```
where $B^+ = A^+A_{0}^{-1}$ and $u_t' = \varepsilon_t'A_{0}^{-1}$ and where $Var(u_t) = (A_{0}A_0')^{-1} = Σ$. 
Thus, the mapping from reduced form to structural coefficients is given by 
```math
A^+ = B^+A_0 = B^+RQ'
```
```math
A0 = inv(cholesky(Σ).R)Q'
```

## Arguments

- `Bplus`: Reduced form coefficient matrix.
- `Σ`: Reduced form error covariance.
- `Q`: Orthogonal matrix. 

## Returns

- `A0`: Structural contemporaneous matrix
- `Aplus`: Structural lagged coefficient matix

## Examples

### Example: Recursive SVAR

```
k = 3
p = 4
k_ex = 1
m = p*k + k_ex
A0 = UpperTriangular(rand(k, k))
Aplus = randn(m, k)
A0inv = inv(A0)
Bplus = Aplus * A0inv
Σ = A0inv' * A0inv
Q = I

A0_reconstructed, Aplus_reconstructed = reduced_to_structural(Bplus, Σ, Q)
```

### Example: Non-recursive SVAR

```
k = 3
p = 4
k_ex = 1
m = p*k + k_ex
A0 = randn(k, k)
P = diagm(ones(1:k))[k:-1:1, :]
F = qr(A0'*P)
L = P*F.R*P
Q = F.Q*P
S = diagm(sign.(diag(L)))
Q = Q*S
L = S*L



Aplus = randn(m, k)
A0inv = inv(A0)
Bplus = Aplus * A0inv
Σ = A0inv' * A0inv

A0_reconstructed, Aplus_reconstructed = reduced_to_structural(Bplus, Σ, Q)
```


"""
function reduced_to_structural(Bplus, Σ, Q)

    Σ = Hermitian(0.5 * (Σ' + Σ))

    R = cholesky(Σ).U
    # normalising diagonal to positive numbers
    for i in 1:size(R, 1)
        if R[i, i] < 0
            R[i, :] *= -1
        end
    end
    R = inv(R)
    A0 = R*Q'
    Aplus = Bplus * A0
    return A0, Aplus
end

@doc raw"""

We assume that the SVAR is described by 
```math
y_t'A_0 = c_t + y_{t-1}'A_1 + ... + y_{t-p}'A_p + epsilon_t'. 
```
Thus, the SVAR can be written as 
```math
y_t'A_0 = x_t'Aplus + epsilon_t', 
```
where $x_t' = [c_t, y_{t-1}', ..., y_{t-p}']$ and $Aplus = [1'; A_1; A_2; ...; A_p]$, and where $Var(epsilon_t)=I$.


## Returns

Returns a three dimensional array with the first dimension corresponding to the endogenous variables, 
the second dimension corresponding to the structural shocks, and the third dimension corresponding to the 
horizon in order of the `horizons` vector. 
"""
function structural_irf(A0::AbstractMatrix{T}, Aplus, p, horizons) where {T}
    k = size(A0, 1)
    m = size(Aplus, 1)
    n_ex = m - p*k
    A0inv = inv(A0)
    Bplus = Aplus * A0inv
    
    max_horizon = Int(maximum(filter(x -> !isinf(x), horizons)))
    R = Vector{Matrix{T}}(undef, max_horizon + 1)
    R[1] = A0inv'
    for i=1:max_horizon
        R[i+1] = zeros(k, k)
        for j = 1:minimum([i, p])
            Aj = Aplus[(n_ex + (j-1)*k + 1):(n_ex + j*k), :]'
            R[i+1] .+= Aj * R[i+1-j]
        end
        R[i+1] = R[1] * R[i+1]
    end
    if isinf(maximum(horizons))
        Rinf = A0'
        for j = 1:p
            Aj = Aplus[(n_ex + (j-1)*k + 1):(n_ex + j*k), :]'
            Rinf .-= Aj
        end
        Rinf = inv(Rinf)
    end
    
    Rout = R[Int.(filter(x->!isinf(x), horizons)) .+ 1]
    if isinf(maximum(horizons))
        push!(Rout, Rinf)
    end
    return cat(Rout...; dims = 3)
end 

function diff(x)
    return vcat(missing, x[2:end] - x[1:(end - 1)])
end

"""
Cholesky identified IRFs.
"""
function orthogonal_irfs(data, p, horizons; kwargs...)
    X, Y = create_XY(data, p; kwargs...)
    Bplus = OLS_Bplus(X, Y)
    errors = Y - X*Bplus
    Σ = cov(errors)
    Σ = Hermitian(0.5 * (Σ' + Σ))
    A0, Aplus = reduced_to_structural(Bplus, Σ, I)
    irfs = structural_irf(A0, Aplus, p, horizons)
    return irfs
end
