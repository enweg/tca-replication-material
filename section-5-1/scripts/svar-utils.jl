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



function simulate_VAR(k, N, Bplus::AbstractMatrix{T}; constant = 0, linear_trend = 0, Σ=I, burnin=100) where {T}
    m = size(Bplus, 1)
    p = Int(m/k)
    y = Matrix{T}(undef, N+burnin+p, k)
    y[1:p, :] .= 0
    
    errors = rand(MultivariateNormal(zeros(k), Σ), N+burnin)
    
    for t = (p+1):(N+burnin+p)
        y[t, :] = constant .+ linear_trend*t .+ vec(y[(t-1):-1:(t-p), :]')'*Bplus + errors[:, t-p]'
    end
    return y[(p+burnin+1):end, :], errors[:, (burnin+1):end]'
end



function simulate_SVAR(k, N, A0::AbstractMatrix{T}, Aplus::AbstractMatrix{T}; 
    constant = 0, linear_trend = 0, shock_dist = MultivariateNormal(zeros(k), I), burnin = 100) where {T}

    A0inv = inv(A0)
    Bplus = Aplus * A0inv
    Σ = A0inv' * cov(shock_dist) * A0inv
    Σ = 0.5 * (Σ' + Σ)

    Y_sim, errors = simulate_VAR(k, N, Bplus; constant = constant, linear_trend = linear_trend, Σ = Σ, burnin = burnin)
    shocks = errors * A0

    return Y_sim, shocks
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
Make the VAR companion matrix

We assume that the SVAR is described by 
```math
y_t'A_0 = c_t + y_{t-1}'A_1 + ... + y_{t-p}'A_p + epsilon_t'. 
```
Thus, the SVAR can be written as 
```math
y_t'A_0 = x_t'Aplus + epsilon_t', 
```
where $x_t' = [c_t, y_{t-1}', ..., y_{t-p}']$ and $Aplus = [1'; A_1; A_2; ...; A_p]$, and where $Var(epsilon_t)=I$.

The reduced form VAR is given by 
```math
y_t' = x_t'Bplus + u_t'
```
where $Bplus = Aplus*A_{0}^{-1}$ and $u_t = epsilon_t'A_{0}^{-1}$, and where $Var(u_t) = (A_{0}'A_0)^{-1} = Σ$.
Taking the transpose on both sides, the more commonly used form emerges, 
```math
y_t = Bplus'x_t + u_t, 
```
which can also be written as
```math
z_t = c_t' + Cz_{t-1} + u_t, 
```
where $z_t = [y_t; y_{t-1}; ...; y_{t-p+1}]$ and where 
```math
C = \begin{bmatrix}
    A_1' & A_2' & A_3' & ... & A_{p_1}' & A_p' \\
    I    & O    & O    & ... & O        & O   \\
    O    & I    & O    & ... & O        & O   \\
    O    & O    & I    & ... & O        & O  \\
    ...  & ...  & ...  & ... & O        & O   \\
    O    & O    & O    & ... & I        & O   \\
\end{bmatrix}. 
```

The companion matrix is then given by $C$.
"""
function make_companion_matrix(Bplus, p)
    k = size(Bplus, 2)
    m = size(Bplus, 1)
    
    Bplus_transpose = Bplus'
    companion = diagm(ones(k*(p-1)))
    companion = hcat(zeros(k*(p-1), m - k*p), companion, zeros(k*(p-1), k))
    companion = vcat(Bplus_transpose, companion)
    return companion
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
function structural_irf_companion(A0, Aplus, p, horizons)
    k = size(A0, 1)
    m = size(Aplus, 1)
    n_ex = m - p*k
    A0inv = inv(A0)
    Bplus = Aplus * A0inv
    
    
    companion = make_companion_matrix(Bplus, p)
    irfs = [(companion^i)[1:k, 1:k] * A0inv' for i in filter(x -> !isinf(x), horizons)]
    if isinf(maximum(horizons))
        irf_inf = A0
        for i=1:p
            irf_inf -= Aplus[(n_ex+(i-1)*k+1):(n_ex+i*k), :]
        end
        push!(irfs, inv(irf_inf)')
    end     
            
    return cat(irfs...; dims = 3)
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




function is_stable_SVAR(A0, Aplus, p; verbose = false)
    Bplus = Aplus*inv(A0)
    companion = make_companion_matrix(Bplus, p)
    maxabs_eigenvalue = maximum(abs.(eigen(companion).values))
    if verbose
        println("Maximum Absolute Eigenvalue = $maxabs_eigenvalue")
    end
    return maxabs_eigenvalue < 1
end



function diff(x)
    return vcat(missing, x[2:end] - x[1:(end - 1)])
end


"""
Assumes y_t'A_0 = x_t'Aplus + epsilon_t'
"""
function to_transition_matrices(A0, Aplus, p, max_horizon)
    A0_wegner = A0'
    Aplus_wegner = Aplus'
    Σ = inv(A0_wegner) * inv(A0_wegner)'
    Σ = Hermitian(0.5 * (Σ' + Σ))
    L = inv(cholesky(Σ).L)
    Q = A0_wegner * inv(L)
    D = inv(diagm(diag(L)))
    DL = D*L
    DQt = D*Q'
    
    k = size(A0, 1)
    
    B = reduce(hcat, [DQt*Aplus_wegner[:, ((i-1)*k+1):(i*k)] for i=p:-1:1])
    B = hcat(B, I-DL)
    Bs = [hcat(
            zeros(k, max(0, (i+1)*k - k*(p+1))),
            B[:, max(1, (end-(i+1)*k+1)):end], 
            zeros(k, (max_horizon+1)*k - (i+1)*k)
          ) for i in 0:max_horizon]
    B = reduce(vcat, Bs)
    
    eye_h = diagm(ones(max_horizon+1))
    Qbb = kron(eye_h, DQt)
    
    return B, Qbb
end


function structural_irf_transmission(A0, Aplus, p, horizons)
    B, Qbb = to_transition_matrices(A0, Aplus, p, maximum(horizons))
    k = size(A0, 1)
    return (inv(I - B)*Qbb)[:, 1:k]
end


function get_random_svar_params(k::Int, p::Int, persistence::Float64)
    # A0 = 10*randn(k, k)
    A0 = 20*rand(k, k) .- 10
    m = k*p
    Aplus = rand(m, k)
    Bplus = Aplus * inv(A0)
    companion = make_companion_matrix(Bplus, p)
    max_abs_eig = maximum(abs.(eigen(companion).values))
    while (max_abs_eig > persistence)
        Bplus .*= 0.99
        companion = make_companion_matrix(Bplus, p)
        max_abs_eig = maximum(abs.(eigen(companion).values))
    end
    Aplus = Bplus*A0
    return A0, Aplus, max_abs_eig
end

"""

y_t'A_0 = ...

A0_ortho and Aplus_ortho are obtained via a cholesky decomposition
"""
function make_structural_B(A0_ortho, Aplus_ortho, p, max_horizon)
    L = A0_ortho'
    D = inv(diagm(diag(L)))
    DL = D*L
    B = D*Aplus_ortho'
    
    B = reduce(hcat, [B[:, ((i-1)*k+1):(i*k)] for i=p:-1:1])
    B = hcat(B, I-DL)
    Bs = [hcat(
            zeros(k, max(0, (i+1)*k - k*(p+1))),
            B[:, max(1, (end-(i+1)*k+1)):end], 
            zeros(k, (max_horizon+1)*k - (i+1)*k)
          ) for i in 0:max_horizon]
    B = reduce(vcat, Bs)
    return B
end

"""
y_t'A_0 = ...

A0_ortho and Aplus_ortho are obtained via a cholesky decomposition
irf0 is the impact impulse response of all variable; it is k times k with k being the number of variables
in the system. If a shock has not been identified, its column is filled with NaN.
"""
function make_structural_Qbb(irf0, A0_ortho, max_horizon)
    
    m, k = size(A0_ortho)
    
    L = A0_ortho'
    D = inv(diagm(diag(L)))
    
    Qt = L*irf0
    DQt = D*Qt
    eye_h = diagm(ones(max_horizon+1))
    Qbb = kron(eye_h, DQt)
    
    return Qbb
end

"""
y_t'A_0 = ...

A0_ortho and Aplus_ortho are obtained via a cholesky decomposition
irf0 is the impact impulse response of all variable; it is k times k with k being the number of variables
in the system. If a shock has not been identified, its column is filled with NaN.
"""
function to_structural_transmission_model(irf0, A0_ortho, Aplus_ortho, p, max_horizon)
    B = make_structural_B(A0_ortho, Aplus_ortho, p, max_horizon)
    Qbb = make_structural_Qbb(irf0, A0_ortho, max_horizon)
    return B, Qbb
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
