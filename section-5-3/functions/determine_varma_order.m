% Determine the most concise VARMA representation possible. 
% All representations are of the form VARMA(kappa+2, kappa+1).
%
% Inputs:
%   M_        - Returned by Dynare
%   options_  - Returned by Dynare
%   oo_       - Returned by Dynare
%   max_kappa - indirectly determines the largest considered VARMA model. If 
%               none of the models up to the VARMA(max_kappa+2, max_kappa+1)
%               model match the DSGE, then no VARMA representation is returned.
%   obs_var   - (Optional) Vector of observed variables. Should be the indices
%               of those variabels in Dynare order.
%
% Outputs:
%   p         - Estimated VAR order of the VARMA model.
%   q         - Estimated MA order of the VARMA model.
%
% Notes:
%   - Based on Corollary 1 and 2 from the referenced paper, the function first checks 
%     if p and q can be immediately determined.
%   - If Corollary 1 and 2 do not apply, the function uses Proposition 1 to determine 
%     p and q iteratively.
%   - An error is raised if a VARMA representation cannot be found with the given max_kappa.
%
% References: 
%   - Morris, S. D. (2016). VARMA representation of DSGE models. 
%     Economics Letters, 138, 30–33. https://doi.org/10.1016/j.econlet.2015.11.027
function [p, q]=determine_varma_order(M_, options_, oo_, max_kappa, obs_var)
  if nargin==4
    [Ax, Ay, B, Cx, Cy, D, vars] = get_ABCD_varma(M_, options_, oo_);
  elseif nargin==5
    [Ax, Ay, B, Cx, Cy, D, vars] = get_ABCD_varma(M_, options_, oo_, obs_var);
  end

  % Main assumption
  if rank(D) < size(D, 1)
    error("D is not invertible.");
  end

  % Corollary 1 and 2
  m = size(Ax, 1);
  n = size(Cx, 1);
  if n == m && rank(Cx) == size(Cx, 1)
    p = 2;
    q = 1;
    if all(abs(Ay) < 1e-10) & all(abs(Cy) < 1e-10)
      p = 1;
      q = 1;
    end
    return; 
  end

  % Proposition 1
  p = Inf;
  q = Inf; 
  C = [Cx Cy];
  A = [[Ax Ay]; [Cx Cy]];
  F_kappa = Cx;
  for kappa=1:max_kappa
    F_kappa = [(C*A^(kappa-1)*[Ax' Cx']')' F_kappa'];
    F_kappa = F_kappa';
    F_kappa_plus = pinv(F_kappa);
    Psi_kappa = F_kappa_plus(:, 1:n);
    if rank(Psi_kappa) == size(Psi_kappa, 2)
      p = kappa + 2;
      q = kappa + 1;
      return;
    end
  end

  error("Could not find a VARMA representation with max_kappa="+max_kappa)
end




