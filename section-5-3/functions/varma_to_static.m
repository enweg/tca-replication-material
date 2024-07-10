% Convert a VARMA into the static form x = B*x + Omega*eps used in Wegner etal (2024).
%
% Inputs:
%   vars    - Vector of variable IDs.
%   A0      - Initial impact matrix of size k x k.
%   Phis    - VAR coefficient matrices of size k x k x p.
%   Psis    - MA coefficient matrices of size k x k x q.
%   horizon - Integer specifying the maximum horizon.
%   T       - (Optional) Transformation matrix for reordering variables.
%             Corresponds to the transmission matrix in Wegner etal (2024)
%
% Outputs:
%   B       - See Wegner etal (2024)
%   Oomega  - See Wegner etal (2024)
%   k       - Number of variables in the system.
%   vars    - Reordered vector of variable IDs.
%
% Notes:
%           - The VARMA coefficients can be obtained from get_varma_coeffs
%
% References: 
%   - Wegner, E., Lieb, L., Smeekes, S., & Wilms, I. (2024). 
%     Transmission Channel Analysis in Dynamic Models. 
%     arXiv preprint arXiv:2405.18987.
%
function [B, Oomega, k, vars]=varma_to_static(vars, A0, Phis, Psis, horizon, T)
  h = horizon;
  k = size(A0, 1);
  p = size(Phis, 3);
  q = size(Psis, 3);

  if nargin > 5
    warning("Reordering variables")
    A0 = A0 * T';
    for i=1:p
      Phis(:, :, i) = Phis(:, :, i) * T';
    end
    vars = vars(T * vec(1:k));
  end

  [Q, L] = ql(A0);
  D = diag(diag(L).^(-1));
  DL = D*L;
  DQt = D * Q';

  Oomega = zeros(k * (h+1), k * (h+1));
  B = zeros(k * (h+1), k * (h+1));
  for c=1:(h+1)
    for r=c:(h+1)
      if r == c
        B((k*(r-1)+1):(r*k), (k*(r-1)+1):(r*k)) = eye(k) - DL;
        Oomega((k*(r-1)+1):(r*k), (k*(r-1)+1):(r*k)) = DQt;
      else
        i = r - c;
        if i <= p
          B((k*(r-1)+1):(k*r), (k*(c-1)+1):(k*c)) = DQt * Phis(:, :, i);
        end
        if i <= q
          Oomega((k*(r-1)+1):(k*r), (k*(c-1)+1):(k*c)) = DQt * Psis(:, :, i);
        end
      end
    end
  end

end



