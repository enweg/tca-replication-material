% Compute the Impulse Response Functions (IRFs) of a structural VARMA model.
%
% Inputs:
%   A0      - Initial impact matrix (n x n), where n is the number of variables.
%   Phis    - 3D array of autoregressive coefficients (n x n x p), where p is the AR order.
%   Psis    - 3D array of moving average coefficients (n x n x q), where q is the MA order.
%   horizon - Integer specifying the number of periods for which to compute the IRFs.
%
% Outputs:
%   irfs    - 3D array of IRFs (n x n x (horizon+1)). The dimensions correspond 
%             to the number of variables, the number of shocks, and the number 
%             of periods respectively.
%
function irfs=irf_varma(As, Psis, horizon)
  p = length(As) - 1;  % also includes A0;
  q = length(Psis);
  A0inv = inv(As{1});
  As = As(2:end);
  n = size(A0inv, 1);

  % calculating irfs
  irfs = zeros(n, n, horizon+1);
  irfs(:, :, 1) = A0inv;
  for h=1:horizon
    for i=1:min(p, h)
      irfs(:, :, h+1) = irfs(:, :, h+1) + As{i}*irfs(:, :, h-i+1);
    end
    if h <= q
      irfs(:, :, h+1) = irfs(:, :, h+1) + Psis{h};
    end
    irfs(:, :, h+1) = A0inv * irfs(:, :, h+1);
  end

end
