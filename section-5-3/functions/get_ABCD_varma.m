% Get the ABCD representation used in Morris (2016).
%
% Inputs:
%   M_        - Returned by Dynare
%   options_  - Returned by Dynare
%   oo_       - Returned by Dynare
%   obs_var   - (Optional) Vector of observed variables. If not provided, 
%               the function will use the observed variables defined in the .mod file.
%
% Outputs:
%   Ax        - See Morris (2016)
%   Ay        - See Morris (2016)
%   B         - See Morris (2016)
%   Cx        - See Morris (2016)
%   Cy        - See Morris (2016)
%   D         - See Morris (2016)
%   vars      - Vector of variable IDs corresponding to the VARMA representation.
%
% References: 
%   - Morris, S. D. (2016). VARMA representation of DSGE models. 
%     Economics Letters, 138, 30–33. https://doi.org/10.1016/j.econlet.2015.11.027
function [Ax, Ay, B, Cx, Cy, D, vars]=get_ABCD_varma(M_, options_, oo_, obs_var)
  if ~isfield(options_,'varobs_id')
    warning('get_ABCD: No observables have been defined using the varobs-command.')
    return;    
  end

  %get state indices
  ipred = M_.nstatic+(1:M_.nspred)';
  %get observable position in decision rules
  if nargin == 3
    warning("Using observed variables from .mod file");
    obs_var=oo_.dr.inv_order_var(options_.varobs_id);
  end

  vars = [ipred; obs_var];

  %get state transition matrices
  [A,B] = kalman_transition_matrix(oo_.dr,ipred,1:M_.nspred);
  %get observation equation matrices
  [C,D] = kalman_transition_matrix(oo_.dr,obs_var,1:M_.nspred);

  % Stack them into a VAR
  ns = size(A, 1);
  nc = size(C, 1);
  A1 = [A zeros(ns, nc); C zeros(nc, nc)];
  Psi = [B; D];
  % Now sort them into observed and unobserved
  % that is, whenever a variable is observed and a state variable 
  % is will be twice in the current system. We need to remove the second 
  % appearance of it, and move the state equation into the observed block.
  for vo=obs_var'
    if any(ismember(vo, ipred))
      % observed variable is also state variable 
      ixs = find(vars==vo);
      ix_s = ixs(1);
      ix_o = ixs(2);

      % copy over row and column of A1
      A1(ix_o, :) = A1(ix_s, :);
      A1(:, ix_o) = A1(:, ix_s);;;;
      % then remove ix_s column and row 
      A1 = A1([1:(ix_s-1), (ix_s+1):size(A1, 1)], [1:(ix_s-1), (ix_s+1):size(A1, 1)]);

      % copy over row in Psi
      Psi(ix_o, :) = Psi(ix_s, :);
      % remove ix_s row
      Psi = Psi([1:(ix_s-1), (ix_s+1):size(Psi, 1)], :);

      % remove entry in vars
      vars = vars([1:(ix_s-1), (ix_s+1):size(vars, 1)]);
    end
  end


  n = size(A1, 1);
  no = size(obs_var, 1);
  ns = n - no;
  % removing all unnecessary state variables
  keep = repelem(true, size(A1, 1));
  for i=1:ns
    keep(i) = any(abs(A1([1:(i-1), (i+1):size(A1, 1)], i)) > 1e-10);
    if ~keep(i)
      ns = ns - 1;
      warning("Removing state variable " + get_variable_ordering(M_, options_, oo_, vars(i)));
    end
  end
  A1 = A1(keep, keep);
  Psi = Psi(keep, :);
  vars = vars(keep);

  D = Psi((ns+1):end, :);
  if size(D, 1) ~= size(D, 2)
    error("D is not square.");
  end
  Ax = A1(1:ns, 1:ns);
  Ay = A1(1:ns, (ns+1):end);
  Cx = A1((ns+1):end, 1:ns);
  Cy = A1((ns+1):end, (ns+1):end);
  B = Psi(1:ns, :);
end
