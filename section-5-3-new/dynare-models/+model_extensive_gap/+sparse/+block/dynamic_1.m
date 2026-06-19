function [y, T] = dynamic_1(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
  y(64)=params(14)*y(31)+x(1);
  y(65)=params(15)*y(32)+x(2);
  y(66)=params(16)*y(33)+x(3);
end
