function [y, T, residual, g1] = static_6(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T)
residual=NaN(8, 1);
  residual(1)=(params(4)*y(19))-(y(22)-params(3)*y(2));
  residual(2)=(y(25))-(params(1)*params(5)*y(28)+params(1)*(1-params(2))*y(25)-y(16));
  residual(3)=(y(2))-(y(2)-1/params(3)*y(16)-y(31)*1/params(3)*(params(14)-1));
  residual(4)=(y(9))-(params(2)*y(5)+y(9)*(1-params(2)));
  residual(5)=(y(28))-(y(13)-y(9));
  residual(6)=(y(22))-(y(13)-y(19));
  residual(7)=(0)-(y(28)*params(9)-y(32)*params(9)+y(22)*(1-params(9)));
  residual(8)=(y(13))-(y(2)*params(6)+y(5)*params(7));
if nargout > 3
    g1_v = NaN(19, 1);
g1_v(1)=params(3);
g1_v(2)=(-params(6));
g1_v(3)=(-(params(1)*params(5)));
g1_v(4)=1;
g1_v(5)=(-params(9));
g1_v(6)=1;
g1_v(7)=1/params(3);
g1_v(8)=(-params(2));
g1_v(9)=(-params(7));
g1_v(10)=1-(1-params(2));
g1_v(11)=1;
g1_v(12)=params(4);
g1_v(13)=1;
g1_v(14)=(-1);
g1_v(15)=1;
g1_v(16)=(-(1-params(9)));
g1_v(17)=(-1);
g1_v(18)=(-1);
g1_v(19)=1;
    if ~isoctave && matlab_ver_less_than('9.8')
        sparse_rowval = double(sparse_rowval);
        sparse_colval = double(sparse_colval);
    end
    g1 = sparse(sparse_rowval, sparse_colval, g1_v, 8, 8);
end
end
