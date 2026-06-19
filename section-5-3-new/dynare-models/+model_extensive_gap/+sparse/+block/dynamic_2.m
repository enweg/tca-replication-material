function [y, T, residual, g1] = dynamic_2(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
residual=NaN(8, 1);
  y(46)=y(35)*params(6)+y(38)*params(7);
  residual(1)=(y(35))-(y(68)-1/params(3)*y(49)-y(64)*1/params(3)*(params(14)-1));
  residual(2)=(y(55))-(y(46)-y(52));
  residual(3)=(0)-(y(61)*params(9)-y(65)*params(9)+y(55)*(1-params(9)));
  residual(4)=(y(42))-(params(2)*y(38)+(1-params(2))*y(9));
  residual(5)=(y(58))-(params(8)*(y(38)-y(5)-params(1)*(y(71)-y(38))));
  residual(6)=(params(4)*y(52))-(y(55)-params(3)*y(35));
  residual(7)=(y(58))-(params(1)*params(5)*y(94)+params(1)*(1-params(2))*y(91)-y(49));
  residual(8)=(y(61))-(y(46)-y(9));
if nargout > 3
    g1_v = NaN(27, 1);
g1_v(1)=(-(1-params(2)));
g1_v(2)=1;
g1_v(3)=params(8);
g1_v(4)=1/params(3);
g1_v(5)=1;
g1_v(6)=1;
g1_v(7)=params(4);
g1_v(8)=1;
g1_v(9)=(-(1-params(9)));
g1_v(10)=(-1);
g1_v(11)=1;
g1_v(12)=(-params(7));
g1_v(13)=(-params(2));
g1_v(14)=(-(params(8)*(1+params(1))));
g1_v(15)=(-params(7));
g1_v(16)=1;
g1_v(17)=(-params(6));
g1_v(18)=params(3);
g1_v(19)=(-params(6));
g1_v(20)=1;
g1_v(21)=1;
g1_v(22)=(-params(9));
g1_v(23)=1;
g1_v(24)=(-(params(8)*(-params(1))));
g1_v(25)=(-1);
g1_v(26)=(-(params(1)*(1-params(2))));
g1_v(27)=(-(params(1)*params(5)));
    if ~isoctave && matlab_ver_less_than('9.8')
        sparse_rowval = double(sparse_rowval);
        sparse_colval = double(sparse_colval);
    end
    g1 = sparse(sparse_rowval, sparse_colval, g1_v, 8, 24);
end
end
