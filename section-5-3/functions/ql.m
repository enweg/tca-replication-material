% Compute the QL decomposition of a matrix A.
%
% The QL decomposition factors a given matrix A into a product of an orthogonal matrix Q 
% and a lower triangular matrix L. This is similar to the QR decomposition but the triangular 
% matrix is lower triangular instead of upper triangular.
%
% Inputs:
%   A - Input matrix to be decomposed.
%
% Outputs:
%   Q - Orthogonal matrix.
%   L - Lower triangular matrix.
%
function [Q, L]=ql(A)
    k = size(A, 2);
    P = eye(k);
    P = P(k:-1:1, :);

    [Q, R] = qr(A*P);
    Q = Q*P;
    L = P*R*P;
end
