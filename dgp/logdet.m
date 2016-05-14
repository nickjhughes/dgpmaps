
function v = logdet(A)
%LOGDET Calculate the log of the determine of the given matrix.
%
% v = logdet(A)
%
% If the matrix is positive definite, a Cholesky decomposition is used.

use_chol = true;
try
    L = chol(A);
catch err
    % Easiest way to check for positive-definiteness
    if (strcmp(err.identifier,'MATLAB:posdef'))
        use_chol = false;
    else
        rethrow(err);
    end
end

if use_chol
    v = 2*sum(log(diag(L)));
else
    [~, U, P] = lu(A);
    du = diag(U);
    c = det(P) * prod(sign(du));
    
    sign_u = sign(prod(du));
    sign_c = sign(c);
    
    v = log(abs(c)) + sum(log(abs(du)));
    v = v*sign_u*sign_c;
end
