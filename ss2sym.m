% convert a state-space model to a symbolic transfer function
% this is very time consuming and should be used with care

function G = ss2sym(S)

    s = sym('s');
    I = eye(length(S.A));
    
    G = S.C*(s*I-S.A)^(-1)*S.B + S.D;
    
end