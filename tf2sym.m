% convert a transfer function model to symbolic
% this is much faster than ss2sym.m

function G = tf2sym(S)

    s = sym('s');
    
    [Num,Den] = tfdata(S);
    G = poly2sym(cell2mat(Num),s)/poly2sym(cell2mat(Den),s);
    
end