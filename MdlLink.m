function Gm = MdlLink(varargin)

    if iscell(varargin{1})
        arg = varargin{1};
    else
        arg = varargin;
    end
    
    Gm = arg{1};
    for n = 2:length(arg)
        Gm = append(Gm,arg{n});
    end

    seq = 0:(length(arg)-1);
    seq2 = zeros(1,2*length(arg));

    for n = 1:length(arg)
        seq2(2*n-1) = 3*seq(n)+2;
        seq2(2*n) = 3*seq(n)+3;
    end

    pvect = [3*seq+1,seq2];

    Gm.B = Gm.B(:,pvect);
    Gm.C = Gm.C(pvect,:);
    Gm.D = Gm.D(pvect,pvect);

end