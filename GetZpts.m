% This is a function to obtain initial boundary points for various standard
% shapes

function [Z,dS] = GetZpts(shapeprop,N)

% Extract parameteres
shape = shapeprop.shape;
h = shapeprop.h;
len = shapeprop.len;
    
% options_fsolve = optimset('Display','on');
options_fsolve = optimset('Display','off');

if strcmp(shape,'parabola')
    
    a = h/sqrt(len);
    f = a^2/4; % focal length of parabola
    qfunc = @(w) sqrt(f^2 + w.^2);
    sfunc = @(w) w.*qfunc(w)./f + f.*log((w+qfunc(w))./f);
    Lhalf = sfunc(h/2) + h;
    dS = 2*Lhalf/N; 
    m = (N-2)/2-floor(h/dS);
    svec = (1:1:m).'.*dS;
    wvec = zeros(m,1);
    wguess = 0;
    for sindx = 1:length(svec)
        wvec(sindx) = fsolve(@(x) sfunc(x)-svec(sindx),wguess,options_fsolve);
        wguess = wvec(sindx);
    end
    yvec_arc = wvec.*2;
    xvec_arc = len.*(1-yvec_arc.^2./(h^2));
    xvec = [len;xvec_arc;zeros(N-m*2-1,1);flipud(xvec_arc)];
    yvec = [0;yvec_arc;dS.*(floor(h/dS):-1:1).';0;-(1:1:floor(h/dS)).'.*dS;...
        -flipud(yvec_arc)];
    Z = xvec + yvec.*1i;
    
elseif strcmp(shape,'ellipse')
    assert(len>=h,'len should be the semi-major axis and h the semi-minor axis!');

    mparam = 1-(len/h)^2;
    sfunc = @(x) h.*(ellipticE(acos(x/len),mparam)-ellipticE(0,mparam));
%     eccentricity = sqrt(1-(h/len)^2);
%     [~,E] = ellipke(eccentricity^2);
%     C = 2*len*E;
%     L = C+2*h;
    L = (sfunc(0)+h)*2; % total length
    dS = L/N;
    
    m = (N-2)/2-floor(h/dS);
    svec = (1:1:m).'.*dS;
    xvec_arc = zeros(m,1);
    xguess = len;
    for sindx = 1:length(svec)
        xvec_arc(sindx) = fsolve(@(x) sfunc(x)-svec(sindx),xguess,options_fsolve);
        xguess = xvec_arc(sindx);
    end
    yvec_arc = h.*sqrt(1-(xvec_arc./len).^2);
    xvec = [len;xvec_arc;zeros(N-m*2-1,1);flipud(xvec_arc)];
    yvec = [0;yvec_arc;dS.*(floor(h/dS):-1:1).';0;-(1:1:floor(h/dS)).'.*dS;...
        -flipud(yvec_arc)];
    Z = real(xvec) + real(yvec).*1i;
    
end


end