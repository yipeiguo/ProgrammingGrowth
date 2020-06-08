% This is a function to respace Z

function [Znew,Snew] = respaceZ(Z)

    N = length(Z);
    
    px = real([Z;Z(1)]);
    py = imag([Z;Z(1)]);
    bpts = interparc(N+1,px,py,'spline');
    Znew = bpts(1:N,1) + bpts(1:N,2).*1i;
    Snew = arclength(px,py,'s');    

end