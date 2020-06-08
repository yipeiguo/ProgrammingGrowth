% This is a function to get Vn due to growth given boundary points of the
% tissue. Allows optional respacing.

function [Vn,theta,Z,S,gpts,cmat,charges,ifintissue,X,Y] = ...
    GetVnFromZ_new(Z,S,secparams,decayparams,gparams,otherparams,gamma,prevData,respace)

% Extract variables
N = length(Z);
q = length(secparams.mumaxvec); % number of chemicals
xscan = otherparams.xscan; yscan = otherparams.yscan;
numx = length(xscan); numy = length(yscan);


% equally space out the points along boundary and recalculate total length
% (if desired)
if respace == true
    [Z,S] = respaceZ(Z);
end
dS = S/N;

% obtain concentration profiles, growth sources positions and strengths
if exist('prevData','var')
    cmat_guess = prevData.cmat;    
else
    cmat_guess = ones(numy,numx,q);
end

[cmat,gpts,charges,Vn,~,~,theta,ifintissue,X,Y,~] = ...
    FindSystemProp_new(Z,dS,gamma,secparams,gparams,decayparams,otherparams,cmat_guess);



end