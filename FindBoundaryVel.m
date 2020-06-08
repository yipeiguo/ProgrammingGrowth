% This is a function to solve for the normal velocities at the boundary of
% the tissue given the boundary points and the point sources of growth.
% This uses boundary integral formulation based on thesis by Ali Khalid.
% INPUTS:
% - Z: (N x 1 vector) complex coordinates (z) of points along tissue contour
% - gpts: (G x 1 vector) complex coordinates (z) of points within the growth zone
% - charges: (G x 1 vector) growth rate at those points
% - dS: length of boundary/N
% - gamma: smoothing parameter between 0 and infinity. The smaller gamma
% is, the higher the smoothing effect.
% OUTPUTS:
% - Vn: (N x 1 vector) normal velocity at each of the points specified in Z
% - Ux: (N x 1 vector) x component of the velocities in Vn (in cartesian
% coords)
% - Vy: (N x 1 vector) y component of the velocities in Vn (in cartesian
% coords)

function [Vn,Ux,Vy,thetavec] = ...
    FindBoundaryVel(charges,gpts,Z,dS,gamma)

% Extract variables
N = length(Z);
G = length(gpts);
assert(length(charges) == G, 'dimensions of charges must be same as that of gpts');

% coordinates of N boundary mid-points
Zmid = (-[Z(end);Z(1:end-1)] + 9.*Z + 9.*[Z(2:end);Z(1)] - [Z(3:end);Z(1:2)])./16;

% dZ/ds
dZds = ([Z(end-1:end);Z(1:end-2)] - 8.*[Z(end);Z(1:end-1)] + ...
    8.*[Z(2:end);Z(1)] - [Z(3:end);Z(1:2)])./(12*dS);

% dZmid/ds
dZmidds = (-[dZds(end);dZds(1:end-1)] + 9.*dZds + ...
    9.*[dZds(2:end);dZds(1)] - [dZds(3:end);dZds(1:2)])./16;

% Create matrix of coefficients M
M = (dS*1i)./(repmat(Z.',N,1)-repmat(Zmid,1,N)).*repmat(dZmidds,1,N);
diag = ([Z(2:end);Z(1)] - Zmid)./([Z(2:end);Z(1)] - Z);
offdiag = (Zmid-Z)./([Z(2:end);Z(1)] - Z);
M(1:N+1:end) = M(1:N+1:end) + diag.';
M(N+1:N+1:end) = M(N+1:N+1:end) + offdiag(1:end-1).';
M(N) = M(N) + offdiag(end);

% Calculate tangential velocities (from growth potential) 
zdiffMat = repmat(Z,1,G)-repmat(gpts.',N,1);
vtauhatMat = repmat(charges.',N,1)./(2*pi).*...
    real(repmat(dZds,1,G).*conj(zdiffMat))./abs(zdiffMat).^2;
vtauhat = sum(vtauhatMat,2);

zdiffMat_mid = repmat(Zmid,1,G)-repmat(gpts.',N,1);
vtauhatMat_mid = repmat(charges.',N,1)./(2*pi).*...
    real(repmat(dZmidds,1,G).*conj(zdiffMat_mid))./abs(zdiffMat_mid).^2;
vtauhat_mid = sum(vtauhatMat_mid,2);

% Create R vector
Rvec = -(pi*1i).*vtauhat_mid + dS.*dZmidds.*...
   sum(repmat(vtauhat.',N,1)./(repmat(Z.',N,1)-repmat(Zmid,1,N)),2);

% Define H matrix (for smoothing purposes)
% H is the third order difference operator
rowInd = repmat((1:N)',5,1);
colInd = [[N-1;N;(1:N-2)'];[N;(1:N-1)'];(1:N)';[(2:N)';1];[(3:N)';1;2]];
eleVec = [ones(N,1);-ones(N,1).*4;ones(N,1).*6;-ones(N,1).*4;ones(N,1)];
H = sparse(rowInd,colInd,eleVec,N,N);

% Create effective A matrix and B vector in order to solve Ax=B problem
Amat = M.'*M + H./gamma;
Bvec = M.'*Rvec;

% solve Ax = B where x = vn_tilde is the normal background velocity
vn_tilde = Amat\Bvec;
vn_tilde = real(vn_tilde);

% Calculate normal velocity that comes from growth potential
vnhatMat = repmat(charges.',N,1)./(2*pi).*...
    imag(repmat(dZds,1,G).*conj(zdiffMat))./abs(zdiffMat).^2;
vnhat = sum(vnhatMat,2);

% sum tilde and hat contributions
Vn = vn_tilde + vnhat;
Ux = imag(Vn.*dZds);
Vy = -real(Vn.*dZds);

div_vntilde = sum(vn_tilde)*dS;
if abs(div_vntilde) > 1e-3*sum(vnhat)
    warning('vntilde divergence exceeds 1e-3 of total divergence');
end

% thetavec: angle between tangent vector at that point and the x-axis
thetavec = atan2(imag(dZds),real(dZds));

end