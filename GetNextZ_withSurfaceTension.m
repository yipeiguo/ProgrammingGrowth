% This is a function that outputs the coordinates of the boundary points at 
% the next time step given Vn due to growth (calculated separately).
% Takes into account surface tension, and contribution to velocity due to
% surface tension is calculated implicitly.

% INPUTS:
% - y: (2N x 1) vector consisting of the x followed by y coordinates of
% each of the boundary points
% - latticeparams: structure containing parameters of the lattice grid that
% will be used to specify the positions of the point growth sources
% - secparams: structure containing the secretion parameters (used to
% calculate steady-state concentrations which determines growth rates)
% - gparams: structure containing growth parameters

function [Znext,Lnext,theta_next,kappa_next,U] = ...
    GetNextZ_withSurfaceTension(Zcurr,theta_curr,Vn,L_curr,zeta,dt,includeT)

% Extract parameters
N = length(Vn);
alphavec = (1:N-1).'; % this should be column vector
alphamax = N;
dS = L_curr/N;

if zeta > 0 || includeT == true
    % define alternative angle vector that is between 0 and 2pi
    theta_alt = theta_curr;
    theta_alt(theta_alt<0) = 2*pi + theta_alt(theta_alt<0); 

    % Calculate velocities due to surface tension
    dthetadalpha = ([theta_curr(end-1:end);theta_curr(1:end-2)] - 8.*[theta_curr(end);theta_curr(1:end-1)] + ...
        8.*[theta_curr(2:end);theta_curr(1)] - [theta_curr(3:end);theta_curr(1:2)])./12;
    dthetadalpha_alt = ([theta_alt(end-1:end);theta_alt(1:end-2)] - 8.*[theta_alt(end);theta_alt(1:end-1)] + ...
        8.*[theta_alt(2:end);theta_alt(1)] - [theta_alt(3:end);theta_alt(1:2)])./12;
    dthetadalpha(abs(dthetadalpha)>0.75*pi) = dthetadalpha_alt(abs(dthetadalpha)>0.75*pi);
    dthetadalpha_mid = (-[dthetadalpha(end);dthetadalpha(1:end-1)] + 9.*dthetadalpha + ...
        9.*[dthetadalpha(2:end);dthetadalpha(1)] - [dthetadalpha(3:end);dthetadalpha(1:2)])./16;
    kappa_mid = dthetadalpha_mid./dS;
    Utilde_mid = -zeta.*kappa_mid;

    % total normal velocity:
    Vn_mid = (-[Vn(end);Vn(1:end-1)] + 9.*Vn + 9.*[Vn(2:end);Vn(1)] - [Vn(3:end);Vn(1:2)])./16;
    U_mid = Vn_mid + Utilde_mid;

    % calculate tangential velocities at each point (set such that spacing
    % between points stay the same during the evolution)
    thetaprimeU_mid = -dthetadalpha_mid.*U_mid;
    thetaprimeU_cumsum = cumsum(thetaprimeU_mid);
    Tvec = [0; thetaprimeU_cumsum(1:end-1)- alphavec./N.*thetaprimeU_cumsum(end)];
    
    % calculate length of boundary at the next time step
    Lnext = L_curr - dt*thetaprimeU_cumsum(end);
    dS_next = Lnext/N;
    
end

if zeta > 0
    
    % solve for theta at next time step
    % Define second order central difference matrix:
    rowInd = repmat((1:N)',3,1);
    colInd = [[N;(1:N-1)'];(1:N)';[(2:N)';1]];
    eleVec = [ones(N,1);-ones(N,1).*2;ones(N,1)];
    M = speye(N) - (zeta*alphamax/(dS_next*Lnext)).*sparse(rowInd,colInd,eleVec,N,N);
    % RHS vector:
    dVndalpha = ([Vn(end-1:end);Vn(1:end-2)] - 8.*[Vn(end);Vn(1:end-1)] + ...
        8.*[Vn(2:end);Vn(1)] - [Vn(3:end);Vn(1:2)])./12;
    RHSvec = theta_curr + (dt*alphamax/L_curr).*(dVndalpha + Tvec.*dthetadalpha);

    theta_next = M\RHSvec;
    theta_next_alt = theta_next;
    theta_next_alt(theta_next_alt<0) = 2*pi + theta_next_alt(theta_next_alt<0); 

    % Actual x and y velocities of every point (with Utilde calculated at new
    % time step)
    dthetadalpha_next = ([theta_next(end-1:end);theta_next(1:end-2)] - 8.*[theta_next(end);theta_next(1:end-1)] + ...
        8.*[theta_next(2:end);theta_next(1)] - [theta_next(3:end);theta_next(1:2)])./12;
    dthetadalpha_next_alt = ([theta_next_alt(end-1:end);theta_next_alt(1:end-2)] - 8.*[theta_next_alt(end);theta_next_alt(1:end-1)] + ...
        8.*[theta_next_alt(2:end);theta_next_alt(1)] - [theta_next_alt(3:end);theta_next_alt(1:2)])./12;
    dthetadalpha_next(abs(dthetadalpha_next)>0.75*pi) = dthetadalpha_next_alt(abs(dthetadalpha_next)>0.75*pi);
    kappa_next = dthetadalpha_next./dS_next;
    Utilde_next = -zeta.*kappa_next;

    U = Vn + Utilde_next;
    Ux = Vn.*cos(theta_curr-pi/2) + Utilde_next.*cos(theta_next-pi/2);
    Vy = Vn.*sin(theta_curr-pi/2) + Utilde_next.*sin(theta_next-pi/2);
else
    U = Vn;
    Ux = Vn.*cos(theta_curr-pi/2);
    Vy = Vn.*sin(theta_curr-pi/2);
end

if includeT == true
    Uxtot = Tvec.*cos(theta_curr) + Ux;
    Vytot = Tvec.*sin(theta_curr) + Vy;
%     Uxtot = Tvec.*cos(theta_curr) + U.*cos(theta_curr-pi/2);
%     Vytot = Tvec.*sin(theta_curr) + U.*sin(theta_curr-pi/2);
elseif includeT == false
    Uxtot = Ux;
    Vytot = Vy;
%     Uxtot = U.*cos(theta_curr-pi/2);
%     Vytot = U.*sin(theta_curr-pi/2);
end

% get Z coordinates at next time point
Znext = Zcurr + (Uxtot + 1i.*Vytot).*dt;

end