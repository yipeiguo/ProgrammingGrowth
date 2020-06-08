% This is a function to get steady-state concentrations on a grid, 
% growth rates and positions of growth sources, and velocties at boundary of
% tissue, given a set of rules and parameters

% 12/22/19. After writing 'Findcss_new.m'/

function [cmat_curr,gpts,charges,Vn,UV,VnRatio,theta,ifintissue,X,Y,bpts] = ...
    FindSystemProp_new(Z,dS,gamma,secparams,gparams,decayparams,otherparams,cmat_guess)

q = length(secparams.mumaxvec); % number of chemicals

% extract variables from secparams and gparams. 
nvec = gparams.nvec;
gmax = gparams.gmax;
if isfield(gparams,'ifgthres')
    ifgthres = gparams.ifgthres;
else
    ifgthres = false;
end
if isfield(gparams,'gzOnlyTip')
    gzOnlyTip = gparams.gzOnlyTip;
else
    gzOnlyTip = false;
end
% matrix representing how threshold concentrations for growth changes
% linearly with concentrations of chemicals
if gparams.ifreggthres == true
    kMat = gparams.kMat;
    kpower = gparams.kpower;
end
% if isfield(gparams,'kMat') 
%     kMat = gparams.kMat;
% else
%     kMat = zeros(q,q);
% end

% remove gpts that are too close to tissue boundary
if isfield(otherparams,'limit')
    limit = otherparams.limit;
else
    limit = 0;
end

% find steady-state concentrations
[cmat_curr,ifintissue,success,~,X,Y,bpts] = ...
    Findcss_new(secparams,decayparams,otherparams,[real(Z),imag(Z)],dS,cmat_guess);
assert(success == true, 'css not found!');

% calculate growth rates based on these concentrations (i.e. 'growth charges')
gratesMat = gmax;
for qindx = 1:q
    if nvec(qindx)~=0
        cthresmat = 1;
        if gparams.ifreggthres == true
            kvec = kMat(qindx,:);        
            qoI = find(kvec~=0);
            if ~isempty(qoI)
                for kk = 1:length(qoI)
                    cthresmat = cthresmat.*(1 + kvec(qoI(kk)).*cmat_curr(:,:,qoI(kk))).^kpower;
%                     cthresmat = cthresmat + kvec(qoI(kk)).*cmat_curr(:,:,qoI(kk)).^kpower;
    %                 cthresmat = cthresmat + kvec(qoI(kk)).*max(cmat_curr(:,:,qoI(kk))-1,0);
    %                 cthresmat = max(cthresmat + kvec(qoI(kk)).*(cmat_curr(:,:,qoI(kk))-1),0.1);
                end
                cthresmat = max(cthresmat,0);
            end
        end
        gfactor = cmat_curr(:,:,qindx).^nvec(qindx)./...
            (cthresmat.^nvec(qindx) + cmat_curr(:,:,qindx).^nvec(qindx));
        gfactor(isnan(gfactor)) = 1;
        gratesMat = gratesMat.*gfactor;
    end
end
gratesMat(~ifintissue)=0;
gzInds = find(gratesMat>=gmax*gparams.gzthres);
if gzOnlyTip == true
    gzInds = gzInds(X(gzInds) > max(X(ifintissue))/2);
end
% romove points that are too close to boundary of tissue:
if limit > 0
    gzpts = [X(gzInds),Y(gzInds)];
    pwdist = pdist2(gzpts,bpts);
    mindist2bound = min(pwdist,[],2);
    gzInds = gzInds(mindist2bound>limit);    
end
gpts = X(gzInds)+Y(gzInds).*1i;
if ifgthres == true
    gratesMat = gmax.*(gratesMat>=gmax*gparams.gzthres);
end
charges = gratesMat(gzInds).*secparams.darea;

% find boundary velocities given growth charges
[Vn,Ux,Vy,theta] = FindBoundaryVel(charges,gpts,Z,dS,gamma);
UV = [Ux,Vy];
VnRatio = Vn(2:end)./Vn(1);

end

