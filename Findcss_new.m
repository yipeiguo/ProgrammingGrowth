% This is a function to solve for steady-state chemical concentrations on a
% grid given the secretion rules.

function [cmat_curr,ifintissue,success,errorvec,X,Y,bpts] = ...
    Findcss_new(secparams,decayparams,otherparams,bpts,dS,cmat_guess)

% Extract variables/parameters
maxiter = otherparams.maxiter;
xscan = otherparams.xscan; yscan = otherparams.yscan;
dx = xscan(2)-xscan(1); dy = yscan(2)-yscan(1); 
delta = min(dx,dy);
numx = length(xscan); numy = length(yscan);

mumaxvec = secparams.mumaxvec;
q = length(mumaxvec);

ifregsec = secparams.ifregsec;
if ifregsec == true
    nijmat = secparams.nijmat;
    Kijmat = secparams.Kijmat;
    [regulated,regulator] = find(nijmat~=0);
    
    if secparams.ifregsecthres == true
        Qmat = secparams.Qmat;
        Qpower = secparams.Qpower;
    end    
end

ifmumaxreg = secparams.ifmumaxreg;
if ifmumaxreg == true
    muregmat = secparams.muregmat;
    muregpower = secparams.muregpower;
end

decayRates = decayparams.decayRates;

% redefine boundary points if necessary (original points too far spaced)
if dS > delta
    S = dS*size(bpts,1);
    N = ceil(S/delta);
    bpts = interparc(N+1,[bpts(:,1);bpts(1,1)],[bpts(:,2);bpts(1,2)],'spline');
    bpts = bpts(1:N,:);
end

% select points that are inside tissue:
[X,Y] = meshgrid(xscan,yscan);
ifintissue = inpolygon(X(:),Y(:),bpts(:,1),bpts(:,2));

% Extract points that are along domain boundary:
Inds_lb = find(X == min(xscan));
Inds_rb = find(X == max(xscan));
Inds_bb = find(Y == min(yscan)); % bottom boundary (smallest y)
Inds_tb = find(Y == max(yscan)); % top boundary
Inds_allb = union(union(Inds_lb,Inds_rb),union(Inds_bb,Inds_tb));
Inds_notb = setdiff((1:numx*numy).',Inds_allb);
numInds_notb = length(Inds_notb);
numInds_b = length(Inds_allb);

% construct matrix of coefficients to solve the linear diffusion eqn
rowinds = repmat(Inds_notb,5,1);
colinds = [Inds_notb;Inds_notb-numy;Inds_notb+numy;Inds_notb-1;Inds_notb+1];

rowinds_BC = Inds_allb;
colinds_BC = Inds_allb;
elems_BC = ones(numInds_b,1);


success = false;
cmat_curr = cmat_guess;
errorvec = zeros(1,maxiter);
for iterindx = 1:maxiter
    
    % Calculate secretion and degradation rates based on concentrations
    mumat = zeros(numy,numx,q);  
    decayRateMat = zeros(numy,numx,q);  
    for qindx = 1:q
        mumat_q = mumaxvec(qindx).*ones(numy,numx);
        if ifmumaxreg == true
            muregvec = muregmat(qindx,:);
            qoI_mumaxreg = find(muregvec~=0);
            if ~isempty(qoI_mumaxreg)
                for kk = 1:length(qoI_mumaxreg)
                    mumat_q = mumat_q.*(1 + muregvec(qoI_mumaxreg(kk)).*cmat_curr(:,:,qoI_mumaxreg(kk))).^muregpower;
%                     mumat_q = mumat_q + muregvec(qoI_mumaxreg(kk)).*cmat_curr(:,:,qoI_mumaxreg(kk)).^muregpower;
                end
                mumat_q = max(mumat_q,0);
            end
        end
            
        if ifregsec == true
            
            nijvec = nijmat(qindx,:);
            Kijvec = Kijmat(qindx,:);
            for qindx2 = 1:q
                if nijvec(qindx2)~=0
                    cthresmat = Kijvec(qindx2);
                    if secparams.ifregsecthres == true
                        [~,whichind] = ismember([qindx,qindx2],[regulated,regulator],'rows');
                        Qvec = Qmat(whichind,:);
                        qoI = find(Qvec~=0);
                        if ~isempty(qoI)
                            for kk = 1:length(qoI)
                                cthresmat = cthresmat.*(1 + Qvec(qoI(kk)).*cmat_curr(:,:,qoI(kk))).^Qpower;
%                                 cthresmat = cthresmat + Qvec(qoI(kk)).*cmat_curr(:,:,qoI(kk)).^Qpower;
    %                             cthresmat = max(cthresmat + Qvec(qoI(kk)).*(cmat_curr(:,:,qoI(kk))-1),1e-3);
                            end
                            cthresmat = max(cthresmat,0);
                        end
                    end
                    factor = cmat_curr(:,:,qindx2).^nijvec(qindx2)./...
                        (cthresmat.^nijvec(qindx2) + ...
                        cmat_curr(:,:,qindx2).^nijvec(qindx2));
                    factor(isnan(factor)) = 1;
                    mumat_q = mumat_q.*factor;
                end
            end
        end
        mumat_q(~ifintissue) = 0;
        mumat(:,:,qindx) = reshape(mumat_q,[numy,numx]);
        decayRateMat(:,:,qindx) = ones(numy,numx).*decayRates(qindx);
    end    
    
    % Calculate concentrations based on these secretion rates
    cmat_next = zeros(numy,numx,q);
    for qindx = 1:q
        decayRateMat_q = decayRateMat(:,:,qindx);
        mumat_q = mumat(:,:,qindx);
        % update coefficient matrix
        elems = [-ones(numInds_notb,1).*(2/dx^2+2/dy^2)-decayRateMat_q(Inds_notb);...
            ones(numInds_notb*2,1)./(dx^2); ones(numInds_notb*2,1)./(dy^2)];
        coeffmat = sparse(rowinds,colinds,elems,numy*numx,numy*numx) + ...
            sparse(rowinds_BC,colinds_BC,elems_BC,numy*numx,numy*numx);
        
        % update RHS vector (source term)
        RHSvec = -mumat_q(:);

        % Solve for concentrations
        cvec_q = coeffmat\RHSvec;

        % convert back to full matrix
        cmat_q = reshape(cvec_q,[numy,numx]);
        cmat_next(:,:,qindx) = cmat_q;
    end
    
    % compare with original cmat
    cdiffmax = max(abs(cmat_next-cmat_curr),[],'all');
    errorvec(iterindx) = cdiffmax;
    
    if cdiffmax <= otherparams.ceps
        success = true;
        break
    else
        cmat_curr = cmat_next;
    end
end
errorvec = errorvec(1:iterindx);
assert(success==true,'concentrations not found!');


end

