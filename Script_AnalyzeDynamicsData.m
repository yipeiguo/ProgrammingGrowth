% This is a script for analyzing stored data of forward growth dynamics.

close all; 
clear variables;

%% load data workspace
% pathname = './DataFromCluster/1ginhib/';
% fn = '1ginhib_mueff_8';

pathname = '.\DataFromCluster\2ginhib\';
fn = '2ginhib_mu1eff_8_mu2eff_200_mu2maxregcoeff_0_k2overk1_1_secthreseff_0pt9_N_800_gamma_100';

% pathname = '.\DataFromCluster\1ginhib_1gthresreg\';
% fn = '1ginhib_1gthresreg_mu1eff_8_mu2eff_0pt1_mu2maxregcoeff_0_k2overk1_0pt2_secthresregsign_0_gthresregsign_-1';
load(strcat(pathname,fn));

%% Calculate how minimum concentrations at tissue tip and gz properties 
% vary as tissue grows
cminMat = zeros(numintervals,q);
gzSizeVec = zeros(numintervals);
for intIndxOI = 1:numintervals
    ZoI = Z_record(:,intIndxOI);
    gptsOI = gpts_record{intIndxOI};
    GoI = length(gptsOI);
    cmatOI = cmat_record{intIndxOI};  
    IndsInOI = Inds_intissue_record{intIndxOI};  
    X_intissue = X(IndsInOI);
    if GoI==0
        break
    end
    tipIndices = find(X_intissue(:)>max(real(ZoI))*0.9);
    for qindx = 1:q
        cmat_q = cmatOI(:,:,qindx);
        c_intissue = cmat_q(IndsInOI);
        cminMat(intIndxOI,qindx) = min(c_intissue(tipIndices));
    end
    gzSizeVec(intIndxOI) = GoI*secparams.darea;
end

if intIndxOI < numintervals
    finalint = intIndxOI - 1;
else
    finalint = numintervals;
end

%% Plot initial state
figure;
plot([yinit(1:N);yinit(1)],[yinit(1+N:end);yinit(1+N)],'k-','LineWidth',2);
hold on
if max(gvec_init) > 0
    scatter(real(gpts),imag(gpts),'ro','filled');
    hold on
%     quiver(real(Z),imag(Z),UV_init(:,1),UV_init(:,2),'Color','b');
else
    warning('no cells are growing!');
end
ylim([-0.6,0.6]);
% xlim([1.2,1.7]);
xlim([-0.1,1.7]);
axis equal 
    
%%
figure;
for qindx = 1:q
    subplot(q,1,qindx);
%     imagesc([xscan(xminInd),xscan(xmaxInd)],[yscan(yminInd),yscan(ymaxInd)],...
%         cmat_init_cut(:,:,qindx)); 
%     hold on
    plot(real([Z;Z(1)]),imag([Z;Z(1)]),'k-','markersize',3,'LineWidth',1.5);
    hold on
%     if max(gvec_init) > 0
%         scatter(real(gpts),imag(gpts),'ro','filled');
%         hold on
%         quiver(real(Z),imag(Z),UV_init(:,1),UV_init(:,2),'Color','g');
%     else
%         warning('no cells are growing!');
%     end
    ylim([-0.6,0.6]);
    axis equal 
    colorbar
end


%% Plot tissue dynamics
% close all
maxmksize = 36;
scalefactor = 1;
colormat = rand(numintervals,3);
colormat(:,1) = 1;
colormat(:,2:3) = 0;
figure;
% plot([yinit(1:N);yinit(1)],[yinit(1+N:end);yinit(1+N)],'k-','LineWidth',1.5);
% hold on  
% scatter(real(gpts),imag(gpts),'k','filled');
% hold on
% quiver(real(Z),imag(Z),UV_init(:,1).*scalefactor,UV_init(:,2).*scalefactor,...
%     'Color','k','AutoScale','on');
% hold on
% for intIndx = [(3:3:min(6,finalint)),min(finalint+1,50)]
for intIndx = 1
% for intIndx = [(5:5:min(50,finalint))]
% for intIndx = [5,10,25,finalint]
    Zvec = Z_record(:,intIndx);
    
    gptsVec = gpts_record{intIndx};
    chargesVec = charges_record{intIndx};
    thetaVec = theta_record(:,intIndx);
    VnVec = Vn_record(:,intIndx);
    VnxVec = VnVec.*cos(thetaVec-pi/2);
    VnyVec = VnVec.*sin(thetaVec-pi/2);
%     Uvec = U_record(:,intIndx);
%     UxVec = Uvec.*cos(thetaVec-pi/2);
%     UyVec = Uvec.*sin(thetaVec-pi/2);
    plot(real([Zvec;Zvec(1)]),imag([Zvec;Zvec(1)]),...
        '-','Color',colormat(intIndx,:),'LineWidth',1.5); % boundary points
    hold on
    scatter(real(gptsVec),imag(gptsVec),(chargesVec+1e-10)./max(chargesVec+1e-10).*maxmksize,...
        colormat(intIndx,:),'filled'); % growth charges
    hold on
    quiver(real(Zvec),imag(Zvec),VnxVec.*scalefactor,VnyVec.*scalefactor,...
        'AutoScale','on'); % velocity vector due to growth
    hold on
end
% quiver(real(Zvec),imag(Zvec),VnxVec.*scalefactor,VnyVec.*scalefactor,...
%     'Color','b','AutoScale','on'); % velocity vector due to growth    
% quiver(real(Zvec),imag(Zvec),UxVec.*scalefactor,UyVec.*scalefactor,...
%     'Color','r','AutoScale','on'); % velocity vector due to growth
hold on
plot(real([Zvec;Zvec(1)]),imag([Zvec;Zvec(1)]),...
    'k-','LineWidth',1.5); % boundary points
hold on
% scatter(real(gptsVec),imag(gptsVec),(chargesVec+1e-10)./max(chargesVec+1e-10).*maxmksize,...
%     'k','filled'); % growth charges
hold on
axis equal
xlim([0.5,2]);

%% plot how concentrations at tip and growth zone size change with time
figure;
for qindx = 1:min(q,3)
    subplot(2,2,qindx);
    plot(1:finalint,cminMat(1:finalint,qindx),'x');
    xlabel('int');
    ylabel('cmin');
    title(strcat('cmin (q=',num2str(qindx),')'));
end
subplot(2,2,4);
plot(1:finalint,gzSizeVec(1:finalint),'x');
xlabel('int');
ylabel('gz area');
title('gz area');

%% Here we visualize the different distinct regions (for each chemical, 
% we plot the region at which the chemical is secreted and the region at
% which conditions of that chemical allow for growth).
intIndxOI = max(min(20,finalint),1);
ZoI = Z_record(:,intIndxOI);
gptsOI = gpts_record{intIndxOI};
GoI = length(gptsOI);
cmatOI = cmat_record{intIndxOI};  
IndsInOI = Inds_intissue_record{intIndxOI};  
X_intissue = X(IndsInOI);
Y_intissue = Y(IndsInOI);
maxmksize = 36;
figure;
for qindx = 1:q
    ng = gparams.nvec(qindx);
    cmat_q = cmatOI(:,:,qindx);
    subplot(q,2,(qindx-1)*2+1);
    cthresmat = 1;
    if gparams.ifreggthres == true
        kvec = gparams.kMat(qindx,:);
        qoI = find(kvec~=0);
        if ~isempty(qoI)
            for kk = 1:length(qoI)
                cthresmat = cthresmat.*(1 + kvec(qoI(kk)).*cmatOI(:,:,qoI(kk))).^gparams.kpower;
%                 cthresmat = cthresmat + kvec(qoI(kk)).*cmat_curr(:,:,qoI(kk)).^gparams.kpower;
            end
            cthresmat = max(cthresmat,0);
        end
    end
    
    GrowthRelevancy = cmat_q.^ng./(cthresmat.^ng+cmat_q.^ng);
    GrowthRelevancy_intissue = GrowthRelevancy(IndsInOI);
    plot(real(ZoI),imag(ZoI),'kx-');
    hold on
    if max(GrowthRelevancy_intissue(:))>gparams.gzthres
        scatter(X_intissue(:),Y_intissue(:),...
            (GrowthRelevancy_intissue+1e-10)./max(1e-10+GrowthRelevancy_intissue).*maxmksize,...
            cmat_q(IndsInOI),'filled');
    end
    axis equal 
    xlabel('x'); ylabel('y');
    colorbar;
    title(strcat('growth region allowable by q = ',num2str(qindx)));
    
    mumax = secparams.mumaxvec(qindx);
    secRates = mumax.*ones(numy,numx);
    if secparams.ifmumaxreg == true
        muregvec = secparams.muregmat(qindx,:);
        qoI_mumaxreg = find(muregvec~=0);
        if ~isempty(qoI_mumaxreg)
            for kk = 1:length(qoI_mumaxreg)
                secRates = secRates.*(1 + muregvec(qoI_mumaxreg(kk)).*...
                    cmatOI(:,:,qoI_mumaxreg(kk))).^secparams.muregpower;
%                 secRates = secRates + muregvec(qoI_mumaxreg(kk)).*...
%                     cmat_curr(:,:,qoI_mumaxreg(kk)).^secparams.muregpower;
            end
            secRates = max(secRates,0);
        end
    end        
    
    if secparams.ifregsec == true
        nsec = secparams.nijmat(qindx,:);
        Ksec = secparams.Kijmat(qindx,:);
        for qindx2 = 1:q
            if nsec(qindx2)~=0
                secRates = secRates.*cmatOI(:,:,qindx2).^nsec(qindx2)./...
                    (Ksec(qindx2).^nsec(qindx2) + cmatOI(:,:,qindx2).^nsec(qindx2));
            end
        end        
    end
    secRates_intissue = secRates(IndsInOI);
        
    subplot(q,2,(qindx-1)*2+2);
    plot(real(ZoI),imag(ZoI),'kx-');
    hold on
    scatter(X_intissue(:),Y_intissue(:),...
        (secRates_intissue+1e-10)./max(secRates_intissue+1e-10).*maxmksize,...
        cmat_q(IndsInOI),'filled'); 
    axis equal 
    xlabel('x'); ylabel('y');
    colorbar;
    title(strcat('secretion region of q = ',num2str(qindx)));
end


%% here we plot spatial profile of threshold concentrations 
% (if there is spatial dependence)
nvec = gparams.nvec;
if gparams.ifreggthres == true
    figure;
    for qindx = 1:q
        if nvec(qindx)~=0
            subplot(q,2,(qindx-1)*2+1);
            plot(real(ZoI),imag(ZoI),'kx-');
            hold on

            gratesMat = gmax;
            kvec = gparams.kMat(qindx,:);
            cthresmat = 1;
            qoI = find(kvec~=0);
            if ~isempty(qoI)
                for kk = 1:length(qoI)
                    cthresmat = max(cthresmat + kvec(qoI(kk)).*cmatOI(:,:,qoI(kk)),0);
%                     cthresmat = max(cthresmat + kvec(qoI(kk)).*(cmatOI(:,:,qoI(kk))-1),0.1);
                end
            end
            gfactor = cmatOI(:,:,qindx).^nvec(qindx)./...
                (cthresmat.^nvec(qindx) + cmatOI(:,:,qindx).^nvec(qindx));
            gfactor(isnan(gfactor)) = 1;
            gratesMat = gratesMat.*gfactor;
            
            scatter(X_intissue(:),Y_intissue(:),[],...
                cthresmat(IndsInOI),'filled');
%             scatter(X_intissue(:),Y_intissue(:),...
%                 gratesMat_intissue./max(gratesMat_intissue).*maxmksize,...
%                 cthresmat_intissue,'filled');
            axis equal 
            xlabel('x'); ylabel('y');
            colorbar;        
            title(strcat('growth threshold for q = ',num2str(qindx)));

            subplot(q,2,qindx*2);
            plot(real(ZoI),imag(ZoI),'kx-');
            hold on
            scatter(X_intissue(:),Y_intissue(:),[],gfactor(IndsInOI),'filled');
            axis equal 
            xlabel('x'); ylabel('y');
            colorbar;   
            title(strcat('growth factor for q = ',num2str(qindx)));
    
        end
        
    end
end
















