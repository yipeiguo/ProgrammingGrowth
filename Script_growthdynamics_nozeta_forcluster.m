% This is a script to carry out growth dynamics using BIM, including a
% surface tension term to smooth out the surface (i.e. hopefully reduce
% surface irregularities) and solving the evolution equation in theta-L
% space.

tic

%% Define parameters
% tissue shape parameters
len = 1.5;
h = 0.5;
shape = 'parabola'; % either 'parabola' or 'ellipse'

% lattice grid
delta = 0.01; 
xmin = -3; xmax = len+4; hmin = -h-3; hmax = h+3;
xscan = xmin+delta/2:delta:xmax;
yscan = hmin+delta/2:delta:hmax-delta/2;
numx = length(xscan); numy = length(yscan);
otherparams.xscan = xscan; otherparams.yscan = yscan;

% Discretize boundary and boundary integral method parameters
N = 800; 
shapeprop.shape = shape;
shapeprop.h = h;
shapeprop.len = len;
[Z,dS] = GetZpts(shapeprop,N);
yinit = [real(Z);imag(Z)];

% secretion and growth rules
% if ruletype is 'unspecified', define parameters in script, otherwise 
% define effective parameters
ruletype = '2ginhib'; 
% ruletype = '1ginhib_1gthresreg'; 
if strcmp(ruletype,'1ginhib')
    effparams.mueff = 15;    
elseif strcmp(ruletype,'1ginhib_1gthresreg')
    effparams.mu1eff = 8;
    effparams.mu2eff = 5;
    effparams.mu2max_regcoeff = 100;
    effparams.k2overk1 = 0.2;
    effparams.secthres_regsign = 1; % must be either 1 or -1
    if effparams.secthres_regsign ~= 0
        effparams.secthres_eff = 1.6;
    end
    effparams.gthres_regsign = -1; % must be either 1 or -1
elseif strcmp(ruletype,'2ginhib')
    effparams.mu1eff = 8;
    effparams.mu2eff = 20;
    effparams.mu2max_regcoeff = 0;
    effparams.k2overk1 = 1;
    effparams.secthres_eff = 1.2;
end

[secparams,gparams,decayparams] = SpecifyRuleParams(ruletype,effparams);
q = length(secparams.mumaxvec);

% common parameters
gmax = 1; gparams.gmax = gmax;    
gparams.ifgthres = false;
gparams.gzthres = gmax*0.5;
gparams.gzOnlyTip = true;
secparams.darea = delta^2;

% simulation method and parameters
errorthres = 1e-3;
xstep_max = 0.01;
numintervals = 50;
numtsteps = 5; % number of iterations in each interval
respace = false; % whether to respace boundary points at every iteraction

% smoothing parameters
gamma = 100; % for smoothness of Vn

% other technical parameters (for obtaining ss concs)
otherparams.maxiter = 100;
otherparams.ceps = 1e-3;
otherparams.limit = 0.01;

% Specify initial guess for concentrations
cmat_guess = ones(numy,numx,q);


toc
disp('hyperparameters specified');

%% Specify file name for saving data
if strcmp(ruletype,'1ginhib')
    filename = strcat(ruletype,'_mueff_',num2str(effparams.mueff));
elseif strcmp(ruletype,'1ginhib_1gthresreg')
    if effparams.secthres_regsign == 0
        filename = strcat(ruletype,'_mu1eff_',num2str(effparams.mu1eff),...
            '_mu2eff_',num2str(effparams.mu2eff),...
            '_mu2maxregcoeff_',num2str(effparams.mu2max_regcoeff),...
            '_k2overk1_',num2str(effparams.k2overk1),...
            '_secthresregsign_',num2str(effparams.secthres_regsign),...
            '_gthresregsign_',num2str(effparams.gthres_regsign));
    else
        filename = strcat(ruletype,'_mu1eff_',num2str(effparams.mu1eff),...
            '_mu2eff_',num2str(effparams.mu2eff),...
            '_mu2maxregcoeff_',num2str(effparams.mu2max_regcoeff),...
            '_k2overk1_',num2str(effparams.k2overk1),...
            '_secthresregsign_',num2str(effparams.secthres_regsign),...
            '_secthreseff_',num2str(effparams.secthres_eff),...
            '_gthresregsign_',num2str(effparams.gthres_regsign));
    end    
elseif strcmp(ruletype,'2ginhib')
    filename = strcat(ruletype,'_mu1eff_',num2str(effparams.mu1eff),...
            '_mu2eff_',num2str(effparams.mu2eff),...
            '_mu2maxregcoeff_',num2str(effparams.mu2max_regcoeff),...
            '_k2overk1_',num2str(effparams.k2overk1),...
            '_secthreseff_',num2str(effparams.secthres_eff));
end

filename = strrep(filename,'.','pt');
filename = strcat(filename,'_N_',num2str(N),'_gamma_',num2str(gamma));

%% Calculate properties of system at initial time
tic
[cmat_curr,gpts,charges,Vn,UV,VnRatio,theta,ifintissue,X,Y,bpts] = ...
    FindSystemProp_new(Z,dS,gamma,secparams,gparams,decayparams,otherparams,cmat_guess);
toc

gvec_init = charges;
UV_init = UV;
Vn_init = Vn;
theta_init = theta;
cmat_init = cmat_curr;
ifin_init = ifintissue;

% Cut out cmat that is within tissue
xminInd = find(xscan<0,1,'last');
xmaxInd = find(xscan>len(1),1,'first');
yminInd = find(yscan<-h(1),1,'last');
ymaxInd = find(yscan>h(1),1,'first');
cmat_init_cut = cmat_init(yminInd:ymaxInd,xminInd:xmaxInd,:);
X_cut = X(yminInd:ymaxInd,xminInd:xmaxInd);
Y_cut = Y(yminInd:ymaxInd,xminInd:xmaxInd);


%% Storage arrays
Z_record = zeros(N,numintervals);
Z_record2 = zeros(N,numintervals);
Vn_record = zeros(N,numintervals);
% U_record = zeros(N,numintervals);
theta_record = zeros(N,numintervals);
L_record = zeros(1,numintervals);
L_record2 = zeros(1,numintervals);
% kappa_record = zeros(N,numintervals);
gpts_record = cell(1,numintervals);
cmat_record = cell(1,numintervals);
trecord = zeros(1,numintervals);
dtmat = zeros(numtsteps,numintervals);
charges_record = cell(1,numintervals);

Inds_intissue_record = cell(1,numintervals);
    
%% Main simulation
tic    
% ifVnexist = false;
% ifVnhalfexist = false;
t = 0;
% initialize variables to be used in main simulation
Zcurr = Z;
Lcurr = dS*N;
prevData.cmat = cmat_init;
% Calculate Vn at initial point
[Vncurr,thetacurr2,Zcurr2,Lcurr2,gptscurr,cmatcurr,chargescurr,ifintissuecurr] = ...
    GetVnFromZ_new(Zcurr,Lcurr,secparams,decayparams,gparams,otherparams,gamma,prevData,respace);
for intIndx = 1:numintervals
    fprintf('intIndx: %d \n',intIndx);        
    for tindx = 1:numtsteps
        fprintf('t: %d \n',tindx);
        
        if max(Vncurr) < 1e-5
            break
        end
        
        prevData.cmat = cmatcurr;
        
        dt = xstep_max/max(Vncurr);
        % calculate Z at dt later (1 step)
        Uxcurr = Vncurr.*cos(thetacurr2-pi/2);
        Vycurr = Vncurr.*sin(thetacurr2-pi/2);
        Znext = Zcurr2 + (Uxcurr + 1i.*Vycurr).*dt;
        [Znext,Lnext] = respaceZ(Znext);
        % calculate Z at dt/2 later
        Zhalf = Zcurr2 + (Uxcurr + 1i.*Vycurr).*(dt/2);
        [Zhalf,Lhalf] = respaceZ(Zhalf);
        
        % calculate Vn at half point
        [Vn_half,theta_half2,Z_half2,L_half2,gpts_half,cmat_half,charges_half,ifintissue_half] = ...
            GetVnFromZ_new(Zhalf,Lhalf,secparams,decayparams,gparams,otherparams,gamma,prevData,respace);

        prevData.cmat = cmat_half;    
        Ux_half = Vn_half.*cos(theta_half2-pi/2);
        Vy_half = Vn_half.*sin(theta_half2-pi/2);
        Znext_2s = Z_half2 + (Ux_half + 1i.*Vy_half).*(dt/2);
        [Znext_2s,Lnext_2s] = respaceZ(Znext_2s);
        
        % compare doing 1 step vs doing 2 steps 
        maxZdiff = max(abs(Znext-Znext_2s));
        if maxZdiff < errorthres
            Zcurr = Znext; Lcurr = Lnext;
            
            % calculate Vn at current point
            [Vncurr,thetacurr2,Zcurr2,Lcurr2,gptscurr,cmatcurr,chargescurr,ifintissuecurr] = ...
                GetVnFromZ_new(Zcurr,Lcurr,secparams,decayparams,gparams,otherparams,gamma,prevData,respace);            
        else
            while maxZdiff > errorthres
                % update dt and other variables
                dt = dt/2;
                Znext = Zhalf; Lnext = Lhalf; 
%                 theta_next = theta_half; kappa_next = kappa_half; UtoNext = UtoHalf;
                Vn_next = Vn_half; theta_next2 = theta_half2; Znext2 = Z_half2; 
                Lnext2 = L_half2; gpts_next = gpts_half; 
                cmat_next = cmat_half; charges_next = charges_half; 
                ifintissue_next = ifintissue_half;
                
                % calculate Z at new dt/2 later
                Uxcurr = Vncurr.*cos(thetacurr2-pi/2);
                Vycurr = Vncurr.*sin(thetacurr2-pi/2);        
                Zhalf = Zcurr2 + (Uxcurr + 1i.*Vycurr).*(dt/2);
                [Zhalf,Lhalf] = respaceZ(Zhalf);
        
                % recalculate Vn at new dt/2 later
                [Vn_half,theta_half2,Z_half2,L_half2,gpts_half,cmat_half,charges_half,ifintissue_half] = ...
                    GetVnFromZ_new(Zhalf,Lhalf,secparams,decayparams,gparams,otherparams,gamma,prevData,respace);
                prevData.cmat = cmat_half;        
        
                % calculate new 2-step Z:
                Ux_half = Vn_half.*cos(theta_half2-pi/2);
                Vy_half = Vn_half.*sin(theta_half2-pi/2);
                Znext_2s = Z_half2 + (Ux_half + 1i.*Vy_half).*(dt/2);
                [Znext_2s,Lnext_2s] = respaceZ(Znext_2s);
        
                % compare doing 1 step vs doing 2 steps
                maxZdiff = max(abs(Znext-Znext_2s));        
            end
            
            % update states
            Zcurr = Znext; Lcurr = Lnext; 
            
            Vncurr = Vn_next; thetacurr2 = theta_next2; Zcurr2 = Znext2; 
            Lcurr2 = Lnext2; gptscurr = gpts_next; cmatcurr = cmat_next;
            chargescurr = charges_next; ifintissuecurr = ifintissue_next;
%             ifVnexist = true;
        end 
        dtmat(intIndx,tindx) = dt;
        t = t + dt;            
    end
    Z_record(:,intIndx) = Zcurr;
    L_record(intIndx) = Lcurr;
%     kappa_record(:,intIndx) = kappa_next;
%     U_record(:,intIndx) = UtoNext;
    
    Z_record2(:,intIndx) = Zcurr2;
    L_record2(intIndx) = Lcurr2;
    Vn_record(:,intIndx) = Vncurr;
    theta_record(:,intIndx) = thetacurr2;    
    gpts_record{intIndx} = gptscurr;
    cmat_record{intIndx} = cmatcurr;
    charges_record{intIndx} = chargescurr;  
    Inds_intissue_record{intIndx} = find(ifintissuecurr);
    trecord(intIndx) = t;
    
    if max(Vncurr) < 1e-5
        break
    end
    if mod(intIndx,5) == 0
        save(filename); % save regularly in case sim crashes halfway
        toc
    end
end

disp('Main simulation ended');
toc

%% how minimum concentrations at tissue tip and gz properties vary as 
% tissue grows
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

%% save final data
save(filename)
