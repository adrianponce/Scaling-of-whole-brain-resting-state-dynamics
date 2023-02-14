function [Output] = metropolis_spin_model(C,D,T,n_steps,CBPRG,varargin)

% Simulates the spin model using the Metropolis algorithm.
%
% Inputs:
% - C : connectivity matrix.
% - D : distance matrix.
% - T : Temperature.
% - n_steps : number of coarse-graining steps.
% - CBPRG : if = 1 we use the connectivity-based PRG, 
%           if = 0 we use the FC-based PRG
% Optional inputs:
% - [...] = metropolis_spin_model(C,D,T,n_steps,CBPRG,Tmax,NSUBSIM,num_initialcond)
% with :
%    - Tmax : number of samples (x N)
%    - NSUBSIM : number of batches
%    - num_initialcond : number of initial conditions
%
% Outputs:
% - Magnetization
% - GrIsing : correlation function (correlations vs. distance)
% - Variance : Variance of coarse-grained variables.
% - Psilence : Silence prob. of coarse-grained variables.
% - Spectrum : Eigenspectrum of coarse-grained variables.
% - Heat_capacity : Heat capacity of the spin model.
% - Susceptibility : Suceptibility of the spin model.
%
% Adrián Ponce-Alvarez 2022
%
%--------------------------------------------------------------------------


N = size(C,1); % nb. of spins

% simulation parameters:
if nargin > 5
Tmax=varargin{1}*N;
NSUBSIM=varargin{2}; % nb of batches
num_initialcond = varargin{3};    
else % default:
Tmax=10000*N;
NSUBSIM=5; % nb of batches
num_initialcond = 5;
end


maxK = 2^(n_steps-1); % largest K
indices = find( triu(ones(N)) - eye(N) ); % indices of upper triangle of matrices of size N-by-N

% distances to calculate the correlation function:
R = D(indices); % all distances between nodes
rmax=max(R);
NR = 400;
dr=rmax/NR;
bins = 0:dr:rmax;
nbins = length(bins);
binsnew = [bins bins(end)+dr];
distances = nan(nbins,1);
for w = 1:nbins  
distances(w) = binsnew(w)+(binsnew(w+1)-binsnew(w))/2; % distances
end



G = 1/T; %inverse temperature
wC = G*C; % Scaled interactions    

% initialize:
HeatC = zeros(1,num_initialcond);
Suscept = zeros(1,num_initialcond);
Magnet = zeros(1,num_initialcond);
Gr_funct = zeros(nbins,num_initialcond);

Var_T  = zeros(n_steps,num_initialcond);
Psil_T  = zeros(n_steps,num_initialcond);
Spect_T = nan(maxK,n_steps,num_initialcond);

% Simulation:
%--------------------------------------------------------------------------

for INIT = 1:num_initialcond %Loop across initial conditions

    % Metropolis algorithm
    % Transient regime (variables evolve but are not stored):
    Spin=rand(N,1);
    Spin=2*ceil(Spin-0.5)-1; %random binary vector of +/1 1
    JV=wC*Spin;       
    for t=1:500000
        ind=ceil(N*rand);              % pick one spin randomly
        DeltaE =-2*Spin(ind)*JV(ind);  % energy difference
        if rand<exp(DeltaE)            % update prob.
            JV=JV-2*Spin(ind)*wC(:,ind);
            Spin(ind)=-Spin(ind);
        end
    end
    
    Gr=zeros(nbins,NSUBSIM);
    heatcsub = 0;
    susceptsub = 0;
    magn = 0;
    
    Var    = zeros(n_steps,NSUBSIM);
    Psil   = zeros(n_steps,NSUBSIM);
    Spect  = nan(maxK,n_steps,NSUBSIM);

    
    for nsub=1:NSUBSIM  % loop over batches
        fprintf('run %g over %g (init: %g)\n',nsub,NSUBSIM,INIT)
        Pattern=zeros(floor(Tmax/N),N,'single');
        Energy=zeros(floor(Tmax/N),1);
        M=zeros(floor(Tmax/N),1);
        
        % Metropolis:
        count = 0;
        for t=1:Tmax
            for n=1:N
                ind=ceil(N*rand);
                DeltaE =-2*Spin(ind)*JV(ind);
                if rand<exp(DeltaE)
                    JV=JV-2*Spin(ind)*wC(:,ind);
                    Spin(ind)=-Spin(ind);
                end
            end
            if mod(t,N)==0
            count = count + 1;    
            Pattern(count,:)=single(Spin); %single precision to save memory
            Energy(count)=Spin'*wC*Spin;
            M(count)=mean(Spin);
            end
        end
        
        m = mean(Pattern); % mean spin value <spin>
        I=find(m==1 | m==-1); % spins that didn't fluctuate
        c = corr(Pattern);
        c(I,:) = 0;
        C(:,I) = 0;        
        
        heatcsub = heatcsub + var(Energy)/NSUBSIM; % heat capacity
        susceptsub = susceptsub + var(M)/NSUBSIM;  % susceptibility
        magn = magn + mean(M)/NSUBSIM;             % magnetization
        
        % Correlation function:       
        Z = c(indices); % vector of all pairwise correlations
        x = nan(nbins,1);
        
        % g(r) : correlation as a function of distance between spins
        for w = 1:nbins  
            ii = find( R>=binsnew(w) & R<binsnew(w+1) );
            y = Z(ii);
            if ~isempty(ii)
            x(w) = nanmean(y);
            end
        end

        Gr(:,nsub) = x; % correlation function       
        
        % Phenomenological Renormalization Group:
        raster = (Pattern+1)/2; % transfor to (0,1) variables
        if CBPRG == 1 % if connectivity-based PRG    
        [~,V,P0,S] = PRG_function(raster,n_steps,n_steps,wC);
        else % if FC-based PRG
        [~,V,P0,S] = PRG_function(raster,n_steps,n_steps);
        end
        Var(:,nsub) = V;     % Variance
        Psil(:,nsub) = P0;   % Prob. silence
        Spect(:,:,nsub) = S; % Eigenspectrum
        clear Pattern
                
    end
    
    HeatC(INIT) = heatcsub;           % heat capacity
    Suscept(INIT) = susceptsub;       % susceptibility
    Magnet(INIT) = magn;              % magnetization  
    Gr_funct(:,INIT) = nanmean(Gr,2); % corrlation function
            
    % Renormalization observables:
    Var_T(:,INIT) = mean(Var,2);
    Psil_T(:,INIT) = mean(Psil,2);
    Spect_T(:,:,INIT) = nanmean(Spect,3);
        
    
end


    Output.Magnetization  = mean(Magnet);
    Output.GrIsing        = nanmean(Gr_funct,2);   
    Output.Variance       = nanmean(Var_T,2);
    Output.Psilence       = nanmean(Psil_T,2);
    Output.Spectrum       = nanmean(Spect_T,3);
    Output.Heat_capacity  = mean(HeatC);
    Output.Susceptibility = mean(Suscept);
    Output.distances      = distances;
    



