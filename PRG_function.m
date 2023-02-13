function [M,V,P0,L,Pnz,bins_nz] = PRG_function(raster,n_steps,varargin)

% This function implements the Phenomenological Renormalization Group
% (PRG) method, introduced by Meshulam et al. (2018, 2019). It also can
% implement the extended version of PRG based on connectivity.
%
% Within this method, the collective activity is iteratively coarse-grained 
% by grouping maximally correlated variables (or maximally coupled variables 
% in the connectivity-based case). At each coarse-graining step 
% k=0,1,…,kmax, clusters of size K=2^k are built, resulting in a system of 
% N/K coarse-grained variables and successively ignoring degrees of
% freedom. This code computes several observables of the coarse-grained 
% variables as a function of K.   
%
% Inputs: 
% - raster: binary data, T-by-N matrix, where T is the number of samples
% and N is the number of variables.
% - n_steps: number of coarse-graining steps.
% - fc_based: = 1 FC-based PRG method; = 0 connectivity-based PRG
% optional:
% - [...] = PRG_function(raster,n_steps,C) including a third input to the
% function indicates that the method is connectivity-based.
% C needs to be a N-by-N symmetrical matrix.
%
% Outputs:
% - M: Mean activity of coarse-grained variables as a function of K.
% - V: Variance of coarse-grained activity as a function of K.
% - P0: Probability of silence activity as a function of K.
% - L: eigenspectrum of the covariance matrix of coarse-grained activity 
% for each cluster of size K .
% - Pnz: distribution of non-zero activity.
% - bins_nz: bins associated to Pnz.
%
% Adrián Ponce-Alvarez 2022
%
% Refs:
% Meshulam, L., Gauthier, J.L., Brody, C.D., Tank, D.W. & Bialek, W. 
% Coarse graining, fixed points, and scaling in a large population of neurons. 
% Phys. Rev. Lett. 123, 178103 (2019).
%
% Meshulam, L., Gauthier, J.L., Brody, C.D., Tank, D.W. & Bialek, W. 
% Coarse-graining and hints of scaling in a population of 1000+ neurons. 
%arXiv, 1812.11904 (2018).

%--------------------------------------------------------------------------

% data size:    
[T,N] = size(raster);

% bins for Pnz:
nbins = 100;
bins_nz = linspace(0,1,nbins); % bins for Pnz
dnz = bins_nz(2)-bins_nz(1);

% initialize
M   = zeros(1,n_steps);  % mean
V   = zeros(1,n_steps);  % variance
P0  = zeros(1,n_steps);  % silence prob.
Pnz = nan(nbins,n_steps); % Non-zero activity prob.

X = raster;
M(1) = mean(mean(X));
V(1) = mean(var(X));
P0(1) = sum(X(:)==0)/length(X(:));
sets=cell(1,n_steps);

if nargin > 2
   disp('Connectivity-based PRG method...') 
   fc_based = 0;
   Connectivity = varargin{1};
else
    disp('FC-based PRG method...')
   fc_based = 1;
end

% % flag to return coarse-grained variables
% if nargin > 3
%     flag_return_CG = 1;
%     
%     CGact = cell(1,last_steps);
% else
%     flag_return_CG = 0;
%     CGact = [];
% end



% Iterative coarse-grain:
for k = 2:n_steps
    
    fprintf('GC step: %g \n',k)    
    w=whos('X');
    if k==2
    fprintf('calculating covariance matrix of data of %g GB ...\n',w.bytes/1e9)  
    end

    if fc_based % FC-based PRG method
        FC = corr(X);
        FC = single(FC);
    else  % connectivity-based PRG
        FC = Connectivity;    
    end

    FC( 1 : size(FC,1)+1 : end ) = 0; % diag = 0
    maxfc = max(FC); % max of each column of the FC

    K = 2^(k-1);
    Nk = floor(N/K);
    Y = zeros(T,Nk);
    temp = zeros(2,Nk);

        % Get coarse-grained time series:
        for n=1:Nk

            [~,j]=max(maxfc);
            j = j(1);
            [~,i]=max( FC(:,j) );
            i = i(1);
            % combine:
            Y(:,n) = single(X(:,i) + X(:,j));
            % exclude:
            FC(:,j) = -1;
            FC(:,i) = -1;
            FC(j,:) = -1;
            FC(i,:) = -1;
            maxfc(i) = -1;
            maxfc(j) = -1;
            % store participating units:
            temp(:,n) = [i;j];

        end

        sets{k} = temp;

% observables:
M(k) = mean(mean(Y));
V(k) = mean(var(Y));
P0(k) = sum(Y(:)==0)/length(Y(:));
c = hist(Y(Y>0)/K,bins_nz); 
Pnz(:,k) = c/sum(Y(:)>0)/dnz;

% % return GC rasters if asked
% if flag_return_CG
%    if k >= (n_steps - last_steps + 1) 
%    CGact{k - (n_steps - last_steps)} = Y;   
%    end
% end
X = Y;
clear FC


    % if connectivity-based method:
    if ~fc_based
        % new connectivity:
        Ck = zeros(Nk);
        for p = 1:Nk-1
            i = temp(1,p);
            j = temp(2,p);
            for q = p+1:Nk
            n = temp(1,q);
            m = temp(2,q);
            Ck(p,q) = 1/4*(Connectivity(i,n) + Connectivity(i,m) + Connectivity(j,n) + Connectivity(j,m));
            end
        end
        Connectivity = Ck; % update connectivity
    end

end


% retreive neuron labels of clusters:
Set = sets{2};
for k = 3:n_steps
J = sets{k};
m = J(1,:);
n = J(2,:);
neurons = [Set(:,m);Set(:,n)];
sets{k} = neurons;
Set = neurons;
end

% Eigen-decomposition clusters:
disp('eigen-decomposition...')
maxK = 2^(n_steps-1); % largest K
L = nan(maxK,n_steps);
for k=2:n_steps
    Nk = size(sets{k},2);
    setk = sets{k};
    lambda = zeros(size(setk));
    for c = 1:Nk
        ii = setk(:,c);
        X = raster(:,ii);
        C = cov(X);
        [~,d]=eig(C);
        ds = sort(diag(d),'descend');
        lambda(:,c) = ds;
    end
    K = 2^(k-1);
    mL = mean(lambda,2);
    L(1:K,k) = mL;
end

return
