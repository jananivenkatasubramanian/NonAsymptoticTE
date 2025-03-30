 %% Main --- Non-stoch targeted exploration
clear all;
clc;
simpars=[];

%%
% fsel=10;
delta=0.05;  % probability of violation
eps=0.5; % young's inequality parameter
lambda=1e0; % regularization parameter

run initialize_finite.m

% Set initial bound on estimate D0 - with structure so that  bounds on \hat{A}_0 
% and \hat{B}_0 can be computed 
D0tilde=1e3*eye(nphi);
D0=kron(D0tilde,eye(nx));
D0inv=inv(D0);

% initial estimate \hat{A}_0 and \hat{B}_0
theta0=thetatr + 4e-4;
theta0_check=(theta0-thetatr)'*D0*(theta0-thetatr);
if theta0_check>1
    disp("Initial uncertainty bound not satisfied!");
    disp(theta0_check);
else
    disp(theta0_check);
end
p0=reshape(theta0,[nx,nphi]);
A0=p0(:,1:nx);
B0=p0(:,nx+1);

% Generate \theta samples to derive scenario bounds
thetas=[];
params=[];
n=0;
Ns=1000; % set number of samples based on required confidence
cbar=chi2inv(1-delta,(nx*(nx+nu))); %to scale D0 - to get samples from normrnd
while n<Ns
    t1=(mvnrnd(theta0,D0inv/cbar))';
    t2=(theta0-t1)'*D0*(theta0-t1);
    if t2<=1
        thetas=[thetas,t1];
        paramt=reshape(t1,[nx,nphi]);
        params=[params;paramt];
        n=n+1;
    end
end

%% Tranfer function estimate, bounds

run tf_samples.m
%% Construct transfer matrices and Scenario bounds \tilde{\Gamma}_v, gstar, \gamma_Au, \gamma_Aw
% Exploration time
T=1e9;

% mul=1;
sw=[0.01]; %size 6

% Disturbance proxy variance
sigmaw=sw(1);

% Desired bound D_des
D_des=1e4*eye(nphi);

% gamma_w, Gamma_u, Gamma_w, C_1, C_2, C_3
run constants_finite.m


%% Set initial candidate \tilde{U} and \bar{\gamma}

gammabar=1e4;
Utilde=[];
ut=1e3*ones(nu,L);
for i=1:L
    Utilde=blkdiag(Utilde,ut(:,i));
end

%% Run this section multiple times

run exploration_finite.m
disp(gammae^2);

%% ... or run this loop
% Non-asymptotic exploration: run iterations to reduce suboptimality from convex-relaxation

gtemp=Inf;
for i=1:20
    run exploration_finite.m
    if(abs((gtemp-gammae)/gtemp))<1e-2
        break;
    end
    gtemp=gammae;
end
disp(Ue);
disp(gammae^2);

%%
simpars=[simpars;T,gammae^2];

