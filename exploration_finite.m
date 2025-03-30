%% Exploration LMIs

cvx_begin SDP
cvx_solver SDPT3
cvx_precision high


% Variables

variable tau(1)
variable u(nu,L)
variable gammae(1)
% variable D_des(nphi,nphi) semidefinite %comment this out to use a pre-determined D_des
% variable gammabar(1) semidefinite %comment this out for variable
% variable lambda(1)

% gammabar=5e3;
C3=(8*(sigmaw^2)*ntheta/T)*log(2*T*gamma_Aw*gamma_w + 2*T*T*(gamma_Au+1)*gammabar + lambda);

Ue=[];
for i=1:L
    Ue=blkdiag(Ue,u(:,i));
end

% Exploration energy and \bar{\gamma}>=\gamma_e

S_energy_1=[gammae, ones(1,L)*(Ue'); Ue*ones(L,1), gammae*eye(nu*L)];
% S_energy_1=[gammae, ones(1,L)*(Ue'); Ue*ones(L,1), eye(nu*L)];

S_energy_2=[gammabar, gammae; gammae, 1];

% S_exploration = S_1 - tau * S_2
% S_1
s1_11=(1-eps)*(Ue*Utilde' + Utilde*Ue' - Utilde*Utilde');
s1_12=zeros(nu*L,nphi);
s1_21=s1_12';
s1_22=-((2*(1-eps)/eps)*gammabar*Gamma_u) - ((2*(1-eps)/eps)*Gamma_w) + ((lambda/T)*eye(nphi)) - ((C2+C3)*D_des);
% s1_22=- ((C2+C3*gammabar)*D_des);

S_1=[s1_11, s1_12; s1_21, s1_22];

% S_2
s2_11= -eye(nu*L);
s2_12= V0';
s2_21=V0;
s2_22=Gamma_v-(V0*V0');

S_2=[s2_11, s2_12; s2_21, s2_22];

S_exp= S_1-(tau*S_2);

% optimiaztion problem

minimize gammae

S_energy_1>=0;
% S_energy_2>=0;
tau>=0;
S_exp>=0;
% D_des>=0;

cvx_end

disp(gammae^2);
% Update for the next iteration
Utilde=Ue;
% gammabar=gammae^2;


