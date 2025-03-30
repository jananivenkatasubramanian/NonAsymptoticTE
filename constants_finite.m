%% Sample-based bound on \|A_w\| <= \gamma_Aw

cvx_begin SDP quiet
cvx_precision high

variable gamma_Aw(1)

for i=1:nx:nx*Ns
    pa=params(i:i+nx-1,1:nx);
%     Awi=eye(T*nx);
%     for row = 2:T
%         for col = 1:row-1
%             Awi((row-1)*nx+1:row*nx, (col-1)*nx+1:col*nx) = pa^(row-col);
%         end
%     end
    nAwi=(1-norm(pa)^T)/(1-norm(A));
%     gamma_Aw - (norm(Awi))^2 >=0;
    gamma_Aw - (nAwi)^2 >=0;
end
cvx_end

%% Sample-based bound on \|A_u\| <= \gamma_Au


cvx_begin SDP quiet
cvx_precision high

variable gamma_Au(1)

for i=1:nx:nx*Ns
    pa=params(i:i+nx-1,1:nx);
    pb=params(i:i+nx-1,nx+1:nphi);
%     Aui=zeros(T*nx,T*nu);
%     for row=1:T
%         for col=1:row
%             Aui((row-1)*nx+1:row*nx, (col-1)*nu+1:col*nu) = pa^(row-col) * pb;
%         end
%     end
    nAui=norm(pb)*(1-norm(pa)^T)/(1-norm(A));
    gamma_Au - (nAui)^2 >=0;
end
cvx_end

%% Define \Gamma_u

Gamma_u=(4*L*(gstar^2)/T)*eye(nphi);

% \Gamma_w: Effect of disturbance
% Bound \|W\|^ \leq \gamma_w

gamma_w=(sigmaw^2)*( (1+log(4))*T*nx + 4*log(1/delta) );

%% \Gamma_w
Gamma_w=(L*gamma_Aw*gamma_w/T)*eye(nphi);

%% Upper bound on D_T
% Constants C_1, C_2 and C_3

C1=(4*(sigmaw^2))*( log((5^(2*nx))/(delta^2)) - ntheta*log(lambda) );

% \bar{\theta}
theta_bar=norm(theta0)+norm(chol(D0inv));

C2=(1/T)*(2*C1 + 2*lambda*(theta_bar^2));

% C3= 8*T*nx*nphi*(sigmaw^2)*gamma_Au;


