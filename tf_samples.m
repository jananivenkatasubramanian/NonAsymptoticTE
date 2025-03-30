%% Construct \hat{V}_x and \hat{V}_\phi with prior estimates \hat{A}_0 and \hat{B}_0

V0=[]; %\hat{V}_0
for i=1:L
    vi=inv(exp(1i*2*pi*freqs(i))*eye(nx)-A0)*B0;
    vi_block=[vi;eye(nu)];
    V0=[V0,vi_block];
end

V_tr=[]; %V_\mathrm{tr}
for i=1:L
    vi=inv(exp(1i*2*pi*freqs(i))*eye(nx)-A)*B;
    vi_block=[vi;eye(nu)];
    V_tr=[V_tr,vi_block];
end

Vtilde=V_tr-V0;
VtVtH=Vtilde*Vtilde';
nVtil=norm(VtVtH);

%% \tilde{\Gamma}_v: Transfer functions V_\phi bounds: Sample-based
% \tilde{V} \tilde{V}^H <= \tilde{\Gamma}_\mathrm{v}

V_s=[];
cvx_begin SDP quiet
% cvx_solver mosek
cvx_precision high

variable Gamma_v((nx+nu),(nx+nu)) semidefinite

minimize trace(Gamma_v)
subject to
for i=1:nx:nx*Ns
    pa=params(i:i+nx-1,1:nx);
    pb=params(i:i+nx-1,nx+1:nphi);
    for j=1:L
        vtemp=inv(exp(1i*2*pi*freqs(j))*eye(nx)-pa)*pb;
        V_s=[V_s,[vtemp;eye(nu)]];
    end
    Gamma_v-((V_s-V0)*(V_s-V0)')>=0;
    V_s=[];
end
cvx_end


%% \Gamma_u: Sample-based bound g_star, Transiet error due to input
% gstar for Gamma_u

cvx_begin SDP quiet
cvx_precision high

variable gstar(1)

minimize gstar
subject to
for i=1:nx:nx*Ns
    pa=params(i:i+nx-1,1:nx);
    pb=params(i:i+nx-1,nx+1:nphi);
    gstar-(norm(pb)* ((1-norm(pa))^(-2)) )>=0;
end
cvx_end


