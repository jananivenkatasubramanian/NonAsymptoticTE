green=[0.4660 0.6740 0.1880];
red=[0.6350 0.0780 0.1840];
blue=[0 0.4470 0.7410];
grey=[0.7 0.7 0.7];

% simtest10_3Dec.mat - error-vs-DO
% simtest11_4Dec.mat - error-vs-T

%% Gamma_e vs. Gamma_w for Non-stochastic


figure(1);hold on;
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xlim([1e10 1e15]);
% ylim([7e3 8.5e3]);
grid on;


p1=plot(simpars(:,1),simpars(:,1).*simpars(:,2),'LineWidth',2,'Color',blue);
p2=plot(simpars(:,1),simpars(:,1).*simpars(:,2),'o','LineWidth',1,'Color',blue);
% axis tight;

xlabel('Exploration time T','Interpreter','latex','FontSize',12);
% ylabel('Average input energy $\gamma_\mathrm{e}^2$','Interpreter','latex','FontSize',12);
ylabel('Input energy $T \gamma_\mathrm{e}^2$','Interpreter','latex','FontSize',12);

