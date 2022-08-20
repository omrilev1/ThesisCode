% plot Results
close all; clear all; clc;

SNR_orig = -7:0.25:35;  % 8              % Direct channel SNR
snrLin = 10.^(SNR_orig/10);
profile2 = 1./(1 + (snrLin).^2);

idx = (SNR_orig < 30);

load('PPM_Profile2.mat');
totalEnergy_PPM = totalEnergy;
load('Linear_Profile2.mat');
totalEnergy_Linear = totalEnergy;

figure;subplot(1,2,1);
plot(SNR_orig(idx),10*log10(profile(idx)),'-.','LineWidth',2);hold on;
plot(SNR_orig(idx),10*log10(D_Linear(idx)),'LineWidth',2); hold on;
plot(SNR_orig(idx),10*log10(D_PPM(idx)),'--k','LineWidth',2); hold on;
xlabel('$$\tilde{E}/N$$ [dB]','FontSize',14,'Interpreter','Latex'); ylabel('Distortion [dB]','FontSize',14,'Interpreter','Latex');
lgd = legend({'Profile','Linear','PPM'},'FontSize',10,'TextColor','Black','Location','Best');
grid on; grid minor;
xlim([0 30]);

subplot(1,2,2);
semilogy(SNR_orig,totalEnergy_Linear,'LineWidth',2); hold on;
semilogy(SNR_orig,totalEnergy_PPM,'-.k','LineWidth',2);hold on;
xlabel('$$\tilde{E}/N$$ [dB]','FontSize',14,'Interpreter','Latex');
ylabel('Accumulated Energy/$$\tilde{E}$$','FontSize',14,'Interpreter','Latex');
xlim([0 30]); ylim([0.6 2.4])
lgd = legend({'Linear','PPM'},'FontSize',10,'TextColor','Black','Location','Best');
grid on; grid minor;
% xlim([7.5 30]);
% sgtitle({'Quadratic Profile - Scalar Simulation'},'FontSize',14);

% load('PPM_Profile3.mat');
% load('Linear_Profile3.mat');
% subplot(1,2,2);
% semilogy(SNR,profile,'-.','LineWidth',2.5);hold on;
% semilogy(SNR,D_Linear,'LineWidth',2.5); hold on;
% semilogy(SNR,D_PPM,'--k','LineWidth',2.5); hold on;
% xlabel('$$\tilde{E}/N$$ [dB]','FontSize',14,'Interpreter','Latex');
% lgd = legend({'Profile','Scalar','PPM'},'FontSize',10,'TextColor','Black','Location','Best');
% grid on; grid minor;
% title({'Third Order Profile', strcat('\Delta = ',num2str(Delta), ', \alpha = ',num2str(alpha))},'FontSize',14);
%


%% Infinite Dimension
SNR_orig = -7:0.25:30;  % 8              % Direct channel SNR
snrLin = 10.^(SNR_orig/10);
profile = 1./(1 + (snrLin).^2);


load('PPMInfDim_Profile2.mat');
totalEnergy_PPM = totalEnergy;
load('LinearInfDim_Profile2.mat');
totalEnergy_Linear = totalEnergy;
load('PpmInfDim_Anal_Profile2.mat');
totalEnergy_PPM_Anal = totalEnergy;


% Tuncel graph
load('Tuncel_InfDim_Performace2.mat');
SNR_Tuncel = SNR;  % Direct channel SNR
totalEnergy_Tuncel = totalEnergy;
D_Tuncel = D_tot;

%% Results from Tuncel+Banniasadi Paper

%% Plot results
figure;
plot(SNR_orig,10*log10(profile),'-.','LineWidth',2);hold on;
plot(SNR_Tuncel,10*log10(D_Tuncel),':','LineWidth',2); hold on;
plot(SNR_orig,10*log10(D_Linear),'-','LineWidth',2); hold on;
plot(SNR_orig,10*log10(D_PPM),'--k','LineWidth',2); hold on;
plot(SNR_orig,10*log10(D_PPM_Anal),'-*','LineWidth',0.25); hold on;
xlabel('$$\tilde{E}/N$$ [dB]','FontSize',14,'Interpreter','Latex'); ylabel('Distortion [dB]','FontSize',14,'Interpreter','Latex');
lgd = legend({'Profile','Baniasadi-Tuncel','Linear','PPM Empiric','PPM Analytic Bound'},'FontSize',10,'TextColor','Black','Location','Best');
grid on; grid minor;
xlim([-5 30]);

figure;
semilogy(SNR_Tuncel,totalEnergy_Tuncel,'--','LineWidth',2);hold on;
semilogy(SNR_orig,totalEnergy_Linear,'-.','LineWidth',2); hold on;
semilogy(SNR_orig,totalEnergy_PPM_Anal,'LineWidth',2);hold on;
semilogy(-7:0.25:30,totalEnergy_PPM,':','LineWidth',2);hold on;
xlabel('$$\tilde{E}/N$$ [dB]','FontSize',14,'Interpreter','Latex');
ylabel('Accumulated Energy/$$\tilde{E}$$','FontSize',14,'Interpreter','Latex');
lgd = legend({'Baniasadi-Tuncel','Linear','PPM Analytic','PPM Empiric'},'FontSize',10,'TextColor','Black','Location','Best');
grid on; grid minor;
xlim([-5 30]); ylim([0.5 2.4])

% sgtitle({'Quadratic Profile - Infinite Dimension Simulation'},'FontSize',14);
