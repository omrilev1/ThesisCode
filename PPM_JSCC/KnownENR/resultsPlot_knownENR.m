close all; clear all; clc;

%% plot Results - Uniform Source
data = load('UniformOptSDR'); 
ENR_AnalogPPM = data.ENR; 
analogPPM_optEmpiricMSE =  data.final_MSE; 
analogPPM_optEmpiricBeta =  data.final_Beta;
analogPPM_optExactMSE =  data.optSDR_Analytic_Exact; 
analogPPM_optExactbeta =  data.opt_beta; 

data = load('UniformBurnashevOptSDR'); 
ENR_Burnashev = data.ENR; 
Burnashev_optEmpiricMSE =  data.final_MSE; 
Burnashev_optEmpiricBeta =  data.final_Beta;
Burnashev_optExactMSE =  data.optSDR_Analytic_Exact; 
Burnashev_optExactbeta =  data.opt_beta; 

data = load('UniformTuncelOptSDR'); 
ENR_Tuncel = data.ENR; 
Tuncel_optEmpiricMSE =  data.final_MSE; 
Tuncel_optEmpiricBeta =  data.final_Beta;
Tuncel_optExactMSE =  data.optSDR_Analytic_Exact; 
Tuncel_optExactbeta =  data.opt_beta; 

% data = load('MAP-MMSE comparison results\Uniform\vectorBetaPPM_Rect_Uniform.mat');
% analogPPM_MMSE_optEmpiricMSE = data.final_MSE_MMSE; 


figure;
hold all;
% plot(ENR_AnalogPPM,10*log10(1/12./analogPPM_MMSE_optEmpiricMSE),'--','LineWidth',1.5);
plot(ENR_AnalogPPM,10*log10(1/12./analogPPM_optEmpiricMSE),'-.ko','LineWidth',1.5);
plot(ENR_AnalogPPM,10*log10(1/12./analogPPM_optExactMSE),':','LineWidth',2);
plot(ENR_Tuncel,10*log10(1/12./Tuncel_optEmpiricMSE),'-bp','LineWidth',1.5);
plot(ENR_Burnashev,10*log10(1/12./Burnashev_optEmpiricMSE),'--rd','LineWidth',1.5);
xlabel('ENR [dB]'); ylabel('SDR [dB]'); 
legend({% 'Analog PPM: Empirical MMSE',
    'Analog PPM: Empirical MAP','Analog PPM: Lower Bound','Sevinc Tuncel: Empirical','Burnashev: Empirical'},'Location','northwest','FontSize',12);
grid on; grid minor;
% figure;
% semilogy(ENR_AnalogPPM,analogPPM_optEmpiricBeta,'-.ko','LineWidth',1.5);hold on;
% semilogy(ENR_AnalogPPM,analogPPM_optExactbeta,':','LineWidth',2);
% xlabel('ENR [dB]'); ylabel('Optimal Beta');
% legend({'Analog PPM Optimal \beta: Empirical','Analog PPM Optimal \beta: Bound'},'Location','northwest','FontSize',12);
% grid on; grid minor;
ENRlin = 10.^(ENR_AnalogPPM/10); 
opt_beta = (312*sqrt(pi))^(1/3) * (ENRlin).^(-5/6) .* exp(ENRlin/6);

figure;
semilogy(ENR_AnalogPPM,opt_beta,'-.ko','LineWidth',1.5);hold on;
xlabel('ENR [dB]'); ylabel('Optimal Beta');
legend({'Analog PPM \beta'},'Location','northwest','FontSize',12);
grid on; grid minor;


%% plot Results - Gaussian Source
data = load('GaussianOptSDR'); 
ENR_AnalogPPM = data.ENR; 
analogPPM_optEmpiricMSE =  data.final_MSE; 
analogPPM_optEmpiricBeta =  data.final_Beta;
analogPPM_optExactMSE =  data.optSDR_Analytic_Exact; 
analogPPM_optExactbeta =  data.opt_beta; 

data = load('GaussianBurnashevOptSDR'); 
ENR_Burnashev = data.ENR; 
Burnashev_optEmpiricMSE =  data.final_MSE; 
Burnashev_optEmpiricBeta =  data.final_Beta;
Burnashev_optExactMSE =  data.optSDR_Analytic_Exact; 
Burnashev_optExactbeta =  data.opt_beta; 

% data = load('GaussianTuncelOptSDR'); 
% ENR_Tuncel = data.ENR; 
% Tuncel_optEmpiricMSE =  data.final_MSE; 
% Tuncel_optEmpiricBeta =  data.final_Beta;
% Tuncel_optExactMSE =  data.optSDR_Analytic_Exact; 
% Tuncel_optExactbeta =  data.opt_beta; 


data = load('GaussianTuncelOptSDR_SigmaOptimalCompander'); 
ENR_Tuncel = data.ENR; 
Tuncel_SigmaOptEmpiricMSE =  data.final_MSE; 
Tuncel_SigmaOptEmpiricBeta =  data.final_Beta;
Tuncel_SigmaOptExactMSE =  data.optSDR_Analytic_Exact; 
Tuncel_SigmaOptExactbeta =  data.opt_beta; 


% data = load('MAP-MMSE comparison results\Gaussian\vectorBetaPPM_Rect_Gaussian.mat');
% analogPPM_MMSE_optEmpiricMSE = data.final_MSE_MMSE; 
% Calculate analog PPM bound for the whole region 
% optimize the lower bound
% exact upper bound
analogPPM_optExactMSE = zeros(size(ENR_Tuncel));
opt_beta = zeros(size(ENR_Tuncel));
ENR_Tuncel_Lin = 10.^(ENR_Tuncel/10); 
for i=1:length(ENR_Tuncel)
    
    beta_opt = (13/8)^(1/3) * (ENR_Tuncel_Lin(i))^(-5/6) * exp(ENR_Tuncel_Lin(i)/6);
    beta_vec = beta_opt; % beta_opt*(0.5:1/64:1.5);
    
    D_t = ((13/8) + sqrt(2./beta_vec).* (sqrt(2*beta_vec.*ENR_Tuncel_Lin(i)) - 1).*exp(-ENR_Tuncel_Lin(i) .* (1 - 1./sqrt(2*beta_vec.*ENR_Tuncel_Lin(i))).^2)) ./ ((sqrt(beta_vec.*ENR_Tuncel_Lin(i)) - 1/sqrt(2)).^4) ...
        + exp(-beta_vec.*ENR_Tuncel_Lin(i))./(beta_vec.^2);
    
    D_L = 2*beta_vec.*sqrt(ENR_Tuncel_Lin(i)).*exp(-ENR_Tuncel_Lin(i)/2) .* (1 + 3*sqrt(2*pi./ENR_Tuncel_Lin(i)) + 12*exp(-1)./(beta_vec.*sqrt(ENR_Tuncel_Lin(i))) ...
        + 8*exp(-1)./(sqrt(8*pi)*beta_vec) + sqrt(8./(pi*ENR_Tuncel_Lin(i))) + 12^(3/2) * exp(-3/2) ./(beta_vec.*sqrt(32*pi*ENR_Tuncel_Lin(i))));
    
    [analogPPM_optExactMSE(i),optIdx] = min(D_t + D_L);
    opt_beta(i) = beta_vec(optIdx);

    
end

figure;
hold all;
% plot(ENR_AnalogPPM,10*log10(1./analogPPM_MMSE_optEmpiricMSE),'--','LineWidth',1.5);
plot(ENR_AnalogPPM,10*log10(1./analogPPM_optEmpiricMSE),'-.ko','LineWidth',1.5);
plot(ENR_Tuncel,10*log10(1./Tuncel_SigmaOptEmpiricMSE),'-bp','LineWidth',1.5);
plot(ENR_Burnashev,10*log10(1./Burnashev_optEmpiricMSE),'--rd','LineWidth',1.5);
plot(ENR_Tuncel,10*log10(1./analogPPM_optExactMSE),':','LineWidth',1.5);
xlabel('ENR [dB]'); ylabel('SDR [dB]'); legend({% 'Analog PPM: Empirical MMSE',
    'Analog PPM: Empirical MAP',...
    'Sevinc Tuncel: Empirical','Burnashev: Empirical','Analog PPM: Lower Bound'},'Location','northwest','FontSize',12);
grid on; grid minor;
xlim([4 15]); ylim([-15 50])
figure; 
semilogy(ENR_AnalogPPM,analogPPM_optEmpiricBeta,'-.ko','LineWidth',1.5);hold on;
xlabel('ENR [dB]'); ylabel('Optimal Beta');
legend({'Analog PPM \beta'},'Location','northwest','FontSize',12);
grid on; grid minor;

 

ENR_TuncelLin = 10.^(ENR_Tuncel/10); 
D_BoundTuncel = exp(1.7) * exp(-ENR_TuncelLin/3); 
figure;
semilogy(ENR_Tuncel,D_BoundTuncel,'-.ko','LineWidth',1.5); hold on; 
semilogy(ENR_Tuncel,Tuncel_optEmpiricMSE,'-bp','LineWidth',1.5);hold on; 
semilogy(ENR_Tuncel,Tuncel_SigmaOptEmpiricMSE,'-rh','LineWidth',1.5);hold on; 
semilogy(ENR_Burnashev,Burnashev_optEmpiricMSE,'--','LineWidth',1.5);hold on; 
xlabel('ENR [dB]'); ylabel('Distortion [dB]'); legend({'Tuncel Bound','Tuncel Empiric','Tuncel Empiric Sigma Opt','Burnashev'},'Location','northwest','FontSize',12);
grid on; grid minor;
