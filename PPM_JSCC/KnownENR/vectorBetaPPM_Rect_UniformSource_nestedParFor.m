%% Analog PPM Performance simulation
close all; clear all; clc;

ENR = (4:1:16);
ENRlin = 10.^(ENR/10);

Nrun = 2*1e4;
final_MSE = zeros(1,length(ENR));
final_Beta = zeros(1,length(ENR));

for i=1:length(ENR)
    % pick beta to find the optimal point of the curve
    beta_opt = (312*sqrt(pi))^(1/3) * (ENRlin(i))^(-5/6) * exp(ENRlin(i)/6);
    beta = max(1,beta_opt*(0.2:(1/64):1.25));
    MSE = zeros(size(beta));
    currENRlin = ENRlin(i);
    for beta_idx = 1:length(beta)
        curr_beta = beta(beta_idx);
        currMSE = calc_perf_uniform_ppm(curr_beta,Nrun,currENRlin);
        MSE(beta_idx) = currMSE;
    end
    %     [optVal,optIdx] = min(abs(MSE - 13/8 ./ (beta.*currENRlin).^2));
    [optVal,optIdx] = min(MSE);
    final_MSE(i) = MSE(optIdx);
    final_Beta(i) = beta(optIdx);
    disp(strcat('Finished ENR = ',num2str(ENR(i))));
end

optSDR_Analytic_Approx = 0.072 * (ENRlin).^(-1/3) .* exp(-ENRlin/3);
optBeta_Analytic = (312*sqrt(pi))^(1/3) * (ENRlin).^(-5/6) .* exp(ENRlin/6);

% Tuncel Bound 
D_Tuncel = exp(-0.9458) * exp(-ENRlin/3); 
SDR_Tuncel = 10*log10(1./D_Tuncel);

% optimize the lower bound
% exact upper bound
optSDR_Analytic_Exact = zeros(size(ENR));
opt_beta = zeros(size(ENR));
for i=1:length(ENR)
    beta_opt = (312*sqrt(pi))^(1/3) * (ENRlin(i))^(-5/6) * exp(ENRlin(i)/6);
    beta_vec =  max(1,beta_opt*(0.0125:(1/128):5));
    D_t = (13/8)./(beta_vec .* ENRlin(i)).^2;
    P_t = beta_vec .* sqrt(ENRlin(i)) .* exp(-1*ENRlin(i)/2) / (16*sqrt(pi));
    [optSDR_Analytic_Exact(i),optIdx] = min(D_t .* (1 + (16/13) * sqrt(ENRlin(i)/2) .* exp(-ENRlin(i)/4)) + ...
        (1/6)*(1 + 2./beta_vec + (2./beta_vec).^2) .* P_t .* (1 + 4 * sqrt(pi./ENRlin(i)) + sqrt(3/4) * exp(-1*ENRlin(i)/6)./sqrt(ENRlin(i))));
    opt_beta(i) = beta_vec(optIdx);
end

save('UniformOptSDR.mat','ENR','final_MSE','optSDR_Analytic_Exact','final_Beta','opt_beta');


figure;subplot(211);
hold all;
plot(ENR,10*log10(1/12./final_MSE),'-.ko','LineWidth',1.5);
plot(ENR,10*log10(1/12./optSDR_Analytic_Exact),'--ms','LineWidth',1.5);
xlabel('ENR [dB]'); ylabel('SDR [dB]'); legend({'Empirical','Lower Bound'},'Location','northwest','FontSize',12);
grid on; grid minor;
subplot(212);
semilogy(ENR,final_Beta,'-.ko','LineWidth',1.5);hold on;
semilogy(ENR,opt_beta,'--ms','LineWidth',1.5);
xlabel('ENR [dB]'); ylabel('Optimal Beta');
grid on; grid minor;


figure;subplot(211);
hold all;
plot(ENR,10*log10(1/12./final_MSE),'-.ko','LineWidth',1.5);
plot(ENR,10*log10(1/12./optSDR_Analytic_Approx),'--ms','LineWidth',1.5);
plot(ENR,10*log10(1/12./optSDR_Analytic_Exact),'->','LineWidth',1.5);
xlabel('ENR [dB]'); ylabel('SDR [dB]'); legend({'Empiric Opt SDR','Lower Bound Opt SDR - Approx','Lower Bound Opt SDR - Exact'},'Location','northwest','FontSize',12);
title('Optimal SDR [dB]');
grid on; grid minor;
subplot(212);
semilogy(ENR,final_Beta,'-.ko','LineWidth',1.5);hold on;
semilogy(ENR,optBeta_Analytic,'--ms','LineWidth',1.5);
xlabel('ENR [dB]'); ylabel('Optimal Beta [log]');legend({'Empiric Optimal \beta','Optimal \beta - Bound'},'Location','northwest','FontSize',12);
title('Optimal \beta [log]');
grid on; grid minor;


figure;
hold all;
plot(ENR,10*log10(1/12./final_MSE),'-.ko','LineWidth',1.5);
plot(ENR,10*log10(1/12./optSDR_Analytic_Approx),'--ms','LineWidth',1.5);
plot(ENR,10*log10(1/12./optSDR_Analytic_Exact),'->','LineWidth',1.5);
xlabel('ENR [dB]'); ylabel('SDR [dB]'); legend({'Empiric Opt SDR','Lower Bound Opt SDR - Approx','Lower Bound Opt SDR - Exact'},'Location','northwest','FontSize',12);
title('Optimal SDR [dB]');
grid on; grid minor;
figure;
semilogy(ENR,final_Beta,'-.ko','LineWidth',1.5);hold on;
semilogy(ENR,optBeta_Analytic,'--ms','LineWidth',1.5);
xlabel('ENR [dB]'); ylabel('Optimal Beta [log]');legend({'Empiric Optimal \beta','Optimal \beta - Bound'},'Location','northwest','FontSize',12);
title('Optimal \beta [log]');
grid on; grid minor;



figure;subplot(211);
hold all;
plot(ENR,10*log10(1/12./final_MSE),'-.ko','LineWidth',1.5);
plot(ENR,10*log10(1/12./D_Tuncel),'--bd','LineWidth',1.5);
plot(ENR,10*log10(1/12./optSDR_Analytic_Exact),'--ms','LineWidth',1.5);
xlabel('ENR [dB]'); ylabel('SDR [dB]'); legend({'Empirical','Tuncel High ENR Approx','Lower Bound'},'Location','northwest','FontSize',12);
grid on; grid minor;
subplot(212);
semilogy(ENR,final_Beta,'-.ko','LineWidth',1.5);hold on;
semilogy(ENR,opt_beta,'--ms','LineWidth',1.5);
xlabel('ENR [dB]'); ylabel('Optimal Beta');
grid on; grid minor;

function y = fconv(x,Lx,Ly,Ly2,H)
% Convolution in frequency domain using power of 2 fft

% input fft
X = fft(x, Ly2);

% multiply with precalculated freq domain signal
Y = X.*H;

% inverse fft and truncation
y = real(ifft(Y, Ly2));
y = y(1:1:Ly);

if mod(Lx,2) == 0
    y = y(Lx/2 + 1 : end - (Lx/2) + 1);
else
    y = y((Lx+1)/2 : end - ((Lx+1)/2) + 1);
end

end

function MSE = calc_perf_uniform_ppm(beta,Nrun,ENRlin)

dt = 1/(100*beta);
t = -1.125:dt:1.125;
tIdx = find(abs(t) <= 0.5);

ppmPulse = zeros(size(t));
ppmPulse(abs(t) < 1/(2*beta)) = sqrt(beta);
ppmPulse = ppmPulse / sqrt(sum(abs(ppmPulse.^2)*dt));

Lx = length(ppmPulse);
Ly = length(ppmPulse)+length(ppmPulse)-1;
Ly2 = pow2(nextpow2(Ly));
PPMfreq = fft(ppmPulse, Ly2);

currMSE = zeros(1,Nrun);
currENR = ENRlin;
parfor n=1:Nrun
    
    % generate source
    S = rand - 0.5;
    
    % PPM Modulation
    TxPulse = zeros(size(t));
    TxPulse(abs(t - S) < 1/(2*beta)) = sqrt(beta);
    TxPulse = sqrt(currENR)*TxPulse/sum(abs(TxPulse.^2)*dt);
    
    % AWGN
    noise = randn(size(t));
    noise = sqrt(1/(2*dt))*noise;% sqrt(2*Fs/W)*noise;
    r = TxPulse + noise;
    
    % PPM Correlator receiver
    PPMcorr = fconv(r,Lx,Ly,Ly2,PPMfreq);
    PPMcorr = PPMcorr(tIdx);
    
    [~,maxIdx] = max(PPMcorr);
    sHat = t(min(tIdx) + maxIdx - 1);
    if sHat > 0.5
        sHat = 0.5;
    elseif sHat <= -0.5
        sHat = -0.5;
    end
    
    currMSE(n) = (S - sHat)^2;
end

MSE = sum(currMSE)/Nrun;
end