%% Analog PPM Performance simulation
close all; clear all; clc;

ENR = (9:1:20);
ENRlin = 10.^(ENR/10);

ENR_opt = 15;
ENR_opt_lin = 10^(ENR_opt/10); 

beta_opt = (13/8)^(1/3) * (ENR_opt_lin)^(-5/6) * exp(ENR_opt_lin/6);
Nrun = 2e5;
final_MSE = zeros(1,length(ENR));
final_Beta = zeros(1,length(ENR));

for i=1:length(ENR)
    
    currENRlin = ENRlin(i);
    final_MSE(i) = simulateGaussianPPM(beta_opt,currENRlin,Nrun,6.4);
        
    disp(strcat('Finished ENR = ',num2str(ENR(i))));
        
end

ENRbound = min(ENR):0.125:max(ENR); 
ENRboundLin = 10.^(ENRbound/10); 

% exact upper bound

D_S = ((13/8) + sqrt(2/beta_opt) * (sqrt(2*beta_opt.*ENRboundLin) - 1).*exp(-ENRboundLin .* (1 - 1./sqrt(2*beta_opt.*ENRboundLin)).^2)) ./...
    ((sqrt(beta_opt.*ENRboundLin) - 1/sqrt(2)).^4) + exp(-beta_opt.*ENRboundLin)./(beta_opt.^2);
D_L = 2*beta_opt.*sqrt(ENRboundLin).*exp(-ENRboundLin/2) .* (1 + 3*sqrt(2*pi./ENRboundLin) + 12*exp(-1)./(beta_opt.*sqrt(ENRboundLin)) ...
    + 8*exp(-1)./(sqrt(8*pi)*beta_opt) + sqrt(8./(pi*ENRboundLin)) + 12^(3/2) * exp(-3/2) ./(beta_opt.*sqrt(32*pi*ENRboundLin)));
optDist_Analytic = (D_S + D_L);

% Asymptotic bound 
D_S_asymp = (13/8)./(ENRboundLin * beta_opt).^2;
D_L_asymp = 2*beta_opt.*sqrt(ENRboundLin).*exp(-ENRboundLin/2);
optDist_Analytic_asymp = (D_S_asymp + D_L_asymp);
 
figure;
hold all;
plot(ENRbound,10*log10(optDist_Analytic),':','LineWidth',1.75);
plot(ENRbound,10*log10(optDist_Analytic_asymp),'--','LineWidth',1.75);
plot(ENR,10*log10(final_MSE),'-.ko','LineWidth',1.5);
xlabel('ENR [dB]'); ylabel('Distortion [dB]'); legend({'Upper Bound','Asymptotic Bound','Empirical'},'Location','northwest','FontSize',12);
grid on; grid minor;

function y = fconv(x,Lx,Ly,Ly2,H)
% Convolution in frequency domain using power of 2 fft
% Since the input signal is real we use known fft identities to accelerate
% the fft

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

function [MSE] = simulateGaussianPPM(beta,SNRlin,Nrun,overload)

dt = 1/(350*beta);

t = -overload:dt:overload;
t = t(:);
% t = parallel.pool.Constant(t); 

ppmPulse = zeros(length(t),1,'single');
ppmPulse(abs(t) < 1/(2*beta)) = sqrt(beta);
ppmPulse = ppmPulse / sqrt(sum(abs(ppmPulse.^2)*dt));
% ppmPulse(abs(t) > 0.5) = 0;

Lx = length(ppmPulse);
Ly = length(ppmPulse)+length(ppmPulse)-1;
Ly2 = pow2(nextpow2(Ly));
PPMfreq = fft(ppmPulse, Ly2);
currMSE_MAP = zeros(1,Nrun);

parfor n=1:Nrun
    
    % generate source - Gaussian Source with edge truncation
    S = randn;
    if abs(S) > overload
        S = overload*sign(S);
    end
    
    TxPulse = zeros(length(t),1,'single');
    TxPulse(abs(t - S) < 1/(2*beta)) = sqrt(beta);
    TxPulse = sqrt(SNRlin)*TxPulse/sum(abs(TxPulse.^2)*dt);
    
    noise = randn(length(t),1,'single');
    noise = sqrt(1/(2*dt))*noise;
    r = TxPulse + noise;
    
    % PPM Correlator receiver
    PPMcorr = fconv(r,Lx,Ly,Ly2,PPMfreq);
    
    [~,maxIdx] = max(sqrt(SNRlin)*PPMcorr*dt - 0.5*(t.^2));
    sHat_MAP = t(maxIdx);
    
    currMSE_MAP(n) = (S - sHat_MAP)^2;
    
end

MSE = sum(currMSE_MAP)/Nrun;


end

