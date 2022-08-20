% Digital PPM: Simulation of Burnashev scheme
% The simulation is for uniform source
clear all; clc; % close all; 
saveResults = 1; 

ENR = (4:1:16);
ENRlin = 10.^(ENR/10);

Nrun = 2e5; 

delta = 1/2 + 1/(2*sqrt(2));
dist = {'Uniform'}; % {'Uniform','Gaussian'}
sigmaOptimal = 1; 

final_MSE = zeros(length(dist),length(ENR));
final_Beta = zeros(length(dist),length(ENR));
final_Mean = zeros(length(dist),length(ENR));
disp('Sevinc Tuncel Sim');
for distIdx = 1:length(dist)
    disp('===============================');
    disp(strcat('Current Distribution: ',dist{distIdx}));
    for i=1:length(ENR)
        % pick beta to find the optimal point of the curve
        
        if ENR(i) <= 13
            Nrun = 1e5; % 6e5 for Gaussian source
        else
            Nrun = 3e6; % 6e6 for uniform source 
        end
        
        %calculate optimal Quantization for the Gaussian case
        if strcmp(dist(distIdx),'Gaussian')
            % equations from Sevinc' Tuncel paper, appendix A, for the calculation
            % of the compander
            dx = 1e-5; x = -6.4:dx:6.4;
            if sigmaOptimal
                sigma_2 = 1 + sqrt(2); 
                c = 1.2794; 
                lambda_x = ((1/sqrt(2*pi*sigma_2)) * exp(-x.^2 / (2*sigma_2)));
            else 
                c = 1.3719; 
                beta = 1.1764; 
                lambda_x = (1/(6^(1/3) * c^(2/3))) * (((1/sqrt(2*pi)) * exp(-x.^2 / 2)).^(1/3))  ./ ((delta*c*x.^2 + beta).^(1/3));
            end
            G = cumsum(lambda_x*dx);
            
        else
            dx = 1e-5; x = -1/2:dx:1/2;
            c = 1;
            lambda_x = (1/(6^(1/3) * c^(2/3))) * 1 ./ ((delta*c*x.^2 + 0.10925).^(1/3));
            G = cumsum(lambda_x*dx);
        end
                    
        beta_opt = c*exp(ENRlin(i)/6);% 1.6*exp(ENRlin(i)/12);
        beta = beta_opt; % max(1,beta_opt*(0.75:(1/4):1.55));
        MSE = zeros(size(beta));Mean = zeros(size(beta));
        currENRlin = ENRlin(i);
        
        for beta_idx = 1:length(beta)
            curr_beta = beta(beta_idx);
            [currMSE,currMean] = calc_perf_Tuncel(curr_beta,Nrun,currENRlin,dist{distIdx},G(:),x(:),dx);
            MSE(beta_idx) = currMSE;
            Mean(beta_idx) = currMean;
        end
        [optVal,optIdx] = min(MSE(1:beta_idx));
        final_MSE(distIdx,i) = MSE(optIdx);
        final_Beta(distIdx,i) = beta(optIdx);
        final_Mean(distIdx,i) = Mean(optIdx);
        disp(strcat('Finished ENR = ',num2str(ENR(i))));
    end
    
    switch dist{distIdx}
        case 'Uniform'
            optSDR_Analytic_Exact = ((delta^(2/3))/4) .* exp(-ENRlin/3);
            opt_beta = c * exp(ENRlin/6);
            if saveResults
                save('UniformTuncelOptSDR.mat','ENR','final_MSE','optSDR_Analytic_Exact','final_Beta','opt_beta');
            end 
            figure;subplot(211);
            hold all;
            plot(ENR,10*log10(1/12./final_MSE(distIdx,:)),'-.ko','LineWidth',1.5);
            plot(ENR,10*log10(1/12./optSDR_Analytic_Exact),'--ms','LineWidth',1.5);
            xlabel('ENR [dB]'); ylabel('SDR [dB]'); legend({'Empirical','Approximation'},'Location','northwest','FontSize',12);
            grid on; grid minor;
            subplot(212);
            semilogy(ENR,final_Beta(distIdx,:),'-.ko','LineWidth',1.5);hold on;
            semilogy(ENR,opt_beta,'--ms','LineWidth',1.5);
            xlabel('ENR [dB]'); ylabel('Optimal Beta');
            grid on; grid minor;
            
        case 'Gaussian'
            optSDR_Analytic_Exact = exp(-ENRlin/3);
            opt_beta = c * exp(ENRlin/6);
            if saveResults
                if sigmaOptimal
                    save('GaussianTuncelOptSDR_SigmaOptimalCompander.mat','ENR','final_MSE','optSDR_Analytic_Exact','final_Beta','opt_beta');
                else
                    save('GaussianTuncelOptSDR.mat','ENR','final_MSE','optSDR_Analytic_Exact','final_Beta','opt_beta');
                end
            end
            figure;subplot(211);
            hold all;
            plot(ENR,10*log10(1./final_MSE(distIdx,:)),'-.ko','LineWidth',1.5);
            plot(ENR,10*log10(1./optSDR_Analytic_Exact),'--ms','LineWidth',1.5);
            xlabel('ENR [dB]'); ylabel('SDR [dB]'); legend({'Empirical','Approximation'},'Location','northwest','FontSize',12);
            grid on; grid minor;
            subplot(212);
            semilogy(ENR,final_Beta(distIdx,:),'-.ko','LineWidth',1.5);hold on;
            semilogy(ENR,opt_beta,'--ms','LineWidth',1.5);
            xlabel('ENR [dB]'); ylabel('Optimal Beta');
            grid on; grid minor;
            
    end
end


function [MSE,Mean] = calc_perf_Tuncel(beta,Nrun,ENRlin,distType,G,x,dx)

beta_round = ceil(beta);
orthMat = sparse(eye(beta_round));

currMSE = zeros(1,Nrun);
currMean = zeros(1,Nrun);
currENR = ENRlin;

if (beta_round == 1)
    switch distType
        case 'Uniform'
            MSE = 1/12;
        case 'Gaussian'
            MSE = 1;
    end
    return
end
% Generate Quantization Table

parfor n=1:Nrun
    
    % generate source and quantize
    switch distType
        case 'Uniform'
            S = rand - 0.5;  % [-0.5,0.5]
            
        case 'Gaussian'
            S = randn;
            
    end
    
    % companding: we take the ivalue from G, and then
    % linearly interpolate with neighbours to get finer resolution
    S_Companding = G(round((S - min(x))/dx) + 1);
    
    % Scalar uniform quantization
    % Improved quantization - no half bin offset
    Q_levels = 1/(2*beta_round) : (1/beta_round) : 1-1/(2*beta_round);
    [~, S_q] = min(abs(S_Companding-Q_levels));
    S_q = S_q - 1;
%     S_q = round((beta_round - 1)* S_Companding);
    
    % Modulation
    TxPulse = sparse(orthMat(S_q + 1,:));
    TxPulse = sqrt(currENR)*TxPulse;
    
    % AWGN
    noise = randn(size(TxPulse));
    noise = sqrt(1/2) * noise;% sqrt(2*Fs/W)*noise;
    r = TxPulse + noise;
    
    % Correlator receiver: multiplication with eye matrix 
    % is equivalent to pick the maximal element of the vector     
    [~,maxIdx] = max(r(:));
    
    % De-quantization + De-companding
    iHat = (maxIdx - 1 + 0.5)/beta_round;
    [~,idx] = min(abs(G - iHat));
    sHat  = x(idx);
    
    currMSE(n) = (S - sHat)^2;
    currMean(n) = (S - sHat);
end

MSE = sum(currMSE)/Nrun;
Mean = sum(currMean)/Nrun;
end