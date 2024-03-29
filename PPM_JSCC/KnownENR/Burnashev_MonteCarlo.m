% Digital PPM: Simulation of Burnashev scheme
% The simulation is for uniform source
clear all; clc; % close all; 
saveResults = 1; 

ENR = (4:1:16);
ENRlin = 10.^(ENR/10);

Nrun = 6e6; % 8e6 for Uniform source, 1e6 for Gaussian source 

delta = 1/2 + 1/(2*sqrt(2));
dist = {'Uniform'}; {'Uniform','Gaussian'};

final_MSE = zeros(length(dist),length(ENR));
final_Beta = zeros(length(dist),length(ENR));

disp('Burnashev Sim');
for distIdx = 1:length(dist)
    disp('===============================');
    disp(strcat('Current Distribution: ',dist{distIdx}));
    for i=1:length(ENR)
        % pick beta to find the optimal point of the curve
        
        switch dist{distIdx}
            case 'Uniform'
                beta_opt = 1 * exp(ENRlin(i)/6);
                beta = beta_opt; % max(1,beta_opt*(0.45:(1/32):3));
            case 'Gaussian'
                beta_opt = 1.3719 * exp(ENRlin(i)/6);
                beta = beta_opt; % max(1,beta_opt*(0.45:(1/32):3));
        end
        MSE = zeros(size(beta));
        currENRlin = ENRlin(i);
        for beta_idx = 1:length(beta)
            curr_beta = beta(beta_idx);
            currMSE = calc_perf_Burnashev(curr_beta,Nrun,currENRlin,dist{distIdx});
            MSE(beta_idx) = currMSE;
        end
        [optVal,optIdx] = min(MSE(1:beta_idx));
        final_MSE(distIdx,i) = MSE(optIdx);
        final_Beta(distIdx,i) = beta(optIdx);
        disp(strcat('Finished ENR = ',num2str(ENR(i))));
    end
    
    switch dist{distIdx}
        case 'Uniform'
            optSDR_Analytic_Exact = ((delta^(2/3))/4) .* exp(-ENRlin/3);
            opt_beta = exp(ENRlin/6);
            if saveResults
                save('UniformBurnashevOptSDR.mat','ENR','final_MSE','optSDR_Analytic_Exact','final_Beta','opt_beta');
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
            opt_beta = exp(ENRlin/6);
            if saveResults
                save('GaussianBurnashevOptSDR.mat','ENR','final_MSE','optSDR_Analytic_Exact','final_Beta','opt_beta');
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


function MSE = calc_perf_Burnashev(beta,Nrun,ENRlin,distType)

beta_round = ceil(beta);
orthMat = sparse(eye(beta_round));

currMSE = zeros(1,Nrun);
currENR = ENRlin;

%calculate optimal Quantization for the Gaussian case
if strcmp(distType,'Gaussian')
    dx = 1e-4; x = 0:dx:8;
    pdf = (1/sqrt(2*pi))*exp(-x.^2 / 2) * dx;
    
    if mod(beta_round,2) == 0
        indices = (1:beta_round/2).';
        Delta = (1/beta_round):1e-1/beta_round:16/beta_round;
        quant_vals = indices * Delta;
        
    else
        indices = (0:(beta_round - 1)/2).';
        Delta = (1/beta_round):(0.5e-1)/beta_round:16/beta_round;
        quant_vals = indices * Delta;
    end
    
    % calculate the quantization error, to find the optimal delta
    quant_MSE = zeros(1,size(Delta,2));
    if size(quant_vals,1) == 1
        quant_MSE = sum((quant_vals(:) - x).^2 .* repmat(pdf,size(quant_vals,1),1),2);
    else
        parfor i=1:size(Delta,2)
            curr_vals =  quant_vals(:,i);
            [~,quants] = quantiz(x,curr_vals(1:end-1) + Delta(i)/2,curr_vals);
            quant_MSE(i) = sum((quants(:) - x(:)).^2 .* pdf(:));
        end
    end
    [~,optDelta_Idx] = min(quant_MSE);
    optDelta = Delta(optDelta_Idx);
    
    if(indices(1) == 0)
        opt_codebook = [-flipud(indices(2:end)) * optDelta;indices * optDelta];
        opt_partition = [flipud(-1*(indices(1:end-1) * optDelta + optDelta/2)); (indices(1:end-1) * optDelta + optDelta/2)];
    else
        opt_codebook =  [-flipud(indices) * optDelta;indices * optDelta];
        opt_partition = [flipud(-1*(indices(1:end-1) * optDelta + optDelta/2)); 0; (indices(1:end-1) * optDelta + optDelta/2)];
    end
    
else
    opt_partition = [];
    opt_codebook = [];
end

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
            S_pos = S + 0.5; % [0,1]

            % Improved quantization - no half bin offset 
            Q_levels = 1/(2*beta_round) : (1/beta_round) : 1-1/(2*beta_round);
            [~, S_q] = min(abs(S_pos-Q_levels));
            S_q = S_q - 1;
%             S_q = round((beta_round - 1)* S_pos);
        case 'Gaussian'
            S = randn;
            S_q = quantiz(S,opt_partition,opt_codebook);
    end
    
    % Modulation
    TxPulse = sparse(orthMat(S_q + 1,:));
    TxPulse = sqrt(currENR)*TxPulse;
    
    % AWGN
    noise = randn(size(TxPulse));
    noise = sqrt(1/2)*noise;% sqrt(2*Fs/W)*noise;
    r = TxPulse + noise;
    
    % Correlator receiver: multiplication with eye matrix 
    % is equivalent to pick the maximal element of the vector 
%     corr = orthMat*r(:);
    
    [~,maxIdx] = max(r(:));
    
    % De-quantization
    switch distType
        case 'Uniform'
            sHat = (maxIdx - 1 + 0.5)/beta_round - 0.5;
        case 'Gaussian'
            sHat = opt_codebook(maxIdx);
    end
    
    currMSE(n) = (S - sHat)^2;
end

MSE = sum(currMSE)/Nrun;
end