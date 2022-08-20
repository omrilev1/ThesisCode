%% Simulation of Successive MLM and Linear transmissions, for the transmissions of an uniform source
%% over a channel with unlimited bandwidth. We optimize the scheme parameters according to the simulated profile

close all; clear all; clc;

profileOrder = 2;
saveResults = 1;

% Init parameters
SNR = -35:(0.5e-5):30;  % Direct channel SNR
snrLin = 10.^(SNR/10);

profile = 1./(1 + (snrLin).^profileOrder);

maxStages = 1e7; % number of levels

%% Scheme parameters:
c = 0.00137; d = 0.999; Delta = c; 

% SNR For intersection with profile 
Q_k = (Delta * (1:maxStages)); 
SNR_k = Q_k; 
A_vec = [Delta d.^(1:(maxStages - 1))*Delta];
beta_vec = 1 + (Q_k.^2) - cumsum(A_vec) .*  (Q_k);

E_uncoded = Delta*(d.^(0:maxStages)) ./ ...
    (1 + (Delta^2) * ((0:maxStages) + 1).^2 -...
    ((0:maxStages) + 1).*(Delta^2).*(1 - d.^((0:maxStages) + 1))/(1 - d));

E_dig = (1./(Delta*(1:maxStages))) .* ...
        log(1 + ...
        Delta^2 * (2*(1:maxStages) + 1 - (1:maxStages).*(d.^(1:maxStages)) -...
        (1 - d.^(1 + (1:maxStages)))/(1 - d))...
        ./(1 + Delta^2 * (1:maxStages).^2)...
        );

disp(strcat('Total Energy = ',num2str(sum(E_uncoded) + sum(E_dig))))


%% Simulation parameters

D_tot = zeros(length(SNR),1);
totalEnergy = ones(length(SNR),1);

currNumOfLevels = 1;
LastStepDist = 1; 
for i=1:length(SNR)
    
    prevNumOfLevels = currNumOfLevels;
    if i > 1
        if D_tot(i - 1) >= profile(i - 1)
            currNumOfLevels = prevNumOfLevels + 1;
        end
    end
    
    if prevNumOfLevels ~= currNumOfLevels
        LastStepDist = beta_vec(currNumOfLevels - 1); 
        A_tot = sum(A_vec(1:currNumOfLevels)); 
    end
    
    if currNumOfLevels == 1
        D_tot(i) = 1/(1 + 2*A_vec(1)*snrLin(i));
        totalEnergy(i) = E_uncoded(1);
    else
        D_tot(i) = 1/(LastStepDist + 2*A_tot*snrLin(i));
        totalEnergy(i) = sum(E_uncoded(1:currNumOfLevels)) +...
            sum(E_dig(1:(currNumOfLevels - 1)));
    end
    
    if mod(SNR(i),1) == 0
        disp(strcat('Finished SNR = ',num2str(SNR(i))));
    end
end

if saveResults
    save(strcat('Tuncel_InfDim_Performace',num2str(profileOrder)),'D_tot','totalEnergy','SNR');
end

figure;
semilogy(SNR,D_tot,'LineWidth',2.5); hold on;
semilogy(SNR,profile,'-.','LineWidth',2.5);
xlabel('ENR [dB]','FontSize',14); ylabel('Distortion [dB]','FontSize',14);
lgd = legend({'Baniasadi Tuncel','Profile'},'FontSize',16,'TextColor','Black');
grid on; grid minor;
title(strcat('Distortion and Profile, order = ',num2str(profileOrder)),'FontSize',14);


figure;
semilogy(SNR,totalEnergy,'LineWidth',2.5); hold on;
xlabel('ENR [dB]','FontSize',14); ylabel('Distortion [dB]','FontSize',14);
lgd = legend({'Baniasadi Tuncel'},'FontSize',16,'TextColor','Black');
grid on; grid minor;
title(strcat('Total Accumulated Energy, order = ',num2str(profileOrder)),'FontSize',14);

