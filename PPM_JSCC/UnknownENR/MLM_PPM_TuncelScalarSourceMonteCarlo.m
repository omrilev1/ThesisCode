%% Simulation of Successive MLM and Linear transmissions, for the transmissions of an uniform source
%% over a channel with unlimited bandwidth. We optimize the scheme parameters according to the simulated profile

close all; clear all; clc;

profileOrder = 2;
saveResults = 0;

% Init parameters and arrays structures
SNR = -7:0.25:35;  % 8              % Direct channel SNR
snrLin = 10.^(SNR/10);

profile = 1./(1 + (snrLin).^profileOrder);

maxStages = 100; % number of Digital levels
ModOrder = [1;zeros(maxStages - 1,1)]; % first layer

%% Modulo parameters:
% interval
d1 = 0.999; alpha1 = 1.5;
d2 = 0.999; alpha2= 1.5;
energyAlloc_Analog = (d1.^(1:maxStages)) * alpha1;
energyAlloc_Analog = energyAlloc_Analog.';

energyAlloc_Digital = (d2.^(1:maxStages)) * alpha2; % temporal values
energyAlloc_Digital = energyAlloc_Digital.';
TargetPe_Digital = 1e-2; 
%% Simulation parameters
N_avg = 2^13;
E_0 = 1;

D_Tot = zeros(length(SNR),1);
totalEnergy = zeros(length(SNR),1);

currNumOfLevels = 1;finished = 0;

for i=1:length(SNR)

    prevNumOfLevels = currNumOfLevels;
    if i > 1
        if D_Tot(i - 1) > profile(i - 1)
            currNumOfLevels = prevNumOfLevels + 1;
        end
    end

    if prevNumOfLevels ~= currNumOfLevels
        % calculate optimal number of quantization levels and optimize the
        % digital energy for the next level
        CurrEnergy_Digital = energyAlloc_Digital(currNumOfLevels - 1) * 10^(SNR(i)/10);
        TotalAnalogSNR = sum(energyAlloc_Analog(1:currNumOfLevels)) * 10^(SNR(i)/10);
        Curr_M = optimizeParams(TotalAnalogSNR,CurrEnergy_Digital,TargetPe_Digital);
        if isnan(Curr_M)
            disp('You need to increase energy');
            break
        end
        ModOrder(currNumOfLevels) = Curr_M;

        % calculate total analog vector
        sigma_x_2 = (1/(12*(ModOrder(currNumOfLevels))^2));
        aVec = sqrt(12*energyAlloc_Analog(1:currNumOfLevels)).*(cumprod(ModOrder(1:currNumOfLevels)));
        C_xy = sigma_x_2 * aVec;
        C_yy = sigma_x_2 * (aVec*aVec') + (1/snrLin(i)) * eye(currNumOfLevels);
        AnalogEstimVecor = C_xy' / C_yy;
    end

    currDist = zeros(1,N_avg);
    currE_tot = zeros(1,N_avg);
    for n = 1 : N_avg

        % Randomize a source, and build the sources for each transmission 
        % layers
        S = rand - 0.5;
        TxSource = zeros(1,maxStages);
        TxSource_Analog = zeros(1,maxStages);
        TxSource_Digital = zeros(1,maxStages);
        dither = zeros(1,maxStages);

        TxSource_Analog(1) = S;
        TxSource(1) = S;

        if currNumOfLevels > 1
            for idx=2:currNumOfLevels
                Q_levels = 1/(2*ModOrder(idx)) : (1/ModOrder(idx)) : 1-1/(2*ModOrder(idx));
                [~, S_q] = min(abs(TxSource(idx - 1) + 0.5 - Q_levels));
                S_q = S_q - 1;
                TxSource_Digital(idx) = S_q;
                TxSource_Analog(idx) = TxSource(idx-1) - ((S_q + 0.5)/ModOrder(idx) - 0.5); % \in [-1/(2*ModOrder),1/(2*ModOrder)]
                TxSource(idx) = ModOrder(idx) * TxSource_Analog(idx); % \in [-1/2,1/2]
            end
        end

        % simulate the successive transmission and reception
        y_Analog = zeros(currNumOfLevels,1);
        y_Digital = zeros(currNumOfLevels,1);
        noiseVec = zeros(currNumOfLevels,1);
        sHat = 0;
        for k=1:currNumOfLevels
            %% channel part
            if k==1
                % First layer: simple linear modulation
                noiseVec(1) = sqrt(1/snrLin(i))*randn;
                y_Analog(1) = sqrt(12) * sqrt(energyAlloc_Analog(1)) * S + noiseVec(1);
                y_Digital(1) = 0;
            else
                % Next layers: orthogonal signaling with M signals

                % Improved quantization - no half bin offset
                orthMat = eye(ModOrder(k));
                TxPulse = orthMat(TxSource_Digital(k) + 1,:);
                TxPulse = sqrt(energyAlloc_Digital(k - 1)*snrLin(i))*TxPulse;

                % AWGN
                noise = randn(size(TxPulse));
                noise = sqrt(1/2)*noise;
                r = TxPulse + noise;

                % Correlator receiver: multiplication with eye matrix
                % is equivalent to pick the maximal element of the vector
                [~,maxIdx] = max(r(:));

                % De-quantization
                noiseVec(k) = sqrt(1/snrLin(i))*randn;
                y_Digital(k) = ((maxIdx - 1 + 0.5)/ModOrder(k) - 0.5);
                y_Analog(k) = sqrt(12) * sqrt(energyAlloc_Analog(k)) * TxSource(k) + noiseVec(k);
            end

            %% Decoder: combine digital and analog parts, to produce the
            %% overall estimate
            if k==1
                sHat = y_Analog(1) / (sqrt(12) * sqrt(energyAlloc_Analog(1)));
                e1 = abs(sHat - S);
            else
                % calculate total estimate from digital parts
                digital_estimate = y_Digital./cumprod([1;ModOrder(1:(currNumOfLevels-1))]);

                % generate vector of analog measurements
                analog_meas = y_Analog - sqrt(12*energyAlloc_Analog(1:currNumOfLevels)).*flipud(cumsum(y_Digital));
                sHat = sum(digital_estimate) + AnalogEstimVecor * analog_meas;
                e2 = abs(sHat - S);
                if e2 > e1
                    disp(strcat('Fuck My Life, S = ',num2str(S),' n = ',num2str(n),' digital = ',num2str(y_Digital(2)), ' sq = ',num2str(S_q)));
                    disp(strcat('Theory analog measurement1: ',num2str(sqrt(12*energyAlloc_Analog(1))*TxSource_Analog(2) + noiseVec(1))));
                    disp(strcat('Theory analog measurement2: ',num2str(sqrt(12*energyAlloc_Analog(2))*TxSource(2) + noiseVec(2))));
                    disp(strcat('Actual measurement1: ',num2str(analog_meas(1))));
                    disp(strcat('Actual measurement2: ',num2str(analog_meas(2))));
                end
            end
        end
        currDist(n) = (sHat - S)^2;
    end

    D_Tot(i) = sum(currDist)/N_avg;
    totalEnergy(i) = sum(energyAlloc_Analog(1:currNumOfLevels)) + ...
        sum(energyAlloc_Digital(1:currNumOfLevels));

    if mod(i,5) == 0
        disp(strcat('Finished SNR = ',num2str(SNR(i))));
    end

end

if saveResults
    save(strcat('BaniasadiTuncel_Scalar_Profile',num2str(profileOrder)),'D_Tot','totalEnergy','Delta','alpha','SNR','profile');
end

figure;
semilogy(SNR,D_Tot,'LineWidth',2.5); hold on;
semilogy(SNR,profile,'-.','LineWidth',2.5);
xlabel('ENR [dB]','FontSize',14); ylabel('Distortion [dB]','FontSize',14);
lgd = legend({'Baniasadi & Tuncel','Profile'},'FontSize',16,'TextColor','Black');
grid on; grid minor;
title(strcat('Distortion and Profile, order = ',num2str(profileOrder)),'FontSize',14);


function M_Opt = optimizeParams(TotalAnalogSNR,CurrEnergy,TargetPe)
delta = 1/2 + 1/(2*sqrt(2));
if delta * exp(-CurrEnergy/4) < TargetPe
    M_Opt = round(TargetPe*exp(CurrEnergy/2)/delta);
else
    M_Opt = round(exp(((sqrt(2*CurrEnergy) - sqrt(2*log(delta/TargetPe)))^2)/2));
end
% M = 2:1:floor(exp(CurrEnergy));
% 
% % calculate probability of error - Gallager's bound
% P_e_bound = zeros(size(M));
% delta = 1/2 + 1/(2*sqrt(2));
% P_e_bound(M <= exp(CurrEnergy/4)) = delta*exp(log(M(M <= exp(CurrEnergy/4))) - CurrEnergy/2);
% P_e_bound(M > exp(CurrEnergy/4)) = delta*exp(-0.5 * (sqrt(2*CurrEnergy) - sqrt(2*log(M(M > exp(CurrEnergy/4))))).^2);
% P_e_bound = min(1,P_e_bound);

% calculate distortion
% dist = (1 - P_e_bound).*(1./(12*(M.^2))) .* 1./(1 + TotalAnalogSNR * (1./(12*(M.^2)))) + P_e_bound/6;
% 
% [~,M_Opt_idx] = min(dist);
% M_Opt = M(M_Opt_idx);

% [~,M_Opt_idx] = min(abs(TargetPe - P_e_bound));
% if isempty(M_Opt_idx)
%     M_Opt = NaN;
% else
%     M_Opt = M(M_Opt_idx);
% end

end



