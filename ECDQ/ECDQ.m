close all; clear all; clc;

%% This Script Calculate the R(D) Curve for ECDQ scheme - 1D case
% We calculate the R(D) For constant dither , and then average and plot
% for the random dither case

% Run Parameters :
% Signal Variance is 1 , but the distribution can be set as input

clear ; clc;

pdfType = {'Gaussian','Laplace','Exp'}; % 'Exp' , 'Gaussian' , 'Laplace'
nPoints = 5;
DeltaVec = [0.05:0.05:nPoints];
variance = 1;
ditherLength = 240;

zeroDitherH = zeros(length(pdfType),length(DeltaVec));
zeroDitherDist = zeros(length(pdfType),length(DeltaVec));

avgH = zeros(length(pdfType),length(DeltaVec));
avgDist = zeros(length(pdfType),length(DeltaVec));

currH = zeros(length(pdfType),ditherLength,length(DeltaVec));
currDist = zeros(length(pdfType),ditherLength,length(DeltaVec));

for i=1:length(pdfType)
    for j=1:length(DeltaVec)
        
        [pdf,cordX,cordY,dx] = pdfGenerator(pdfType(i),variance,1);
        Delta = DeltaVec(j);
        
        if min(cordX >= 0)
            
            dither = linspace(-1*Delta/2,Delta/2,ditherLength);
            effectiveLength = ceil(9*sqrt(variance) / Delta);
            
            % make number of Delta segments odd
            if mod(effectiveLength,2) == 0
                effectiveLength = effectiveLength + 1;
            end
            codebook  = 0  : Delta : effectiveLength*Delta;
            partition = Delta/2 : Delta : effectiveLength*Delta;
        else
            dither = linspace(-Delta/2,Delta/2,ditherLength);
            effectiveLength = ceil(7*sqrt(variance) / Delta);
            % make number of Delta segments odd
            if mod(effectiveLength,2) == 0
                effectiveLength = effectiveLength + 1;
            end
            
            codebook = -((effectiveLength+1)/2)*Delta : Delta : ((effectiveLength+1)/2)*Delta;
            partition = -((effectiveLength-1)/2)*Delta - Delta/2 : Delta : ((effectiveLength-1)/2)*Delta + Delta/2;
        end
        for k=1 : length(dither)
            % quantize distribution and calculate the probabilites : Pr{s+u in lambda}
            quantIndex = quantiz(cordX + dither(k),partition,codebook);
            segmentsVal = unique(quantIndex);
            
            currH(i,k,j) = 0;
            for valIdx = 1 : length(segmentsVal)
                quantSegment = (quantIndex == segmentsVal(valIdx));
                prob = sum(pdf(quantSegment))*dx;
                currDist(i,k,j) = currDist(i,k,j) + dx*sum(pdf(quantSegment) .* (codebook(valIdx) - (cordX(quantSegment) + dither(k))).^2);
                
                if prob > 0
                    currH(i,k,j) = currH(i,k,j) - prob*log2(prob);
                end
            end
        end
        
        entropyForAvg = currH(i,:,j);
        distForAvg = currDist(i,:,j);
        
        avgH(i,j) = sum(entropyForAvg)/ditherLength;
        avgDist(i,j) = sum(distForAvg)/ditherLength;
        
        [~,zeroDitherIdx] = min(abs(dither - 0));
        zeroDitherH(i,j) = entropyForAvg(zeroDitherIdx);
        zeroDitherDist(i,j) = distForAvg(zeroDitherIdx);
    end
end

deltaToPlot = [1:1:4];
pdfToPlot = {'Gaussian','Laplace','Exp'};
[indices] = resultPlot(deltaToPlot,DeltaVec,pdfType,pdfToPlot,currH,currDist,ditherLength,zeroDitherH,zeroDitherDist,variance);

function [indices] = resultPlot(deltaToPlot,deltaVec,pdfToPlot,pdfVec,H,D,ditherLength,zdH,zdD,variance)

% Calculate the Avg over the dither
avgH = reshape(mean(H,2),size(H,1),[]);
avgDist = reshape(mean(D,2),size(D,1),[]);

% find Delta's indices in the delta vector
indices = zeros(1,length(deltaToPlot));
for i=1:length(deltaToPlot)
    [~,indices(i)] = min(abs(deltaVec - deltaToPlot(i)));
end

% Figure titles
pdfIdx = zeros(1,length(pdfToPlot));
currTitle = [];
for i=1:length(pdfToPlot)
    for j=1:length(pdfVec)
        count = 1;
        if strcmp(pdfVec(j),pdfToPlot(i))
            pdfIdx(i) = count;
            break;
        else
            count = count + 1;
        end
    end
    currTitle = [currTitle pdfToPlot(i)];
end

% iterate over the distributions and plot the resulting curves
for i=1:length(pdfToPlot)
    
    if strcmp(pdfToPlot(i),'Gaussian')
        % R(D) Shannon Gaussian curve
        shannonD = 0.001:0.001:1;
        Rd = max(0.5*log(1./shannonD),0);
    else
        % R(D) Shannon SLB curve
        shannonD = 0.001:0.001:1;
        SourceEntropy = EntropyCalc(pdfToPlot(i),variance); 
        Rd = max(SourceEntropy - 0.5*log(2*pi*exp(1)*shannonD),0);
    end
    
    % convex hull of R(D|U)
    currH = reshape(H(i,:,:),1,[]);
    currD = reshape(D(i,:,:),1,[]);
    convHull_Idx = convhull([50 currD],[50 currH]);
    
    % take only the points on the relevant region
%     [~,relevant_ConvHull_Idx] = find(currD(convHull_Idx) < 4);
%     convHull_Idx_ToPlot = convHull_Idx(relevant_ConvHull_Idx);
    convHull_Idx_ToPlot = convHull_Idx(:);
    
    % Plot
    figure; hold all
    plot(shannonD,Rd,'-','LineWidth',2.5)
    plot(currD(convHull_Idx_ToPlot),currH(convHull_Idx_ToPlot),'--','LineWidth',2.5)
    plot(avgDist(i,:),avgH(i,:),'-','LineWidth',2.5)
    plot(zdD(i,:),zdH(i,:),'-.','LineWidth',2.5)
    currLegend = {'Shannon Lower Bound for R(D)','Optimum Dither Convex Hull','ECDQ','ECQ'};
    for j=1:length(deltaToPlot)
        plot(D(i,:,indices(j)),H(i,:,indices(j)),'--','LineWidth',2.5)
        currLegend = [currLegend strcat('\Delta =',num2str(deltaToPlot(j)))];
    end
    plot(1*ones(1,100),linspace(0,0.005,100),'-k*','LineWidth',1.5)
    
    grid on; grid minor;
    xlabel('D'); ylabel('R [bits]')
    legend(currLegend);
    title(currTitle(i));
    
    xlim([0 2]); ylim([0 3.5])
    
    % Plot the optimal Delta of the points on the convex hull
    [I,J] = ind2sub([size(H,2) size(H,3)],convHull_Idx_ToPlot);
    optimalDither = zeros(1,length(I));
    for k = 1 : length(I)
        currDither = linspace(-deltaVec(J(k))/2,deltaVec(J(k))/2,ditherLength);
        optimalDither(k) =  currDither(I(k))/deltaVec(J(k));
    end
    figure;
    plot(sqrt(deltaVec(J)),optimalDither,'gp','LineWidth',1.5);
    grid on; grid minor;
    xlabel('\surd{\Delta}'); ylabel('dither/\Delta');
    title(strcat('optimal Dither Value - ',currTitle(i)));
    
    
end
end





