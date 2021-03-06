close all; clear all; clc;

%% This Script Calculate the R(D) Curve for ECDQ scheme - 2D case
% We calculate the R(D) For constant dither , and then average and plot
% for the random dither case

clear ; close all; clc;

pdfType = {'Gaussian','Laplace','Exp'}; % 'Exp' , 'Gaussian' , 'Laplace' , {'Gaussian','Laplace','Exp'};
nPoints = 4;
DeltaVec = 0.25:0.25:nPoints;
variance = 0.25;
ditherLength = 2000; % 5e2;
effectiveDitherLength = ditherLength;

avgH = zeros(length(pdfType),length(DeltaVec));
avgDist = zeros(length(pdfType),length(DeltaVec));


zeroDitherH = zeros(length(pdfType),length(DeltaVec));
zeroDitherDist = zeros(length(pdfType),length(DeltaVec));

currH = zeros(length(pdfType),effectiveDitherLength,length(DeltaVec));
currDist = zeros(length(pdfType),effectiveDitherLength,length(DeltaVec));


% lattice generator matrix
G = [0 sqrt(3); 2 1];

% generate lattice with the regular matrix , will be normalized according
% to delta in the main loop
[basic_latticePoints] = latticeGen(G,50);

for i=1:length(pdfType)
    for j=1:length(DeltaVec)
        
        [pdf,cordX,cordY,dx] = pdfGenerator(pdfType(i),variance,2);
        Delta = DeltaVec(j);
        
        % generate 2D dither
        dither = generate2D_dither(Delta,ditherLength,2,Delta*basic_latticePoints);
        
        % generate 2D lattice points (codebook)
        latticePoints = Delta*basic_latticePoints;
        
        rowStack_x = cordX(:);
        rowStack_y = cordY(:);
        
        parfor k=1 : size(dither,2)
            % quantize distribution and calculate the probabilites : Pr{s+u in lambda}
            currDither = dither(:,k);
            
            tempSym = [rowStack_x.' ; rowStack_y.'] + currDither;
            latticePoints_temp = reshape(latticePoints,2,1,[]);
            error = tempSym - latticePoints_temp ;
            errorAbs = (error(1,:,:)).^2 + (error(2,:,:)).^2;
            [~,quantIndex] = min(errorAbs,[],3);
            
            segmentsVal = unique(quantIndex);
            
            currH(i,k,j) = 0;
            for valIdx = 1 : length(segmentsVal)
                quantSegment = (quantIndex == segmentsVal(valIdx));
                
                indx_2D = reshape(quantSegment,size(cordX));
                prob = sum(sum(pdf(indx_2D)))*dx;
                
                % calculate distortion efficiently in 2D
                valid_XY = [rowStack_x(quantSegment) rowStack_y(quantSegment)].';
                XY_plusDither = valid_XY + dither(:,k);
                total_Error = latticePoints_temp(:,1,segmentsVal(valIdx)) - XY_plusDither ;
                
                %                 errorNorm = sqrt((total_Error(1,:)).^2 + (total_Error(2,:)).^2);
                errorNorm = ((total_Error(1,:)).^2 + (total_Error(2,:)).^2)/2;
                
                total_Error2D = zeros(size(pdf));
                total_Error2D(indx_2D) = errorNorm;
                
                currDist(i,k,j) = currDist(i,k,j) + dx*sum(sum(pdf .* total_Error2D));
                
                if prob > 0
                    currH(i,k,j) = currH(i,k,j) - prob*log2(prob);
                end
            end
        end
        
        entropyForAvg = currH(i,:,j);
        distForAvg = currDist(i,:,j);
        
        avgH(i,j) = sum(entropyForAvg)/effectiveDitherLength;
        avgDist(i,j) = sum(distForAvg)/effectiveDitherLength;
        
        [~,zeroDitherIdx] = min(sum(abs(dither),1));
        zeroDitherH(i,j) = entropyForAvg(zeroDitherIdx);
        zeroDitherDist(i,j) = distForAvg(zeroDitherIdx);
        
        display(strcat('finished Delta = ',num2str(DeltaVec(j)),' For PDF',pdfType(i)))
    end
end

deltaToPlot = 1:1:4;
pdfToPlot = {'Gaussian','Laplace','Exp'};
[indices] = resultPlot(deltaToPlot,DeltaVec,pdfType,pdfToPlot,currH,currDist,ditherLength,zeroDitherH,zeroDitherDist,variance);

function [indices] = resultPlot(deltaToPlot,deltaVec,pdfToPlot,pdfVec,H,D,ditherLength,zeroDitherH,zeroDitherDist,variance)

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
    
    if strcmp(pdfToPlot,'Gaussian')
        % R(D) Shannon curve
        shannonD = 0.001:0.001:1;
        Rd = max(0.5*log(variance./shannonD),0);
    else
        % R(D) Shannon curve
        shannonD = 0.001:0.001:1;
        SourceEntropy = EntropyCalc(pdfToPlot(i),variance); 
        Rd = max(SourceEntropy - 0.5*log(2*pi*exp(1)*shannonD),0);
%         lambda = 1/sqrt(0.25);
%         Rd = max(-0.5*log(lambda*shannonD),0);
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
    plot(shannonD,max(Rd,0),'-','LineWidth',2)
    %     plot(currD(convHull_Idx_ToPlot),currH(convHull_Idx_ToPlot)/2,'cp')
    plot(avgDist(i,:),avgH(i,:)/2,'-','LineWidth',2)
    plot(zeroDitherDist(i,:),zeroDitherH(i,:)/2,'--','LineWidth',2)
    %     currLegend = {'Shannon : R(D) = 0.5*log(1/D)','Optimum Dither Convex Hull','ECDQ'};
    currLegend = {'Shannon : R(D) = 0.5*log(1/D)','ECDQ','Zero Dither'};
    
    curveStyle = ['ro';'ks';'cp';'mo';'gs'];
    for j=1:length(deltaToPlot)
        plot(D(i,:,indices(j)),H(i,:,indices(j))/2,curveStyle(j,:),'LineWidth',1.7)
        currLegend = [currLegend strcat('\Delta =',num2str(deltaToPlot(j)))];
    end
    grid on; grid minor;
    xlabel('D'); ylabel('R [bits]')
    legend(currLegend);
    title(currTitle(i));
     
end
end



