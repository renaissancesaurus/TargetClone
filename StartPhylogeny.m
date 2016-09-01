
%% run per sample

%don't know how to read the file.. It either replaces everything with
%strings and str2num does not work because of the NAN
%without that there is no result

fnam='all3209.txt';
data=read_mixed_csv(fnam,'\t');
data = cellfun(@(s) {str2double(s)},data);
% data(:,17) = [];
% data(:,10) = [];
%manual input of data
%data = all6107;

%load the segmentation
fileID = fopen('Tmean_3209_seg.txt');
C = textscan(fileID,'%s');
fclose(fileID);

segmentation = C{1};


cloneNum = 2;

allMu = permn(0:1:100, cloneNum);
%allMu(1:11,:) = [];
validMu = arrayfun(@(x)isequal (x, 100), sum(allMu,2));
validMuInd = find(validMu == 1);
% validMuInd = validMuInd(1:11);
%allMu = [allMu zeros(length(allMu),1)];

kmin = 0;
kmax = 6;

[alleleMap, allAlleleDistanceMatrix] = computeAllAlleleDistances(kmin, kmax);
allPermutations = permn([kmin+1:kmax], cloneNum-1);

allSolutions = cell(1, size(data,2));
allSolutionsMu = cell(1,size(data,2));
allSolutionsPhasing = cell(1,size(data,2));

uniqueSegments = unique(segmentation);

%given the 
bestCSamples = cell.empty;
bestCSamplesScores = cell.empty;
%size(data,2)

data = {1/2, 1/2, 1/5, 0};
data = {1/2, 1/2, 0.4444, 0};

segmentation = {'1A'};
uniqueSegments = {'1A'};
for sample = 1:size(data,2)
   
    sampleData = cell2mat(data(:, sample));

    nanInd = isnan(sampleData);
    sampleData(nanInd) = [];
    newSegmentation = segmentation;
    newSegmentation(nanInd) = [];

    segmentedLAF = ones(1, length(uniqueSegments));
    for i = 1:length(uniqueSegments)

        %process per segment individually (strange results)

        segmentLAFInd = strcmp(newSegmentation, uniqueSegments(i));
        segmentedLAF(i) = median(sampleData(segmentLAFInd));

    end
    segmentedLAF = round((segmentedLAF*1000)) / 1000.0;
    %LAF = round((sampleData*1000)) / 1000.0;
    
    %run the model on these LAF

    allBest = cell(1, length(validMuInd));
    allBestScores = zeros(1, length(validMuInd));
    allBestPhasings = {};
    all3BestScores = cell.empty;
    parfor muIndex = 1:length(validMuInd)
        positionBest = cell(1,3);
        bestPhasing = cell(1,3);
        positionBestScores = zeros(length(uniqueSegments), 3);
        for segmentInd = 1:length(uniqueSegments)
            
            if isnan(segmentedLAF(segmentInd))

                positionBestScores(segmentInd,:) = 1; %do not include in weighing
                positionBest{1} = [positionBest{1}; 3];
                positionBest{2} = [positionBest{2}; 2];
                positionBest{3} = [positionBest{3}; 1];
                bestPhasing{1} = [bestPhasing{1} 'ABB'];
                bestPhasing{2} = [bestPhasing{2} 'AB'];
                bestPhasing{3} = [bestPhasing{3} 'B'];
                continue;
            end
          
            allScores = [];
            allPhasings = {};
            for permutation = 1:length(allPermutations)
                c = allPermutations(permutation,:);
                [score, bestLAFPhasing, likelyPhylogeny] = fitness(c, allMu(validMuInd(muIndex),:)/100, segmentedLAF(segmentInd), alleleMap, allAlleleDistanceMatrix);
                allScores = [allScores; score];
                allPhasings{permutation} = bestLAFPhasing;
                %iterNum = iterNum + 1;
            end
            [sortedScores, ind] = sort(allScores(:,1),1, 'descend');
            positionBestScores(segmentInd,:) = sortedScores(1:3,1);
            positionBest{1} = [positionBest{1}; allPermutations(ind(1),:)];
            positionBest{2} = [positionBest{2}; allPermutations(ind(2),:)];
            positionBest{3} = [positionBest{3}; allPermutations(ind(3),:)];
            bestPhasing{1} = [bestPhasing{1} allPhasings{ind(1)}];
            bestPhasing{2} = [bestPhasing{2} allPhasings{ind(2)}];
            bestPhasing{3} = [bestPhasing{3} allPhasings{ind(3)}];
            %iterNum = iterNum + 1;
        end

        allBest{muIndex} = positionBest; 
        if isempty(positionBest{1})
            allBestScores(muIndex) = -Inf;
            allBestPhasings{muIndex} = cell.empty;
        else

            if length(find(positionBest{1} == 0)) > 1
                allBestScores(muIndex) = -Inf;
            else
                allBestScores(muIndex)= prod(positionBestScores(:,1));
            end

            allBestPhasings{muIndex} = bestPhasing;
        end
        
        %separately save the total score
        all3BestScores{muIndex} = positionBestScores;
    end
    
    [~, ind] = max(allBestScores);
    allSolutions{sample} = allBest{ind};
    allSolutionsMu{sample} = allMu(validMuInd(ind),:);
    allSolutionsPhasing{sample} = allBestPhasings{ind};
    bestCSamples{sample} = allBest{ind};
    bestCSamplesScores{sample} = all3BestScores{ind};
end
cell2mat(allSolutionsMu)
%obtain the best 2 per sample
best2Ploidy = cell2mat(cellfun(@(x) [x{1} x{2}], bestCSamples, 'UniformOutput', 0));
best2Scores = cell2mat(cellfun(@(x) [x(:,1) x(:,2)], bestCSamplesScores, 'UniformOutput', 0));
%find the corresponding alleles
best2Alleles = cellfun(@(x) [x{1}' x{2}'], allSolutionsPhasing, 'UniformOutput', 0);

best2AllelesConverted = cell(length(uniqueSegments), (length(bestCSamples)*2));
for segmentInd = 1:length(uniqueSegments)
    multipleSampleInd = 1;
    for sample = 1:length(bestCSamples)

        sampleAlleles = best2Alleles{sample}(segmentInd,:);
        %alleleSolutionsConverted{sample} = sampleAlleles;
        best2AllelesConverted(segmentInd,multipleSampleInd) = sampleAlleles(1);
        best2AllelesConverted(segmentInd,multipleSampleInd+1) = sampleAlleles(2);
        multipleSampleInd = multipleSampleInd+2;
    end
end

profile on
%at every position, make the combinations
solutions = zeros(length(uniqueSegments), length(bestCSamples));
solutionsAlleles = cell(length(uniqueSegments), 1);
solutionScores = zeros(length(uniqueSegments), 1);
for segmentInd = 1:length(uniqueSegments)

    %combvecString = 'combvec(';
    allCombString = 'allcomb(';
    for sample = 1:length(bestCSamples)
        
        %obtain the 1,2,3 best for this segment
        positionBest = bestCSamples{sample};
        
        
        %combvecString = [combvecString, '[', int2str(positionBest{1}(segmentInd)), ', ', int2str(positionBest{2}(segmentInd)), ', ', int2str(positionBest{3}(segmentInd)), ']'];
        %combvecString = [combvecString, '[', int2str(positionBest{1}(segmentInd)), ', ', int2str(positionBest{2}(segmentInd)), ']'];
        allCombString = [allCombString, '[', int2str(positionBest{1}(segmentInd)), ', ', int2str(positionBest{2}(segmentInd)), ']'];
        if(sample < length(bestCSamples))
            
            %combvecString = [combvecString, ', '];
            allCombString = [allCombString, ', '];
        end
        
    end
   % combvecString = [combvecString, ')'];
    allCombString = [allCombString, ')'];
    %combinations = eval(combvecString);
    combinations = eval(allCombString);
    
    %make combinations
    %uniqueCombinations = unique(combinations', 'rows');
     normalClone = 2;
%     normalCloneRepeated = normalClone(ones(size(uniqueCombinations,1),1),:)';
%     fullCombinations = [normalCloneRepeated' uniqueCombinations];
    
    %score for this segment: associate the segment with the correct
    %combination in C
    
    currentBestScore = 0;
    currentBestSolution = zeros(1,length(bestCSamples));
    currentBestAlleles = [];
    %combinationAlleles = cell(length(combinations), 1);
    blockSize = 500000;
    totalIterations = length(combinations) / blockSize;
    iteration = 1;
    blockIndex = 1;
    for iterationIndex = 1:totalIterations+1

        if ((blockIndex+blockSize)) > length(combinations)
            blockCombinations = combinations(blockIndex:length(combinations),:);
        else
            %take these combinations from the array
            blockCombinations = combinations(blockIndex:(blockIndex+blockSize),:);
        end
        combinationScores = zeros(length(blockCombinations),1 );
        %compute the scores per block
         parfor combination = blockIndex:length(blockCombinations)
             LAFPs = zeros(1, length(bestCSamples));
%             %obtain the LAFP for this combination and sample
%             %sampleAlleles = cell.empty;

            %quickly obtain the correct score
            
            for sample = 1:length(bestCSamples) %ignore the first
                
                %what is the score for this sample?
                cn = blockCombinations(combination, sample);
                 sol1 = bestCSamples{sample}{1}(segmentInd);
                 sol2 = bestCSamples{sample}{2}(segmentInd);
                ind = 1;
                if(cn == sol2)
                   ind = 2; 
                end
                LAFPs(sample) = bestCSamplesScores{sample}(segmentInd, ind);
                %segmentC = cellfun(@(x) x(segmentInd), bestCSamples{sample}, 'UniformOutput', 0);
               % correctInd = cell2mat(segmentC) == cn;
                %segmentLAFP = bestCSamplesScores{sample}(segmentInd, correctInd);
                %LAFPs(sample) = segmentLAFP;
            end
             fullCombination = [normalClone combinations(combination,:)];
%             %compute the overall distance
             distance = 1 / (sum(pdist(fullCombination', 'cityblock')) +1);
%             %         %compute the score of the solution for a sample
             combinationScore = LAFPs * distance;
             combinationScores(combination) = prod(combinationScore);
            %dlmcell('solutionAlleles.txt', sampleAlleles, '-a'); %store as file, save memory
            
         end
        %find the best CN of this segment for each sample
        [currentBest, ind] = max(combinationScores);
        if currentBest > currentBestScore
           currentBestScore = currentBest;
           currentBestSolution = blockCombinations(ind,:);
           
           sampleAlleles = cell.empty;
           for sample = 1:length(bestCSamples) %ignore the first
               
               %what is the score for this sample?
               cn = blockCombinations(ind, sample);
               segmentC = cellfun(@(x) x(segmentInd), bestCSamples{sample}, 'UniformOutput', 0);
               correctInd = cell2mat(segmentC) == cn;
               if(sum(allSolutionsMu{sample} == [100 0]) > 1)
                   sampleAlleles{sample} = 'AB';
               else
                   sampleAlleles{sample} = allSolutionsPhasing{sample}{correctInd}{segmentInd};
               end
               
               
           end
           currentBestAlleles = sampleAlleles;

        end
        %bestCombination = blockCombinations(ind,:);
        %find the maximum inside this block
        %then use this to find the best. 
        blockIndex = blockIndex + (blockSize+1);
        iteration = iteration + 1;
    end
    clear combinations 
    
    %re-compute the alleles for this combination, is the fastest way
   
    %these are the solutions for all samples.
    solutions(segmentInd, :) = currentBestSolution;
    solutionsAlleles{segmentInd} = currentBestAlleles;
    solutionScores(segmentInd) = currentBestScore;
end

thresholdedSomVar = [0 1 0 1];
somVarDistanceMatrix = zeros(size(thresholdedSomVar,2), size(thresholdedSomVar,2));
% compute a distance matrix, 1 to 0 not possible. 
for sample1 = 1:size(thresholdedSomVar,2)
   sample1nan = isnan(thresholdedSomVar(:,sample1));
   for sample2 = 1:size(thresholdedSomVar,2)
        firstSample = thresholdedSomVar(:,sample1);
        secondSample = thresholdedSomVar(:,sample2);
        samplenan = isnan(thresholdedSomVar(:,sample2)) | isnan(thresholdedSomVar(:,sample1));
        firstSample(samplenan) = [];
        secondSample(samplenan) = [];
        totalDistance = 0;
        for varInd = 1:length(firstSample)
            if(firstSample(varInd) == 1 && secondSample(varInd) < 1)
                totalDistance = Inf;
                break;
            end
            distance = abs(firstSample(varInd) - secondSample(varInd));
            totalDistance = totalDistance + distance;
        end
        somVarDistanceMatrix(sample1, sample2) = (totalDistance / length(firstSample));
   end
end

%multiply matrices
updatedSomVarDistance = somVarDistanceMatrix;
updatedSomVarDistance(updatedSomVarDistance == 0) = 1;
updatedSomVarDistance(updatedSomVarDistance ~= 0 & ~isinf(updatedSomVarDistance)) = 1;
updatedSomVarDistance = updatedSomVarDistance - diag(diag(updatedSomVarDistance));
%convert alleles to columns as samples
alleleSolutionsConverted = cell(length(uniqueSegments), length(bestCSamples));
for segmentInd = 1:length(uniqueSegments)
    for sample = 1:length(bestCSamples)

        sampleAlleles = cellfun(@(x) x{sample}, solutionsAlleles, 'UniformOutput', 0);
        %alleleSolutionsConverted{sample} = sampleAlleles;
        alleleSolutionsConverted(:,sample) = sampleAlleles;
    end
end

T6107Alleles = alleleSolutionsConverted;
alleleDistanceMatrix = zeros(size(T6107Alleles,2), size(T6107Alleles,2));
for sampleInd = 1:size(T6107Alleles,2)
    for sampleInd2 = 1:size(T6107Alleles,2)
        sample1Phasing = T6107Alleles(:,sampleInd);
        %sample1StringPhasings = cellfun(@(x) char(x), sample1Phasing{1}, 'UniformOutput', 0);
        
        sample2Phasing = T6107Alleles(:,sampleInd2);
        %sample2StringPhasings = cellfun(@(x) char(x), sample2Phasing{1}, 'UniformOutput', 0);
        
        %compute the pairwise distances
        
        AMatches1 = strfind(sample1Phasing, 'A');
        BMatches1 = strfind(sample1Phasing, 'B');
        
        ALengths1 = cellfun(@length, AMatches1);
        BLengths1 = cellfun(@length, BMatches1);
        
        AMatches2 = strfind(sample2Phasing, 'A');
        BMatches2 = strfind(sample2Phasing, 'B');
        
        ALengths2 = cellfun(@length, AMatches2);
        BLengths2 = cellfun(@length, BMatches2);
        
        %compute sum of pairwise distances
        totalDistance = 0;
        for position = 1:length(ALengths1)
            
            if ALengths1(position) == 0 && ALengths2(position) > 0
                distance = Inf;
            elseif BLengths1(position) == 0 && BLengths2(position) > 0
                distance = Inf;
            else
                distance = abs(ALengths1(position) - ALengths2(position)) + ...
                    abs(BLengths1(position) - BLengths2(position)) + 1;
            end
            totalDistance = totalDistance + distance;
        end
        alleleDistanceMatrix(sampleInd, sampleInd2) = totalDistance;
    end
end
alleleDistanceMatrix(isinf(alleleDistanceMatrix)) = Inf;
alleleDistanceMatrix = alleleDistanceMatrix - diag(diag(alleleDistanceMatrix));

totalDistanceMatrix = alleleDistanceMatrix + updatedSomVarDistance;
totalDistanceMatrix(isnan(totalDistanceMatrix)) = Inf;
bestFit = totalDistanceMatrix;
bestFit(isinf(bestFit)) = 0;
bestFit(find(bestFit > 0)) = 1;
bestFitCell = {bestFit};
distanceCell = {totalDistanceMatrix};
bestTree = edmondsPhylogeny(bestFitCell, distanceCell);

% 
% 
%     fid = fopen('solutionAlleles.txt');
%     textFormat = repmat('%s', 1, length(bestCSamples));
%     line = textscan(fid, textFormat, 1, 'delimiter', '\t', 'HeaderLines', 5);
%     alleles = cellfun(@(x) x, line);



%store in file:
% all solutions, alleles, best 2 ploidy, best 2 alleles, best 2 scores. 
dlmcell('T6107AllSolutions.txt', alleleSolutionsConverted);

fnam='6107SomVar.txt';
SomVar=read_mixed_csv(fnam,'\t');
SomVar = cellfun(@(s) {str2double(s)},SomVar);
% make this a black and white model
thresholdedSomVar = cell2mat(SomVar);
% if > 0.01, accept
% if NA, make no comparison
thresholdedSomVar(thresholdedSomVar > 0.01) = 1;

somVarDistanceMatrix = zeros(length(thresholdedSomVar), length(thresholdedSomVar));
% compute a distance matrix, 1 to 0 not possible. 
for sample1 = 1:length(thresholdedSomVar)
   sample1nan = isnan(thresholdedSomVar(:,sample1));
   for sample2 = 1:length(thresholdedSomVar)
        firstSample = thresholdedSomVar(:,sample1);
        secondSample = thresholdedSomVar(:,sample2);
        samplenan = isnan(thresholdedSomVar(:,sample2)) | isnan(thresholdedSomVar(:,sample1));
        firstSample(samplenan) = [];
        secondSample(samplenan) = [];
        totalDistance = 0;
        for varInd = 1:length(firstSample)
            if(firstSample(varInd) == 1 && secondSample(varInd) < 1)
                totalDistance = Inf;
                break;
            end
            distance = abs(firstSample(varInd) - secondSample(varInd));
            totalDistance = totalDistance + distance;
        end
        somVarDistanceMatrix(sample1, sample2) = (totalDistance / length(firstSample));
   end
end

treeWeights = {somVarDistanceMatrix};
treeEdges = {somVarDistanceMatrix ~= 0};
[phylogenies, scores] = edmondsPhylogeny(treeEdges, treeWeights);

%make a simple tree for the distances
phylotree = seqlinkage(phylogenies{1})
view(phylotree)

newSolutions = allSolutions;
newSolutionsMu = allSolutionsMu;
newSolutionsPhasing = allSolutionsPhasing;
%write to a text file

allData = zeros(size(newSolutions{1},1), length(newSolutions));
allPhasingData = cell(1, length(newSolutions));
%collect all data in one big matrix
for i = 1:length(newSolutions)
    if ~isempty(newSolutions{i})
        allData(:,i) = newSolutions{i};
    end
    phasings = allSolutionsPhasing(i);
    stringPhasings = cellfun(@(x) char(x), phasings{1}, 'UniformOutput', 0);
    
    allPhasingData{i} = stringPhasings;
    
end
dlmwrite('2subclones_T3209_AllPositions.txt',allData,'delimiter','\t');

a = [];
for i = 1:length(allPhasingData)
    a = [a allPhasingData{i}']
end

dlmcell('2Subclones_T3209_AllPositions_alleles.txt',a); 

%% use the same model now, but instead of keeping 2 always keep 3 and compute the best based on that

