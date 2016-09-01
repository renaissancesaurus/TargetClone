function [  ] = startTargetClone( lafFile, vafFile, segmentationFile, sampleNamesFile, cOutputFile, aOutputFile, muOutputFile, treeOutputSuffix, kmin, kmax, metric, blockSize )
%STARTTARGETCLONE TargetClone method
%   Infers the most likely copy numbers, alleles and frequency per sample
%   from lesser allele frequency measurements. Also reconstructs a
%   subclonal evolution tree for the samples based on the copy numbers,
%   alleles, somatic variants, or a combination of the copy numbers and
%   somatic variants or the alleles and somatic variants. 

%% Step 1. Parse data from input files

%load the LAF
lafData = read_mixed_csv(lafFile,'\t');
lafData = cellfun(@(s) {str2double(s)},lafData);

%load the segmentation
fileID = fopen(segmentationFile);
C = textscan(fileID,'%s');
fclose(fileID);

segmentation = C{1};
uniqueSegments = unique(segmentation);

sampleNames = read_mixed_csv(sampleNamesFile, '\t');

%% Step 2. Initialize variables required by the method

cloneNum = 2; %default for now

%get all mu that sum to 100
allMu = permn(0:1:100, cloneNum); 
validMu = arrayfun(@(x)isequal (x, 100), sum(allMu,2)); 
validMuInd = find(validMu == 1);

%predefine all the possible alleles and distances, saves computational time
%per iteration
[alleleMap, allAlleleDistanceMatrix] = computeAllAlleleDistances(kmin, kmax);
%all possible copy numbers that need to be tested per segment
allPossibleCn = permn(kmin:kmax, cloneNum-1);

%% Step 3. Exhaustive search, first within each sample

%variables to store the solutions
allCSolutions = cell(1, size(lafData,2));
allMuSolutions = cell(1,size(lafData,2));
allAlleleSolutions = cell(1,size(lafData,2));

best2CPerSample = cell.empty;
best2CPerSamplesScores = cell.empty;

%find the best C, A and mu per sample
for sample = 1:size(lafData,2)
    %get the lesser allele frequencies per sample
    sampleLafData = cell2mat(lafData(:, sample));
    %if a measurement is missing, remove it and the segmentation entry too
    nanInd = isnan(sampleLafData);
    sampleLafData(nanInd) = [];
    newSegmentation = segmentation;
    newSegmentation(nanInd) = [];
    
    %group LAF measurements by segment
    segmentedLAF = cell(1, length(uniqueSegments));
    for i = 1:length(uniqueSegments)
        segmentLAFInd = strcmp(newSegmentation, uniqueSegments(i));
        segmentedLAF{i} = sampleLafData(segmentLAFInd);

    end
    
    %use the mean/median per segment or use all measurements
    %division is necessary for 3 decimal rounding
    LAF = cell(1, size(segmentedLAF, 2));
    for segmentInd = 1:size(segmentedLAF, 2)
        if strcmp(metric, 'individual') == 1
            LAF{segmentInd} = round((segmentedLAF{segmentInd}*1000)) / 1000.0;
        elseif strcmp(metric, 'mean') == 1
            LAF{segmentInd} = round((mean(segmentedLAF{segmentInd})*1000)) / 1000.0;
        else %median by default 
            LAF{segmentInd} = round((median(segmentedLAF{segmentInd})*1000)) / 1000.0;
        end
    end

    %Start exhaustive search
    
    %initialize cells to store solutions 
    best2CnPerMu = cell(1, length(validMuInd)); %allBest
    bestCnPerMuScores = zeros(1, length(validMuInd)); %allBestScores
    best2AllelesPerMu = cell(1, length(validMuInd)); %allBestPhasings
    best2CnPerMuScores = cell(1, length(validMuInd)); %all3BestScores
    parfor muIndex = 1:length(validMuInd)
        %store the two best copy numbers and alleles per segment 
        best2CnOfSegment = cell(1,2);
        best2AllelesOfSegment = cell(1,2);
        best2CnOfSegmentScores = zeros(length(uniqueSegments), 2);
        for segmentInd = 1:length(uniqueSegments)
            %if the entire segment has no measurements (missing values), 
            %use a default of 2 or 3 copies. 
            if isnan(LAF{segmentInd}) == size(LAF{1}, 2)
                best2CnOfSegmentScores(segmentInd,:) = 1; %do not include in weighing
                best2CnOfSegment{1} = [best2CnOfSegment{1}; 3];
                best2CnOfSegment{2} = [best2CnOfSegment{2}; 2];
                best2AllelesOfSegment{1} = [best2AllelesOfSegment{1} 'ABB'];
                best2AllelesOfSegment{2} = [best2AllelesOfSegment{2} 'AB'];
                continue;
            end
            %compute the probability of each possible copy number for a
            %segment in combination with our current mu given the LAF
            %measurement of that segment. 
            allPossibleCnScores = [];
            cnAssociatedAlleles = {};
            for possibleCn = 1:length(allPossibleCn)
                cn = allPossibleCn(possibleCn,:);
                [cnScore, bestAlleles] = computePCMuGivenLAF(cn, allMu(validMuInd(muIndex),:)/100, LAF{segmentInd}, alleleMap, allAlleleDistanceMatrix);
                allPossibleCnScores = [allPossibleCnScores; cnScore];
                cnAssociatedAlleles{possibleCn} = bestAlleles;
            end
            %sort the scores for each copy number to get the 2 best
            [sortedScores, ind] = sort(allPossibleCnScores(:,1),1, 'descend');
            best2CnOfSegmentScores(segmentInd,:) = sortedScores(1:2,1);
            %save the 2 best copy numbers and corresponding alleles
            best2CnOfSegment{1} = [best2CnOfSegment{1}; allPossibleCn(ind(1),:)];
            best2CnOfSegment{2} = [best2CnOfSegment{2}; allPossibleCn(ind(2),:)];
            best2AllelesOfSegment{1} = [best2AllelesOfSegment{1} cnAssociatedAlleles{ind(1)}];
            best2AllelesOfSegment{2} = [best2AllelesOfSegment{2} cnAssociatedAlleles{ind(2)}];
        end
        %save the best 2 copy numbers per mu
        best2CnPerMu{muIndex} = best2CnOfSegment; 
        %handle issues (should formally throw an error here)
        if isempty(best2CnOfSegment{1})
            bestCnPerMuScores(muIndex) = -Inf;
            best2AllelesPerMu{muIndex} = cell.empty;
        else
            %compute the score of the best copy number and this mu, this
            %score is used to choose a mu per sample
            if length(find(best2CnOfSegment{1} == 0)) > 1
                bestCnPerMuScores(muIndex) = -Inf;
            else
                bestCnPerMuScores(muIndex)= prod(best2CnOfSegmentScores(:,1));
            end
            %save the best 2 alleles at each segment for this mu
            best2AllelesPerMu{muIndex} = best2AllelesOfSegment;
        end
        
        %save the likelihood scores for the two best C and mu
        best2CnPerMuScores{muIndex} = best2CnOfSegmentScores;
    end
    %determine mu based on the best copy number for each mu, the best C and
    %mu that maximize the likelihood of observing the LAF measurements is
    %used to infer the best mu for this sample
    [~, ind] = max(bestCnPerMuScores);
    allCSolutions{sample} = best2CnPerMu{ind};
    allMuSolutions{sample} = allMu(validMuInd(ind),:);
    %save the additional 2 best solutions per sample
    allAlleleSolutions{sample} = best2AllelesPerMu{ind};
    best2CPerSample{sample} = best2CnPerMu{ind};
    best2CPerSamplesScores{sample} = best2CnPerMuScores{ind};
end

%% Step 4. Update results by including information across samples
% Makes all possible combinations with the 2 copy numbers per segment
% across all samples and recomputes P(C) based on these combinations.
% The final solution maximizes P(C,mu | LAF) per sample based on the
% recomputed P(C). 

C = zeros(length(uniqueSegments), length(best2CPerSample));
A = cell(length(uniqueSegments), 1);
solutionScores = zeros(length(uniqueSegments), 1);
for segmentInd = 1:length(uniqueSegments)

    %make a string for the allcomb function such that we generate all
    %possible combinations of copy numbers across samples based on the 2
    %best copy numbers per segment. Each combination is [Ci]. Allcomb is 
    %faster than the MATLAB builtin combvec function, but does the same thing. 
    allCombString = 'allcomb(';
    for sample = 1:length(best2CPerSample)
        
        best2C = best2CPerSample{sample};
        allCombString = [allCombString, '[', int2str(best2C{1}(segmentInd)), ', ', int2str(best2C{2}(segmentInd)), ']'];
        if(sample < length(best2CPerSample))
            allCombString = [allCombString, ', '];
        end
        
    end
   
    allCombString = [allCombString, ')'];
    %!!! Can use a lot of RAM for many samples !!! Not tested with more
    %than 25 samples. 
    combinations = eval(allCombString);
    
    %Compute the likelihood for each combination
    currentBestScore = 0;
    currentBestC = zeros(1,length(best2CPerSample));
    currentBestAlleles = [];
    
    % we use the blockSize parameter to ensure that this loop does not
    % use up all memory by computing the likelihood scores for smaller sets
    % of combinations and saving only the best solution per smaller set. 
    totalIterations = length(combinations) / blockSize;
    iteration = 1;
    blockIndex = 1;
    for iterationIndex = 1:totalIterations+1
        %we do not have blockSize more combinations, use everything until
        %the last combination
        if ((blockIndex+blockSize)) > length(combinations)
            blockCombinations = combinations(blockIndex:length(combinations),:);
        else
            %select a new set of combinations
            blockCombinations = combinations(blockIndex:(blockIndex+blockSize),:);
        end
        combinationScores = zeros(length(blockCombinations),1 );
        %compute the likelihood scores per block
        parfor combination = blockIndex:length(blockCombinations)
            LAFPs = zeros(1, length(best2CPerSample)); %P(LAF | C,mu)
            %compute P(LAF | C, mu) * P(C) per sample but use the whole
            %combination
            for sample = 1:length(best2CPerSample)
                
                %there are two possible copy numbers for each sample
                %computed in the previous step, find which one it is in
                %this current combination
                cn = blockCombinations(combination, sample);
                sol1 = best2CPerSample{sample}{1}(segmentInd);
                sol2 = best2CPerSample{sample}{2}(segmentInd);
                ind = 1;
                if(cn == sol2)
                    ind = 2;
                end
                %obtain the P(LAF|C,mu) for this copy number 
                LAFPs(sample) = best2CPerSamplesScores{sample}(segmentInd, ind);
                
            end
            %With all P(LAF|C,mu) known, we can compute P(C) for the
            %combination
            
            %Add a copy number of 2, corresponding to the normal component.
            %This ensures that we optimize the distance to the normal
            %sample too as much as possible. 
            fullCombination = [2 combinations(combination,:)];
            %compute P(C)
            distance = 1 / (sum(pdist(fullCombination', 'cityblock')) +1);
            %compute P(C,mu | LAF)
            combinationScore = LAFPs * distance;
            %combination score is for all samples, so take the product
            %to compute a P(C,mu | LAF) across all samples
            combinationScores(combination) = prod(combinationScore);
        end
        %find the best CN of this segment for each sample
        [currentBest, ind] = max(combinationScores);
        %if it is better than for the previous subset of combinations, keep
        %this C and A as the best. 
        if currentBest > currentBestScore
           currentBestScore = currentBest;
           currentBestC = blockCombinations(ind,:);
           
           
           %Associate the C with the correct alleles (depending on if 1st or 2nd is chosen)
           sampleAlleles = cell.empty;
           for sample = 1:length(best2CPerSample)
               
               cn = blockCombinations(ind, sample);
               segmentC = cellfun(@(x) x(segmentInd), best2CPerSample{sample}, 'UniformOutput', 0);
               correctInd = cell2mat(segmentC) == cn;
               %if only normal cells are inferred to be present in the
               %sample, set the 'tumor subclone' alleles to 'AB' as well
               if(sum(allMuSolutions{sample} == [100 0]) > 1)
                   sampleAlleles{sample} = 'AB';
               else
                   sampleAlleles{sample} = allAlleleSolutions{sample}{correctInd}{segmentInd};
               end
               
               
           end
           currentBestAlleles = sampleAlleles;

        end
        blockIndex = blockIndex + (blockSize+1);
        iteration = iteration + 1;
    end
    %remove the combinations to clear up memory
    clear combinations 
    
    %these are the solutions for all samples for the current segment.
    C(segmentInd, :) = currentBestC;
    A{segmentInd} = currentBestAlleles;
    solutionScores(segmentInd) = currentBestScore;
end

%% Step 5. Generate subclonal evolution trees
generateSubclonalTrees(C, A, vafFile, sampleNames, treeOutputSuffix); %the function itself generates the trees

%% Step 6. Generate output
%C
%add sample names to header
fid = fopen(cOutputFile, 'w');
if fid == -1; error('Cannot open file: %s', outfile); end
fprintf(fid, '%s\t', sampleNames{:});
fprintf(fid, '\n');
fclose(fid);
dlmwrite(cOutputFile,C,'delimiter','\t','-append');

%Convert A to better format for output
allelesPerSample = cell(length(A), length(A{1}));
for i = 1:length(A)
    segmentA = A{i};
    for j = 1:length(segmentA)
        allelesPerSample{i,j} =  A{i}{j};
    end
end

fid = fopen(aOutputFile, 'w');
if fid == -1; error('Cannot open file: %s', outfile); end
fprintf(fid, '%s\t', sampleNames{:});
fprintf(fid, '\n');
fclose(fid);
dlmcell(aOutputFile,allelesPerSample,'delimiter','\t','-a');
%convert mu as well
muPerSample = zeros(2, size(allMuSolutions,2));
for i = 1:length(allMuSolutions)
    for j = 1:length(allMuSolutions{i})
        muPerSample(j,i) = allMuSolutions{i}(j);
    end
end

fid = fopen(muOutputFile, 'w');
if fid == -1; error('Cannot open file: %s', outfile); end
fprintf(fid, '%s\t', sampleNames{:});
fprintf(fid, '\n');
fclose(fid);
dlmwrite(muOutputFile,muPerSample,'delimiter','\t', '-append');
end

