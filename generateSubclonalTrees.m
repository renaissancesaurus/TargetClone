function [ ] = generateSubclonalTrees( C, A, somVarFile, sampleNames, treeOutputSuffix )
%GENERATESUBCLONALTREES Generates subclonal evolution trees
%   The trees are automatically generated for all combinations of C, A and
%   somatic variants

%% Requirements
SomVar = read_mixed_csv(somVarFile, '\t');
SomVar = cell2mat(cellfun(@(s) {str2double(s)}, SomVar));
% make sure that the somatic variants are always in binary format
thresholdedSomVar = SomVar;
thresholdedSomVar(thresholdedSomVar > 0.01) = 1;


%% 1. Distance matrix for copy numbers
cDistanceMatrix = zeros(size(C,2), size(C,2));
for sampleInd = 1:size(C,2)
    for sampleInd2 = 1:size(C,2)
        sample1 = C(:, sampleInd);
        
        sample2 = C(:,sampleInd2);
        
        %compute sum of pairwise distances
        totalDistance = 0;
        for position = 1:length(sample1)
            
            if sample1(position) == 0 && sample2(position) > 0
                distance = Inf;
            else
                distance = abs(sample1(position) - sample2(position));
            end
            totalDistance = totalDistance + distance;
        end
        cDistanceMatrix(sampleInd, sampleInd2) = totalDistance;
    end
end
cDistanceMatrix(isinf(cDistanceMatrix)) = Inf;

%% 2. Distance matrix for alleles

%convert the alleles to the correct format so that we can obtain it per
%sample for all segments
allelesPerSample = cell(length(A), length(A{1}));
for i = 1:length(A)
    segmentA = A{i};
    for j = 1:length(segmentA)
        allelesPerSample{i,j} =  A{i}{j};
    end
end

alleleDistanceMatrix = zeros(size(A,1), size(A,1));
for sampleInd = 1:size(A,1)
    for sampleInd2 = 1:size(A,1)
        sample1Phasing = {allelesPerSample{:,sampleInd}};
        sample2Phasing = {allelesPerSample{:,sampleInd2}};
        
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
                    abs(BLengths1(position) - BLengths2(position));
            end
            totalDistance = totalDistance + distance;
        end
        alleleDistanceMatrix(sampleInd, sampleInd2) = totalDistance;
    end
end
alleleDistanceMatrix(isinf(alleleDistanceMatrix)) = Inf;

%% 3. Distance matrix for somatic variants

somVarDistanceMatrix = zeros(size(thresholdedSomVar,2), size(thresholdedSomVar,2));
% compute a distance matrix, 1 to 0 (loss of somatic variant) not allowed. 
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

updatedSomVarDistance = somVarDistanceMatrix;
updatedSomVarDistance(updatedSomVarDistance == 0) = 1;
updatedSomVarDistance(updatedSomVarDistance ~= 0 & ~isinf(updatedSomVarDistance)) = 1;

%% Reconstructing the actual trees

%% C only
bestFit = cDistanceMatrix;
bestFit(isinf(bestFit)) = 0;
bestFit(find(bestFit > 0)) = 1;
bestFitCell = {bestFit};
distanceCell = {cDistanceMatrix};
bestTree = edmondsPhylogeny(bestFitCell, distanceCell, sampleNames);

%print to file
child_handles = allchild(0);
names = get(child_handles,'Name');
k = find(strncmp('Biograph Viewer', names, 15));
print(child_handles(k),'-dpng', '-r0',strcat(treeOutputSuffix{1}, 'cOnly', treeOutputSuffix{2}, '.png'));
close(child_handles(k))

%% A only 
bestFit = alleleDistanceMatrix;
bestFit(isinf(bestFit)) = 0;
bestFit(find(bestFit > 0)) = 1;
bestFitCell = {bestFit};
distanceCell = {alleleDistanceMatrix};
bestTree = edmondsPhylogeny(bestFitCell, distanceCell, sampleNames);

%print to file
child_handles = allchild(0);
names = get(child_handles,'Name');
k = find(strncmp('Biograph Viewer', names, 15));
print(child_handles(k),'-dpng', '-r0',strcat(treeOutputSuffix{1}, 'aOnly', treeOutputSuffix{2}, '.png'));
close(child_handles(k))

%% Somatic variants only
bestFit = somVarDistanceMatrix;
bestFit(isinf(bestFit)) = 0;
bestFit(find(bestFit > 0)) = 1;
bestFitCell = {bestFit};
distanceCell = {somVarDistanceMatrix};
bestTree = edmondsPhylogeny(bestFitCell, distanceCell, sampleNames);

%print to file
child_handles = allchild(0);
names = get(child_handles,'Name');
k = find(strncmp('Biograph Viewer', names, 15));
print(child_handles(k),'-dpng', '-r0',strcat(treeOutputSuffix{1}, 'somVarOnly', treeOutputSuffix{2}, '.png'));
close(child_handles(k))

%% C + somatic variants
totalCSDistanceMatrix = cDistanceMatrix .* updatedSomVarDistance;

bestFit = totalCSDistanceMatrix;
bestFit(isinf(bestFit)) = 0;
bestFit(find(bestFit > 0)) = 1;
bestFitCell = {bestFit};
distanceCell = {totalCSDistanceMatrix};
bestTree = edmondsPhylogeny(bestFitCell, distanceCell, sampleNames);

%print to file
child_handles = allchild(0);
names = get(child_handles,'Name');
k = find(strncmp('Biograph Viewer', names, 15));
print(child_handles(k),'-dpng', '-r0',strcat(treeOutputSuffix{1}, 'cAndSomVar', treeOutputSuffix{2}, '.png'));
close(child_handles(k))

%% A + somatic variants

totalASDistanceMatrix = alleleDistanceMatrix .* updatedSomVarDistance;
bestFit = totalASDistanceMatrix;
bestFit(isinf(bestFit)) = 0;
bestFit(find(bestFit > 0)) = 1;
bestFitCell = {bestFit};
distanceCell = {totalASDistanceMatrix};
bestTree = edmondsPhylogeny(bestFitCell, distanceCell, sampleNames);

%print to file
child_handles = allchild(0);
names = get(child_handles,'Name');
k = find(strncmp('Biograph Viewer', names, 15));
print(child_handles(k),'-dpng', '-r0',strcat(treeOutputSuffix{1}, 'aAndSomVar', treeOutputSuffix{2}, '.png'));
close(child_handles(k))

end

