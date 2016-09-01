function [] = validateTrees(finalC, finalA, realTrees, realSomVar, finalSSC, finalSSA)


noiseLevels = [0, 0.02, 0.04, 0.06, 0.08, 0.1];
%% Copy numbers only

noiseScore = cell.empty;
for i = 1:length(finalC) %for every noise level
    
    noiseC = finalC{i}; %all iterations
    
    iterationScores = zeros(length(noiseC), 1);
    for iteration = 1:length(noiseC) %per iteration
        realPhylogeny = realTrees{i}{iteration};

        %compute the paiwise distances between all samples
        
        %compute the distance per position
        %then sum the distances
        iterationC = noiseC{iteration};
        iterationScore = 0;
        
        distanceMatrix = zeros(size(iterationC,2), size(iterationC,2));
        for sampleInd = 1:size(iterationC,2)
            for sampleInd2 = 1:size(iterationC,2)
                sample1C = iterationC(:, sampleInd);
                sample2C= iterationC(:,sampleInd2);
              
                %compute sum of pairwise distances
                totalDistance = 0;
                for position = 1:length(sample1C)
                    
                    if sample1C(position) == 0 && sample2C(position) > 0
                        distance = Inf;
                    else
                        distance = abs(sample1C(position) - sample2C(position));
                    end
                    totalDistance = totalDistance + distance;
                end
                distanceMatrix(sampleInd2, sampleInd) = totalDistance;
            end
        end
        distanceMatrix(isinf(distanceMatrix)) = Inf;
        
        %If no relations are possible between samples
        testTree = distanceMatrix;
        testTree(isinf(testTree)) = 0;
        if (max(testTree(:)) < 1)
            iterationScore = iterationScore + 0;
            iterationScores(iteration) = score;
            continue; 
        end
        
        bestFit = distanceMatrix;
        bestFit(isinf(bestFit)) = 0;
        bestFit(find(bestFit > 0)) = 1;
        bestFitCell = {bestFit};
        distanceCell = {distanceMatrix};
        bestTree = edmondsPhylogeny(bestFitCell, distanceCell, {});

        child_handles = allchild(0);
        names = get(child_handles,'Name');
        k = find(strncmp('Biograph Viewer', names, 15));
        close(child_handles(k))
        bestTree = bestTree{1};
        % check if the tree is of the same size
        bestTree(bestTree > 0) = 1;
        if size(bestTree, 1) < size(realPhylogeny)
            iterationScore = iterationScore + 0;
            iterationScores(iteration) = score;
           continue; 
        end
        %Compute how often parents are switched
        difference = bestTree - realPhylogeny;
        distance = difference(difference < 0);
        totalDistance = sum(abs(distance));
        score = 1 - (totalDistance / size(realPhylogeny,2));
        iterationScore = iterationScore + score;
        iterationScores(iteration) = score;
    end
    noiseScore{i} = iterationScores;
end

% Plot 
boxplot(cell2mat(noiseScore), noiseLevels, 'colors', [0 0 0])

title('C');
xlabel('Noise (?)');
ylabel('Tree reconstruction accuracy');
figure;

%% A only

noiseScore = cell.empty;
for i = 1:length(finalA) %for every noise level
    
    noiseC = finalA{i}; %all iterations
    
    iterationScores = zeros(length(noiseC), 1);
    for iteration = 1:length(noiseC)
        realPhylogeny = realTrees{i}{iteration};

        iterationPhasing = noiseC{iteration};
        allSamplePhasing = cell(length(iterationPhasing{1}), 1);
        iterationScore = 0;
        for sample = 1:length(iterationPhasing{1})
            
            samplePhasings = cellfun(@(x) cell2mat(x(sample)), iterationPhasing, 'UniformOutput', 0);
            allSamplePhasing{sample} = samplePhasings;
        end
        distanceMatrix = zeros(size(iterationPhasing,2), size(iterationPhasing,2));
        for sampleInd = 1:length(allSamplePhasing)
            for sampleInd2 = 1:length(allSamplePhasing)
                sample1Phasing = allSamplePhasing{sampleInd};
                sample2Phasing = allSamplePhasing{sampleInd2};
                
                %compute the pairwise distances based on alleles

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
                distanceMatrix(sampleInd2, sampleInd) = totalDistance;
            end
        end
        distanceMatrix(isinf(distanceMatrix)) = Inf;
        
        %check if any relations are possible
        testTree = distanceMatrix;
        testTree(isinf(testTree)) = 0;
        if (max(testTree(:)) < 1)
            iterationScore = iterationScore + 0;
            iterationScores(iteration) = score;
            continue; 
        end
        
        bestFit = distanceMatrix;
        bestFit(isinf(bestFit)) = 0;
        bestFit(find(bestFit > 0)) = 1;
        bestFitCell = {bestFit};
        distanceCell = {distanceMatrix};
        bestTree = edmondsPhylogeny(bestFitCell, distanceCell, {});
        
        child_handles = allchild(0);
        names = get(child_handles,'Name');
        k = find(strncmp('Biograph Viewer', names, 15));
        close(child_handles(k))
        bestTree = bestTree{1};
        % check if the trees are of same length
        bestTree(bestTree > 0) = 1;
        if size(bestTree, 1) < size(realPhylogeny)
            iterationScore = iterationScore + 0;
            iterationScores(iteration) = score;
           continue; 
        end
        %compute how many parents match
        difference = bestTree - realPhylogeny;
        distance = difference(difference < 0);
        totalDistance = sum(abs(distance));
        score = 1 - (totalDistance / size(realPhylogeny,2));
        iterationScore = iterationScore + score;
        iterationScores(iteration) = score;
    end
    noiseScore{i} = iterationScores;
end

boxplot(cell2mat(noiseScore), noiseLevels, 'colors', [0 0 0])

title('A');
xlabel('Noise (?)');
ylabel('Tree reconstruction accuracy');
figure;

%% somatic variants only

noiseScore = cell.empty;
for i = 1:length(realSomVar) %for every noise level
    
    noiseSomVar = realSomVar{i}; %all iterations
    
    iterationScores = zeros(length(noiseSomVar), 1);
    for iteration = 1:length(noiseSomVar)
        
        iterationSamples = noiseSomVar{iteration};
        
        %compute the all-to-all distance for these samples
        iterationScore = 0;
        
        distanceMatrix = zeros(size(iterationSamples,2), size(iterationSamples,2));
        for sampleInd = 1:size(iterationSamples,2)
            for sampleInd2 = 1:size(iterationSamples,2)
                firstSample = iterationSamples{sampleInd};
                secondSample = iterationSamples{sampleInd2};
                totalDistance = 0;
                for varInd = 1:length(firstSample)
                    if(firstSample(varInd) == 1 && secondSample(varInd) < 1)
                        totalDistance = Inf;
                        break;
                    end
                    distance = abs(firstSample(varInd) - secondSample(varInd));
                    totalDistance = totalDistance + distance;
                end
                distanceMatrix(sampleInd, sampleInd2) = totalDistance;
            end
        end

        treeWeights = {distanceMatrix};
        treeEdges = {distanceMatrix ~= 0};
        [phylogenies, scores] = edmondsPhylogeny(treeEdges, treeWeights, {});
        
        child_handles = allchild(0);
        names = get(child_handles,'Name');
        k = find(strncmp('Biograph Viewer', names, 15));
        close(child_handles(k))
        bestTree = phylogenies{1};
        % check if the trees are of the same size
        bestTree(bestTree > 0) = 1;
        if size(bestTree, 1) < size(realPhylogeny)
            iterationScore = iterationScore + 0;
            iterationScores(iteration) = score;
           continue; 
        end
        %Compute how many parents match
        difference = bestTree - realPhylogeny;
        distance = difference(difference < 0);
        totalDistance = sum(abs(distance));
        score = 1 - (totalDistance / size(realPhylogeny,2));
        iterationScore = iterationScore + score;
        iterationScores(iteration) = score;
    end
    noiseScore{i} = iterationScores;
end
boxplot(cell2mat(noiseScore), noiseLevels, 'colors', [0 0 0])

title('Somatic variants');
xlabel('Noise (?)');
ylabel('Tree reconstruction accuracy');
figure;
%% alleles + somatic variants

%multiply the matrices

noiseScore = cell.empty;
for i = 1:length(finalA) %for every noise level
    
    noiseC = finalA{i}; %all iterations
    
    noiseSomVar = realSomVar{i}; %all iterations
    iterationScores = zeros(length(noiseC), 1);
    for iteration = 1:length(noiseC)
        realPhylogeny = realTrees{i}{iteration};

        iterationSamples = noiseSomVar{iteration};
        %compute the paiwise distances between all samples
       
        %compute the distance per position
        %then sum the distances
        iterationPhasing = noiseC{iteration};
        %put alleles in different format for simplicity
        allSamplePhasing = cell(length(iterationPhasing{1}), 1);
        iterationScore = 0;
        for sample = 1:length(iterationPhasing{1})
            samplePhasings = cellfun(@(x) cell2mat(x(sample)), iterationPhasing, 'UniformOutput', 0);
            allSamplePhasing{sample} = samplePhasings;
        end
        distanceMatrix = zeros(size(iterationPhasing,2), size(iterationPhasing,2));
        for sampleInd = 1:length(allSamplePhasing)
            for sampleInd2 = 1:length(allSamplePhasing)
                sample1Phasing = allSamplePhasing{sampleInd};
                sample2Phasing = allSamplePhasing{sampleInd2};
                
                % Compute allele distances
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
                distanceMatrix(sampleInd2, sampleInd) = totalDistance;
            end
        end
        distanceMatrix(isinf(distanceMatrix)) = Inf;
        
        %compute distance matrix for somatic variants
        somvarDistanceMatrix = zeros(size(iterationSamples,2), size(iterationSamples,2));
        for sampleInd = 1:size(iterationSamples,2)
            for sampleInd2 = 1:size(iterationSamples,2)
                firstSample = iterationSamples{sampleInd};
                secondSample = iterationSamples{sampleInd2};
                totalDistance = 0;
                for varInd = 1:length(firstSample)
                    if(firstSample(varInd) == 1 && secondSample(varInd) < 1)
                        totalDistance = Inf;
                        break;
                    end
                    distance = abs(firstSample(varInd) - secondSample(varInd));
                    totalDistance = totalDistance + distance;
                end
                somvarDistanceMatrix(sampleInd, sampleInd2) = totalDistance;
            end
        end
        
        %multiply matrices
        totalDistanceMatrix = distanceMatrix .* somvarDistanceMatrix;
        
        %Check if any relations are possible
        testTree = totalDistanceMatrix;
        testTree(isinf(testTree)) = 0;
        if (max(testTree(:)) < 1)
            iterationScore = iterationScore + 0;
            iterationScores(iteration) = score;
            continue; 
        end
        
        bestFit = totalDistanceMatrix;
        bestFit(isinf(bestFit)) = 0;
        bestFit(find(bestFit > 0)) = 1;
        bestFitCell = {bestFit};
        distanceCell = {totalDistanceMatrix};
        bestTree = edmondsPhylogeny(bestFitCell, distanceCell, {});
        
        child_handles = allchild(0);
        names = get(child_handles,'Name');
        k = find(strncmp('Biograph Viewer', names, 15));
        close(child_handles(k))
        bestTree = bestTree{1};
        % check if trees are of the same length
        bestTree(bestTree > 0) = 1;
        if size(bestTree, 1) < size(realPhylogeny)
            iterationScore = iterationScore + 0;
            iterationScores(iteration) = score;
           continue; 
        end
        %Score how many parents are the same
        difference = bestTree - realPhylogeny;
        distance = difference(difference < 0);
        totalDistance = sum(abs(distance));
        score = 1 - (totalDistance / size(realPhylogeny,2));
        iterationScore = iterationScore + score;
        iterationScores(iteration) = score;
    end
    noiseScore{i} = iterationScores;
end

boxplot(cell2mat(noiseScore), noiseLevels, 'colors', [0 0 0])

title('A and somatic variants');
xlabel('Noise (?)');
ylabel('Tree reconstruction accuracy');
figure;

%% C + somatic variants
noiseScore = cell.empty;
for i = 1:length(finalC) %for every noise level
    
    noiseC = finalC{i}; %all iterations
    
    noiseSomVar = realSomVar{i}; %all iterations
    iterationScores = zeros(length(noiseC), 1);
    for iteration = 1:length(noiseC)
        realPhylogeny = realTrees{i}{iteration};

        iterationSamples = noiseSomVar{iteration};
        iterationC = noiseC{iteration};

        iterationScore = 0;

        distanceMatrix = zeros(size(iterationC,2), size(iterationC,2));
        for sampleInd = 1:size(iterationC,2)
            for sampleInd2 = 1:size(iterationC,2)
                sample1C = iterationC(:, sampleInd);
                sample2C = iterationC(:,sampleInd2);

                %compute sum of pairwise distances
                totalDistance = 0;
                for position = 1:length(sample1Phasing)
                    
                    if sample1C(position) == 0 && sample2C(position) > 0
                        distance = Inf;
                    else
                        distance = abs(sample1C(position) - sample2C(position));
                    end
                    totalDistance = totalDistance + distance;
                end
                distanceMatrix(sampleInd2, sampleInd) = totalDistance;
            end
        end
        distanceMatrix(isinf(distanceMatrix)) = Inf;
        
        %compute distance matrix for somatic variants
        somvarDistanceMatrix = zeros(size(iterationSamples,2), size(iterationSamples,2));
        for sampleInd = 1:size(iterationSamples,2)
            for sampleInd2 = 1:size(iterationSamples,2)
                firstSample = iterationSamples{sampleInd};
                secondSample = iterationSamples{sampleInd2};
                totalDistance = 0;
                for varInd = 1:length(firstSample)
                    if(firstSample(varInd) == 1 && secondSample(varInd) < 1)
                        totalDistance = Inf;
                        break;
                    end
                    distance = abs(firstSample(varInd) - secondSample(varInd));
                    totalDistance = totalDistance + distance;
                end
                somvarDistanceMatrix(sampleInd, sampleInd2) = totalDistance;
            end
        end
        
        %multiply matrices
        totalDistanceMatrix = distanceMatrix .* somvarDistanceMatrix;
        
        %Check if any relations are possible
        testTree = totalDistanceMatrix;
        testTree(isinf(testTree)) = 0;
        if (max(testTree(:)) < 1)
            iterationScore = iterationScore + 0;
            iterationScores(iteration) = score;
            continue; 
        end
        
        bestFit = totalDistanceMatrix;
        bestFit(isinf(bestFit)) = 0;
        bestFit(find(bestFit > 0)) = 1;
        bestFitCell = {bestFit};
        distanceCell = {totalDistanceMatrix};
        bestTree = edmondsPhylogeny(bestFitCell, distanceCell, {});
        
        child_handles = allchild(0);
        names = get(child_handles,'Name');
        k = find(strncmp('Biograph Viewer', names, 15));
        close(child_handles(k))
        bestTree = bestTree{1};
        % check if the trees are of the same length
        bestTree(bestTree > 0) = 1;
        if size(bestTree, 1) < size(realPhylogeny)
            iterationScore = iterationScore + 0;
            iterationScores(iteration) = score;
           continue; 
        end
        %Score how many parents are the same
        difference = bestTree - realPhylogeny;
        distance = difference(difference < 0);
        totalDistance = sum(abs(distance));
        score = 1 - (totalDistance / size(realPhylogeny,2));
        iterationScore = iterationScore + score;
        iterationScores(iteration) = score;
    end
    noiseScore{i} = iterationScores;
    %noiseScore(i) = iterationScore / length(noisePhasing);
end

boxplot(cell2mat(noiseScore), noiseLevels, 'colors', [0 0 0])

title('Somatic variants and ploidy');
xlabel('Noise (?)');
ylabel('Tree reconstruction accuracy');


%% SINGLE SAMPLE
%% reconstruct based on just copy numbers

noiseScore = cell.empty;
for i = 1:length(finalSSC) %for every noise level
    
    noiseC = finalSSC{i}; %all iterations
    
    iterationScores = zeros(length(noiseC), 1);
    for iteration = 1:length(noiseC)
        realPhylogeny = realTrees{i}{iteration};

        %compute the paiwise distances between all samples
        iterationC = noiseC{iteration};
        iterationScore = 0;
        distanceMatrix = zeros(size(iterationC,2), size(iterationC,2));
        for sampleInd = 1:size(iterationC,2)
            for sampleInd2 = 1:size(iterationC,2)
                sample1C = iterationC(:, sampleInd);
               
                sample2C = iterationC(:,sampleInd2);
                %compute sum of pairwise distances
                totalDistance = 0;
                for position = 1:length(sample1C)
                    
                    if sample1C(position) == 0 && sample2C(position) > 0
                        distance = Inf;
                    else
                        distance = abs(sample1C(position) - sample2C(position));
                    end
                    totalDistance = totalDistance + distance;
                end
                distanceMatrix(sampleInd2, sampleInd) = totalDistance;
            end
        end
        distanceMatrix(isinf(distanceMatrix)) = Inf;
        
        testTree = distanceMatrix;
        testTree(isinf(testTree)) = 0;
        if (max(testTree(:)) < 1)
            iterationScore = iterationScore + 0;
            iterationScores(iteration) = score;
            continue; 
        end
        
        bestFit = distanceMatrix;
        bestFit(isinf(bestFit)) = 0;
        bestFit(find(bestFit > 0)) = 1;
        bestFitCell = {bestFit};
        distanceCell = {distanceMatrix};
        bestTree = edmondsPhylogeny(bestFitCell, distanceCell, {});
        
        child_handles = allchild(0);
        names = get(child_handles,'Name');
        k = find(strncmp('Biograph Viewer', names, 15));
        close(child_handles(k))
        bestTree = bestTree{1};
        
        bestTree(bestTree > 0) = 1;
        if size(bestTree, 1) < size(realPhylogeny)
            iterationScore = iterationScore + 0;
            iterationScores(iteration) = score;
           continue; 
        end
        difference = bestTree - realPhylogeny;
        distance = difference(difference < 0);
        totalDistance = sum(abs(distance));
        score = 1 - (totalDistance / size(realPhylogeny,2));
        iterationScore = iterationScore + score;
        iterationScores(iteration) = score;
    end
    noiseScore{i} = iterationScores;
end

%Plot

boxplot(cell2mat(noiseScore), noiseLevels, 'colors', [0 0 0])

title('C');
xlabel('Noise (?)');
ylabel('Tree reconstruction accuracy');
figure;

%% A only

noiseScore = cell.empty;
for i = 1:length(finalSSA) %for every noise level
    
    noiseC = finalSSA{i}; %all iterations
    
    iterationScores = zeros(length(noiseC), 1);
    for iteration = 1:length(noiseC)
        realPhylogeny = realTrees{i}{iteration};

        iterationSamples = noiseSomVar{iteration};
        %compute the paiwise distances between all samples
        
        iterationPhasing = noiseC{iteration};
        allSamplePhasing = iterationPhasing;
        
        distanceMatrix = zeros(size(iterationPhasing,2), size(iterationPhasing,2));
        for sampleInd = 1:length(allSamplePhasing)
            for sampleInd2 = 1:length(allSamplePhasing)
                sample1Phasing = allSamplePhasing{sampleInd}{1}{1};
                sample2Phasing = allSamplePhasing{sampleInd2}{1}{1};
                
                %compute the pairwise distances based on alleles
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
                distanceMatrix(sampleInd2, sampleInd) = totalDistance;
            end
        end
        distanceMatrix(isinf(distanceMatrix)) = Inf;

       
        bestFit = distanceMatrix;
        bestFit(isinf(bestFit)) = 0;
        bestFit(find(bestFit > 0)) = 1;
        bestFitCell = {bestFit};
        distanceCell = {totalDistanceMatrix};
        bestTree = edmondsPhylogeny(bestFitCell, distanceCell, {});
        
        child_handles = allchild(0);
        names = get(child_handles,'Name');
        k = find(strncmp('Biograph Viewer', names, 15));
        close(child_handles(k))
        if isempty(bestTree)
            iterationScore = iterationScore + 0;
            iterationScores(iteration) = score;
            continue; 
        end
        bestTree = bestTree{1};
        % check how many positions overlap that are not inf
        bestTree(bestTree > 0) = 1;
        if size(bestTree, 1) < size(realPhylogeny)
            iterationScore = iterationScore + 0;
            iterationScores(iteration) = score;
           continue; 
        end
        difference = bestTree - realPhylogeny;
        distance = difference(difference < 0);
        totalDistance = sum(abs(distance));
        score = 1 - (totalDistance / size(realPhylogeny,2));
        iterationScore = iterationScore + score;
        iterationScores(iteration) = score;
    end
    noiseScore{i} = iterationScores;
end

boxplot(cell2mat(noiseScore), noiseLevels, 'colors', [0 0 0])

title('A');
xlabel('Noise (?)');
ylabel('Tree reconstruction accuracy');
figure;

%% somatic variants only

noiseScore = cell.empty;
for i = 1:length(realSomVar) %for every noise level
    
    noiseSomVar = realSomVar{i}; %all iterations
    
    iterationScores = zeros(length(noiseSomVar), 1);
    for iteration = 1:length(noiseSomVar)
        
        iterationSamples = noiseSomVar{iteration};
        
        %compute the all-to-all distance for these samples
        iterationScore = 0;
        
        distanceMatrix = zeros(size(iterationSamples,2), size(iterationSamples,2));
        for sampleInd = 1:size(iterationSamples,2)
            for sampleInd2 = 1:size(iterationSamples,2)
                firstSample = iterationSamples{sampleInd};
                secondSample = iterationSamples{sampleInd2};
                totalDistance = 0;
                for varInd = 1:length(firstSample)
                    if(firstSample(varInd) == 1 && secondSample(varInd) < 1)
                        totalDistance = Inf;
                        break;
                    end
                    distance = abs(firstSample(varInd) - secondSample(varInd));
                    totalDistance = totalDistance + distance;
                end
                distanceMatrix(sampleInd, sampleInd2) = totalDistance;
            end
        end

        treeWeights = {distanceMatrix};
        treeEdges = {distanceMatrix ~= 0};
        [phylogenies, scores] = edmondsPhylogeny(treeEdges, treeWeights, {});
        
        child_handles = allchild(0);
        names = get(child_handles,'Name');
        k = find(strncmp('Biograph Viewer', names, 15));
        close(child_handles(k))
        bestTree = phylogenies{1};
        % check how many positions overlap that are not inf
        bestTree(bestTree > 0) = 1;
        if size(bestTree, 1) < size(realPhylogeny)
            iterationScore = iterationScore + 0;
            iterationScores(iteration) = score;
           continue; 
        end
        difference = bestTree - realPhylogeny;
        distance = difference(difference < 0);
        totalDistance = sum(abs(distance));
        score = 1 - (totalDistance / size(realPhylogeny,2));
        iterationScore = iterationScore + score;
        iterationScores(iteration) = score;
    end
    noiseScore{i} = iterationScores;
end
boxplot(cell2mat(noiseScore), noiseLevels, 'colors', [0 0 0])

title('Somatic variants');
xlabel('Noise (?)');
ylabel('Tree reconstruction accuracy');
figure;
%% A + somatic variants

noiseScore = cell.empty;
for i = 1:length(finalSSA) %for every noise level
    
    noiseC = finalSSA{i}; %all iterations
    
    noiseSomVar = realSomVar{i}; %all iterations
    iterationScores = zeros(length(noiseC), 1);
    for iteration = 1:length(noiseC)
        realPhylogeny = realTrees{i}{iteration};

        iterationSamples = noiseSomVar{iteration};
        
        %compute the paiwise distances between all samples
        
        iterationPhasing = noiseC{iteration};
        allSamplePhasing = iterationPhasing;
        distanceMatrix = zeros(size(iterationPhasing,2), size(iterationPhasing,2));
        for sampleInd = 1:length(allSamplePhasing)
            for sampleInd2 = 1:length(allSamplePhasing)
                sample1Phasing = allSamplePhasing{sampleInd}{1}{1};
                sample2Phasing = allSamplePhasing{sampleInd2}{1}{1};
                
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
                distanceMatrix(sampleInd2, sampleInd) = totalDistance;
            end
        end
        distanceMatrix(isinf(distanceMatrix)) = Inf;
        
        %compute distance matrix for somatic variants
        somvarDistanceMatrix = zeros(size(iterationSamples,2), size(iterationSamples,2));
        for sampleInd = 1:size(iterationSamples,2)
            for sampleInd2 = 1:size(iterationSamples,2)
                firstSample = iterationSamples{sampleInd};
                secondSample = iterationSamples{sampleInd2};
                totalDistance = 0;
                for varInd = 1:length(firstSample)
                    if(firstSample(varInd) == 1 && secondSample(varInd) < 1)
                        totalDistance = Inf;
                        break;
                    end
                    distance = abs(firstSample(varInd) - secondSample(varInd));
                    totalDistance = totalDistance + distance;
                end
                somvarDistanceMatrix(sampleInd, sampleInd2) = totalDistance;
            end
        end
        
        %multiply matrices
        totalDistanceMatrix = distanceMatrix .* somvarDistanceMatrix;
        
        testTree = totalDistanceMatrix;
        testTree(isinf(testTree)) = 0;
        if (max(testTree(:)) < 1)
            iterationScore = iterationScore + 0;
            iterationScores(iteration) = score;
            continue; 
        end
        
        bestFit = totalDistanceMatrix;
        bestFit(isinf(bestFit)) = 0;
        bestFit(find(bestFit > 0)) = 1;
        bestFitCell = {bestFit};
        distanceCell = {totalDistanceMatrix};
        bestTree = edmondsPhylogeny(bestFitCell, distanceCell, {});
        
        child_handles = allchild(0);
        names = get(child_handles,'Name');
        k = find(strncmp('Biograph Viewer', names, 15));
        close(child_handles(k))
        if isempty(bestTree)
            iterationScore = iterationScore + 0;
            iterationScores(iteration) = score;
            continue; 
        end
        bestTree = bestTree{1};
        % check how many positions overlap that are not inf
        bestTree(bestTree > 0) = 1;
        if size(bestTree, 1) < size(realPhylogeny)
            iterationScore = iterationScore + 0;
            iterationScores(iteration) = score;
           continue; 
        end
        difference = bestTree - realPhylogeny;
        distance = difference(difference < 0);
        totalDistance = sum(abs(distance));
        score = 1 - (totalDistance / size(realPhylogeny,2));
        iterationScore = iterationScore + score;
        iterationScores(iteration) = score;
    end
    noiseScore{i} = iterationScores;
    %noiseScore(i) = iterationScore / length(noisePhasing);
end

boxplot(cell2mat(noiseScore), noiseLevels, 'colors', [0 0 0])

title('A and somatic variants');
xlabel('Noise (?)');
ylabel('Tree reconstruction accuracy');
figure;
%% C + somatic variants
noiseScore = cell.empty;
for i = 1:length(finalSSC) %for every noise level
    
    noiseC = finalSSC{i}; %all iterations
    
    noiseSomVar = realSomVar{i}; %all iterations
    iterationScores = zeros(length(noiseC), 1);
    for iteration = 1:length(noiseC)
        realPhylogeny = realTrees{i}{iteration};

        iterationSamples = noiseSomVar{iteration};
        
        %compute the paiwise distances between all samples
        
        %compute the distance per position
        %then sum the distances
        iterationC= noiseC{iteration};

        iterationScore = 0;

        distanceMatrix = zeros(size(iterationC,2), size(iterationC,2));
        for sampleInd = 1:size(iterationC,2)
            for sampleInd2 = 1:size(iterationC,2)
                sample1C = iterationC(:, sampleInd);

                sample2C = iterationC(:,sampleInd2);

                %compute sum of pairwise distances
                totalDistance = 0;
                for position = 1:length(sample1C)
                    
                    if sample1C(position) == 0 && sample2C(position) > 0
                        distance = Inf;
                    else
                        distance = abs(sample1C(position) - sample2C(position));
                    end
                    totalDistance = totalDistance + distance;
                end
                distanceMatrix(sampleInd2, sampleInd) = totalDistance;
            end
        end
        distanceMatrix(isinf(distanceMatrix)) = Inf;
        
        %compute distance matrix for somatic variants
        somvarDistanceMatrix = zeros(size(iterationSamples,2), size(iterationSamples,2));
        for sampleInd = 1:size(iterationSamples,2)
            for sampleInd2 = 1:size(iterationSamples,2)
                firstSample = iterationSamples{sampleInd};
                secondSample = iterationSamples{sampleInd2};
                totalDistance = 0;
                for varInd = 1:length(firstSample)
                    if(firstSample(varInd) == 1 && secondSample(varInd) < 1)
                        totalDistance = Inf;
                        break;
                    end
                    distance = abs(firstSample(varInd) - secondSample(varInd));
                    totalDistance = totalDistance + distance;
                end
                somvarDistanceMatrix(sampleInd, sampleInd2) = totalDistance;
            end
        end
        
        %multiply matrices
        totalDistanceMatrix = distanceMatrix .* somvarDistanceMatrix;
        
        testTree = totalDistanceMatrix;
        testTree(isinf(testTree)) = 0;
        if (max(testTree(:)) < 1)
            iterationScore = iterationScore + 0;
            iterationScores(iteration) = score;
            continue; 
        end
        
        bestFit = totalDistanceMatrix;
        bestFit(isinf(bestFit)) = 0;
        bestFit(find(bestFit > 0)) = 1;
        bestFitCell = {bestFit};
        distanceCell = {totalDistanceMatrix};
        bestTree = edmondsPhylogeny(bestFitCell, distanceCell, {});
        
        child_handles = allchild(0);
        names = get(child_handles,'Name');
        k = find(strncmp('Biograph Viewer', names, 15));
        close(child_handles(k))
        bestTree = bestTree{1};
        % check how many positions overlap that are not inf
        bestTree(bestTree > 0) = 1;
        if size(bestTree, 1) < size(realPhylogeny)
            iterationScore = iterationScore + 0;
            iterationScores(iteration) = score;
           continue; 
        end
        difference = bestTree - realPhylogeny;
        distance = difference(difference < 0);
        totalDistance = sum(abs(distance));
        score = 1 - (totalDistance / size(realPhylogeny,2));
        iterationScore = iterationScore + score;
        iterationScores(iteration) = score;
    end
    noiseScore{i} = iterationScores;
end

boxplot(cell2mat(noiseScore), noiseLevels, 'colors', [0 0 0])

title('C and somatic variants');
xlabel('Noise (?)');
ylabel('Tree reconstruction accuracy');



end

