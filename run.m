function [ bestCSamplesScores, bestCSamples, allSolutionsMu, allSolutionsPhasing ] = run( measurements, segments, allPermutations, allMu, validMuInd, alleleMap, allAlleleDistanceMatrix, method)
%RUN Run TargetClone for validations

    %% set variables
    % set as parameters instead (possible permutations and mu, based on M)

    numSegments = length(segments);
    
    %the actual method

    LAF = cell(1, numSegments);
    for segmentInd = 1:numSegments
        if strcmp(method, 'median') == 1
            LAF{segmentInd} = round((median(measurements{segmentInd})*1000)) / 1000.0;
        elseif strcmp(method, 'mean') == 1
            LAF{segmentInd} = round((mean(measurements{segmentInd})*1000)) / 1000.0;
        else
            LAF{segmentInd} = round((measurements{segmentInd}*1000)) / 1000.0;
        end
    end
    %run the model on these LAF
    
    %determine the number of subclones we should run the model for

    allBest = cell(1, length(validMuInd));
    allBestScores = zeros(1, length(validMuInd));
    allBestPhasings = {};
    all3BestScores = cell.empty;

    parfor muIndex = 1:length(validMuInd)
        
       positionBest = cell(1,3);
        bestPhasing = cell(1,3);
        positionBestScores = zeros(numSegments, 3);
        
        for segmentInd = 1:numSegments
            [allScores, allPhasings] = arrayfun(@(c)computePCMuGivenLAF(c, allMu(validMuInd(muIndex),:)/100, LAF{segmentInd}, alleleMap, allAlleleDistanceMatrix), allPermutations, 'UniformOutput', 0);
            [sortedScores, ind] = sort(cell2mat(allScores),1, 'descend');
            positionBestScores(segmentInd,:) = sortedScores(1:3,1);
            positionBest{1} = [positionBest{1}; allPermutations(ind(1),:)];
            positionBest{2} = [positionBest{2}; allPermutations(ind(2),:)];
            positionBest{3} = [positionBest{3}; allPermutations(ind(3),:)];
            
            bestPhasing{1} = [bestPhasing{1};allPhasings{ind(1)}];
            bestPhasing{2} = [bestPhasing{2}; allPhasings{ind(2)}];
            bestPhasing{3} = [bestPhasing{3}; allPhasings{ind(3)}];

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
    allSolutions{1} = allBest{ind};
    allSolutionsMu{1} = allMu(validMuInd(ind),:);
    allSolutionsPhasing{1} = allBestPhasings{ind};
    bestCSamples{1} = allBest{ind};
    bestCSamplesScores{1} = all3BestScores{ind};

end

