function [ PLAFGivenCMu, alleles] = computePLAFGivenCMu( C, mu, LAF, distanceAlleleMap, allAlleleDistanceMatrix)
%COMPUTEPLAFGIVENCMU Computes P(LAF | C, mu)
%   Computes P(LAFi | Ci, mu) given a Ci (at one segment) in one sample and the
%   corresponding LAF measurements. Also infers the corresponding alleles
%   for this Ci. 

    %% Step 1: compute all possible alleles that can be made with this Ci
    %For a copy number of 2 in the tumor the number of 'A' can be 0, 1 or
    %2.
    combvecString = 'combvec(';
    for i = 2:size(C,2)
        combvecString = [combvecString, '0:C(', int2str(1), ',', int2str(i), ')'];

        if i < size(C,2)
            combvecString = [combvecString, ','];
        end
    end
    combvecString = [combvecString, ')'];
    combinations = eval(combvecString);

    %Use the possible numbers of A and B to make all combinations of
    %alleles
    combinationAlleles = cell.empty;
    cloneIndices = 2:size(C,2);
    for combination = 1:size(combinations,2)

        BCounts = (C(1,cloneIndices)' - combinations(:,combination));
        ACounts = C(1,cloneIndices)' - BCounts;

        for subcloneInd = 1:length(BCounts)
            alleles = [repmat('A', 1, ACounts(subcloneInd)), repmat('B', 1, BCounts(subcloneInd))];

            if isempty(alleles)
                combinationAlleles{subcloneInd, combination} = '';
            else
                combinationAlleles{subcloneInd, combination} = alleles;
            end
        end

    end
    alleleCombinations = combinationAlleles;

    %compute all possible LAF for the alleles and the probability density given
    %this Ci and mu
    [mixtureDensity, alleleLocationMap, inferredAlleles] = computeLAFP(alleleCombinations, mu,  distanceAlleleMap, allAlleleDistanceMatrix);

    %for every LAF measurement provided, find the most likely alleles based on the
    %location of the LAF measurement in the probability density
    lafProbabilities = zeros(1, size(LAF,1));
    allelesPerLAF = {};
    for lafInd = 1:length(LAF)
        %find the x of our LAF
        match = find(round((mixtureDensity(:,1)*1000)) / 1000.0 == LAF(lafInd));

        alleles = '?'; %the LAF is not in de density, never occurred in tests
        if(isempty(match))
            lafProbabilities(lafInd) = 0;
            break;
        end

        %Assign the alleles based on the x of our LAF measurement in the
        %density. If there are two options that generate the same LAF measurement
        %(i.e. AA and BB), we select the answer with the most Bs. 
        for i = 1:size(alleleLocationMap,1)
            if match >= alleleLocationMap(i,1) && match <= alleleLocationMap(i,2)
                allelesPerLAF{lafInd} = inferredAlleles{i}(1);
            end
        end
        %obtain the corresponding P(LAF | Ci, mu) for this LAF. 
        prob = mixtureDensity(match,2);
        lafProbabilities(lafInd) = prob;
    end

    %convert the alleles to string format
    stringAlleles = cellfun(@(x) char(x), allelesPerLAF, 'UniformOutput', 0);
    [unique_strings, ~, string_map]=unique(stringAlleles);
    alleles=unique_strings(mode(string_map));
    %if a segment has multiple LAF measurements, we compute the total
    %probability for this Ci and mu by multiplying the scores of all
    %individual measurements. 
    PLAFGivenCMu = prod(lafProbabilities);
    
end

