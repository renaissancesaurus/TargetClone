function [ distanceAlleleMap, allAlleleDistanceMatrix ] = computeAllAlleleDistances( kmin, kmax )
%COMPUTEALLALLELEDISTANCES Computes all possible alleles for a k from kmin
%to kmax and makes an all to all matrix containing event distances for all
%these alleles

    distanceAlleleMap = cell.empty;
    for i = kmin:kmax %for every possible copy number
        %determine the possible alleles
        for j = kmin:i
           ACount = j;
           BCount = i - ACount;
           alleles = [repmat('A', 1, ACount), repmat('B', 1, BCount)];
           if isempty(alleles)
               alleles = {''}; %no alleles, k = 0
           end
           distanceAlleleMap = [distanceAlleleMap alleles];
        end

    end

    allAlleleDistanceMatrix = eventDistance(distanceAlleleMap);

end

