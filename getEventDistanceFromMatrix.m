function [ alleleDistanceMatrix ] = getEventDistanceFromMatrix( alleles, distanceAlleleMap, allAlleleDistanceMatrix )
%GETEVENTDISTANCEFROMMATRIX Provides a distance matrix for all allele cominations and quickly looks up the distance given a set of alleles

    %find the positions of all alleles in the distance matrix
    matchLocations = cellfun(@(allele)strcmp(allele, distanceAlleleMap), alleles, 'UniformOutput', 0);

    %for all given alleles,
    %obtain the distances (all to all)
    alleleDistanceMatrix = zeros(size(alleles,2), size(alleles,2));
    for i=1:size(alleles,2)
        for j=1:size(alleles,2)
            %find the matching allele combination in the map
            alleleDist = allAlleleDistanceMatrix(matchLocations{i}, matchLocations{j});
            alleleDistanceMatrix(i,j) = alleleDist;
        end
    end
    


end

