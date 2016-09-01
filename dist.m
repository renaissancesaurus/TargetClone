function [ distance ] = dist( copyNumberProfile)
%DIST Compute the manhattan distance between all segments across samples

    distance = 1 / (sum(pdist(copyNumberProfile', 'cityblock')) +1);

end

