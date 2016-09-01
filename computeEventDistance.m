function [ distanceMatrix ] = computeEventDistance( ACounts, BCounts )
%COMPUTEEVENTDISTANCE Compute the event distance for alleles based on their A and B counts

    distanceMatrix = zeros(size(ACounts,2), size(BCounts,2));
    for i = 1:size(ACounts, 2)
       for j = 1:size(BCounts,2)
           if(BCounts(i) == 0 && BCounts(j) > 0)
               distance = Inf;
           elseif (ACounts(i) == 0 && ACounts(j) > 0)
               distance = Inf;
           else
               alleles1 = [ACounts(i) BCounts(i)];
               alleles2 = [ACounts(j) BCounts(j)];
               distance = sum(abs(alleles1 - alleles2));
           end
            distanceMatrix(i,j) = distance;
       end
        
    end


end

