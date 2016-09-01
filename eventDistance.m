function [ distanceMatrix ] = eventDistance( alleles)
%EVENTDISTANCE Computes the event distance between the input alleles

    distanceMatrix = zeros(size(alleles,2), size(alleles,2));
    for i=1:size(alleles,2)
        for j=1:size(alleles,2)
            %If A or B is lost, it cannot be regained
            if(isempty(cell2mat(strfind(alleles(i), 'B'))) && ~isempty(cell2mat(strfind(alleles(j), 'B'))))
                distance = Inf;
            elseif(isempty(cell2mat(strfind(alleles(i), 'A'))) && ~isempty(cell2mat(strfind(alleles(j), 'A'))))
                distance = Inf;
            elseif(strcmp(alleles(i), 'AABB') && strcmp(alleles(j), 'AB')) % AABB to AB is one event (division)
                distance = 1;
            elseif(strcmp(alleles(i), 'AB') && strcmp(alleles(j), 'AABB'))
                distance = 1;
            else
                distance = strdist(char(alleles(i)), char(alleles(j)), 2);
                distance = distance(2);
            end
            distanceMatrix(i,j) = distance + 1; %prevent division by 0   
        end
    end  
end

