function [ phylogenies, scores ] = edmondsPhylogeny( uniqueTreeEdges, uniqueTreeWeights, sampleNames )
%EDMONDSPHYLOGENY Applies Edmond's algorithm to infer subclonal evolution
%trees

    %For every input tree, infer the minimum spanning arborescene 
    phylogenies = cell.empty;
    scores = [];
    for treeInd = 1:length(uniqueTreeEdges)
       
        weightedTree = uniqueTreeEdges{treeInd} .* (1 ./ uniqueTreeWeights{treeInd});
        weightedTree(isnan(weightedTree)) = 0;
        
        V = [1:size(weightedTree, 1)];
        E=incidence_to_3n(weightedTree);
        GT= edmonds(V,E);

        % Would need the bioinformatics toolbox for visualization
        bg=biograph(weightedTree, sampleNames);
        % shows the maximum weight graph
        [GTM, TREEMAX]=reconstruct_2(GT, sampleNames);

        phylogenies = [phylogenies full(GTM)];
        currentTree = full(GTM);
        scores = [scores sum(currentTree(:))];
    end
 
end

