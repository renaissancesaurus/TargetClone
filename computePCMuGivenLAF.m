function [ score, alleles ] = computePCMuGivenLAF( x, mu, LAF, distanceAlleleMap, allAlleleDistanceMatrix)
%COMPUTEPCMUGIVENLAF Computes P(Ci,mu | LAFi) for a single sample and
%segment i

    %We assume that the normal component always has copy number 2
    fullX = [2 x];
  
    %check if we do not have no mu and no copies, then the solution is
    %impossible
    copyTotal = bsxfun(@times, fullX, mu);
    if max(copyTotal) == 0
       score = -Inf;
       likelyPhylogeny = -Inf;
       likelyPhasing = -Inf;
       alleles = '';
       return
    end

    %Compute P(C)
    PC = dist(fullX);
    %Compute P(LAF | C, mu)
    [PLAFGivenCMu, alleles] = computePLAFGivenCMu(fullX, mu, LAF, distanceAlleleMap, allAlleleDistanceMatrix);
    %Compute P(C, mu | LAF)
    score = PLAFGivenCMu * PC;
    
    
end
