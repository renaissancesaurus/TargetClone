

%% choose the number of samples to use
kmin = 0;
kmax = 6;
sampleNum = 5;

% decide how many subclones each sample will get (probabilistically)
% the sample order is determined beforehand to know which samples
% have multiple subclones after running the method, otherwise this 
% will become very random. 

%select how many tumor subclones we have per sample, only 1 allowed at the
%moment
possibleSubcloneNum = [1 2]; %each sample will at least have a normal too
cloneNumProb = [(1-2/sampleNum) 2/sampleNum]; %on average we have 2-3 samples with more than 2
cloneNumProb = [1 0];
a = 1:length(cloneNumProb);
sampleSubcloneNum = randsample(a,sampleNum,true,cloneNumProb);
%The number of leaves in the tree is the number of samples, 1 sample is one
%subclone
treeLeafNum = sum(sampleSubcloneNum);

%define possible permutations of C
all2ClonePermutations = permn(kmin+1:kmax, 1);

%define the possible mu to have

all2Mu = permn(1:1:100, 2);
valid2Mu = arrayfun(@(x)isequal (x, 100), sum(all2Mu,2));
valid2MuInd = find(valid2Mu == 1);

%% build the actual tree

% initialize required components

numSegments = 35;
numMeasurements = 500;
noiseLevels = [0, 0.02, 0.04, 0.06, 0.08, 0.1];
somVarNum = 10;

normalGenome = repmat(2, 1, numSegments);

%define the probabilities of changing alleles between kmin and kmax

[alleleMap, allAlleleDistanceMatrix] = computeAllAlleleDistances(kmin, kmax);

%with these distances, we can find the likelihood of switching to a certain
%allele combination given that we started at AB, the normal genome. 
copies = cellfun(@(allele)length(allele), alleleMap);

iterations = 100;

nEvents = 5; %how many events to introduce when a new subclone is formed

solutions = cell(1, iterations);

allClones = cell(1,iterations);
allClonesAlleles = cell(1,iterations);
allClonesLAF = cell(1,iterations);
allClonesMu = cell(1,iterations);
allSolutions = cell(1,iterations);
allSolutionsMu = cell(1,iterations);
allSolutionPhasings = cell(1,iterations);

finalPloidyAccuracy = cell.empty;
finalPurityAccuracy = cell.empty;
finalAlleleAccuracy = cell.empty;
finalPhasingAccuracy = cell.empty;

finalPloidy = cell.empty;
finalPurity = cell.empty;
finalPhasing = cell.empty;
finalGroundTruth = cell.empty;

finalSSPloidyAccuracy = cell.empty;
finalSSAlleleAccuracy = cell.empty;
finalSSPhasingAccuracy = cell.empty;

finalSSPloidy = cell.empty;
finalSSPhasing = cell.empty;

finalSSGroundTruth = cell.empty; %Ambiguity score

finalSomVar = cell.empty;
finalPhylogenies = cell.empty;

profile on

for noiseInd = 1:length(noiseLevels)
    
    noiseLevel = noiseLevels(noiseInd);
    noisePloidyAccuracy = cell.empty;
    noisePurityAccuracy = cell.empty;
    noiseAlleleAccuracy = cell.empty;
    noisePhasingAccuracy = cell.empty;
    
    noiseSSPloidyAccuracy = cell.empty;
    noiseSSAlleleAccuracy = cell.empty;
    noiseSSPhasingAccuracy = cell.empty;
    
    noiseSSPloidySolutions = cell.empty;
    noiseSSAlleleSolutions = cell.empty;
    noiseSSGroundTruth = cell.empty;
    
    noisePuritySolutions = cell.empty;
    noisePloidySolutions = cell.empty;
    noiseAlleleSolutions = cell.empty;
    
    noiseGroundTruth = cell.empty;
    
    noiseSomVar = cell.empty;
    noisePhylogeny = cell.empty;
    for i = 1:iterations
        
        %distribute somatic variants across samples
        n = sampleNum - 1;
        lb = 1;
        ub = somVarNum;
        m=lb:1:ub; %ensure that the segments have a lower bound length of 3
        a=m(sort(randperm(length(m),n)));
        b=diff(a);
        b(end+1)=ub-sum(b);
        somVarDistribution = b;
        
        % make a random subclonal evolution tree, introduce events and
        % somatic variants
        
        currentLeavesSomVar = cell(1,treeLeafNum);
        currentLeavesSomVar{1} = zeros(somVarNum, 1);
        currentLeaves = cell(1,treeLeafNum);
        currentLeaves{1} = normalGenome;
        currentLeavesAlleles = cell(1, treeLeafNum);
        currentLeavesAlleles{1} = repmat({'AB'}, 1, length(normalGenome));
        currentLeafFrequency = zeros(1, treeLeafNum);
        currentLeafFrequency(1) = 1;
        phylogeny = zeros(treeLeafNum, treeLeafNum);
        for leaf = 2:treeLeafNum
            
            %select node to edit based on frequency
            a = 1:length(currentLeafFrequency);
            from = randsample(a,1,true,currentLeafFrequency);
            
            newSubclone = currentLeaves{from};
            cloneAlleles = currentLeavesAlleles{from};
            
            cloneSomVar = currentLeavesSomVar{from};
            somVarEvents = somVarDistribution(leaf-1);
            
            possibleLocations = ones(1, length(cloneSomVar));
            for somVarEvent = 1:somVarEvents
                a = 1:length(possibleLocations);

                R = randsample(a,1,true,possibleLocations);
                cloneSomVar(R) = 1;
                
                possibleLocations = cloneSomVar == 0;
            end
            
            %compute the probabilities given where the subclone started
            for event = 1:nEvents
                modLocation = randperm(length(newSubclone), 1);
                match = cellfun(@(allele)strcmp(allele,currentLeavesAlleles{from}(modLocation)), alleleMap);
                distances = allAlleleDistanceMatrix(match,:);
                editProbabilities = 1./distances;
                %introduce CNV in random places given the parent
                a = 1:length(editProbabilities);             %# possible numbers

                R = randsample(a,1,true,editProbabilities);
                %update the original genome and show the LAF

                newSubclone(modLocation) = copies(R);
                
                cloneAlleles(modLocation) = alleleMap(R);
            end

            currentLeaves{leaf} = newSubclone;
            phylogeny(from, leaf) = 1;
            
            %update the frequency of the leaves
            newFrequency = 0 + (currentLeafFrequency(from)-0).*rand(1);
            currentLeafFrequency(from) = currentLeafFrequency(from) - newFrequency;
            currentLeafFrequency(leaf) = newFrequency;
            currentLeavesAlleles{leaf} = cloneAlleles;
            currentLeavesSomVar{leaf} = cloneSomVar;
        end
        noiseSomVar{i} = currentLeavesSomVar;
        noisePhylogeny{i} = phylogeny;
        %make the phylogenetic tree a proper graph
        dg = sparse(phylogeny);
        ug = tril(dg + dg');
        % then take samples, each sample is 1 subclone
        availableClones = 1:(length(currentLeaves));
        availableMeasurements = numMeasurements; %the measurements to distribute across segments
        defaultAdjacency = zeros(length(currentLeaves), length(currentLeaves));
        
        sampleMeasurements = cell(1, sampleNum);
        sampleGenomes = cell(1, sampleNum);
        sampleGenomesAlleles = cell(1, sampleNum);
        sampleFrequencies = cell(1, sampleNum);
        sampleSegmentation = cell(1, sampleNum);
        %randomly assign a number of measurements to each segment
        n = numSegments;
        lb = 1;
        ub = availableMeasurements;
        m=lb:3:ub; %ensure that the segments have a lower bound length of 3
        a=m(sort(randperm(length(m),n)));
        b=diff(a);
        b(end+1)=ub-sum(b);
        segmentDistribution = b;
        
        for sampleInd = 1:sampleNum
            
            cloneNum = sampleSubcloneNum(sampleInd);
            clones = cell(1, cloneNum);
            
            %randomly select first clone to sample
            randIndex = randi(length(availableClones), 1);
            randIndexAdjacency = defaultAdjacency;
            randIndexAdjacency(randIndex,randIndex) = 1;
            clones{1} = availableClones(randIndex);
           
            %now we have selected our clones for the sample.
            
            %determine the mu for the sample.
            n = cloneNum + 1;
            lb = 0;
            ub = 100;
            m=lb:ub;
            a=m(sort(randperm(length(lb:ub),n)));
            b=diff(a);
            b(end+1)=ub-sum(b);
            mu = b / 100;
            
            %always add the normal genome (sample 1) to the sample
            clones = [1 clones];
            
            ACounts = zeros(1, length(normalGenome));
            BCounts = zeros(1, length(normalGenome));
            totalGenome = cell(1, length(clones));
            totalGenomeAlleles = cell(1, length(clones));
            
            %assign measurements to the segments
            %we decide beforehand what the CN of the segment will be, then
            %we choose the possible LAF that could have been measured.
            %the possible LAF distribution needs to be made dynamically for
            %increasing noise levels.
            
            
            %compute the LAF for the sample
            for cloneInd = 1:length(clones)
                clone = clones{cloneInd};
                totalGenomeAlleles{cloneInd} = currentLeavesAlleles{clone};
                totalGenome{cloneInd} = currentLeaves{clone};
                
                ALocations = cellfun(@(alleles)find(alleles == 'A'), currentLeavesAlleles{clone}, 'UniformOutput', 0);
                ACount = cellfun(@length, ALocations);
                
                muACounts = bsxfun(@times, ACount, mu(cloneInd));
                ACounts = ACounts + muACounts;
                
                BLocations = cellfun(@(alleles)find(alleles == 'B'), currentLeavesAlleles{clone}, 'UniformOutput', 0);
                BCount = cellfun(@(count)length(count), BLocations);
                
                mBCounts = bsxfun(@times, BCount, mu(cloneInd));
                BCounts = BCounts + mBCounts;
                
            end
            
            segmentLAF = min(ACounts, BCounts) ./ (ACounts + BCounts);
            
            
            %given these numbers of measurements per segment, define
            %the actual LAF measurements
            
            %the LAF can be sampled from the original + noise factor
            
            allMeasurements = cell(1,numSegments);
            %introduce multiple measurements per sample
            for segmentInd = 1:numSegments
                
                lafDistribution = makedist('Normal', 'mu', segmentLAF(segmentInd), 'sigma', noiseLevel);
                if noiseLevel ~= 0
                    lafDistribution = truncate(lafDistribution,0,0.5);
                end
                measurements = random(lafDistribution,segmentDistribution(segmentInd),1);
                allMeasurements{segmentInd} = measurements;
                
            end
            %Update and remove samples that have already been taken from
            %the available sample list
            sampleMeasurements{sampleInd} = allMeasurements;
            sampleSegmentation{sampleInd} = segmentDistribution;
            ind = cellfun(@(x)find(x == availableClones), clones(2:length(clones)), 'UniformOutput', 0);
            availableClones(cell2mat(ind)) = [];
            sampleGenomes{sampleInd} = totalGenome;
            sampleGenomesAlleles{sampleInd} = totalGenomeAlleles;
            sampleFrequencies{sampleInd} = mu;
            
        end
        
        %for these samples, we now wish to run targetclone and validate the
        %outcome.
        
        ploidyAccuracy = 0;
        purityAccuracy = 0;
        alleleAccuracy = 0;
        phasingAccuracy = 0;
        samplePhasings = cell.empty;
        samplePloidy = cell.empty;
        samplePurity = cell.empty;
        sampleScores = cell.empty;
        for sample = 1:sampleNum
            
            if sampleSubcloneNum(sample) == 1
                validMuInd = valid2MuInd;
                allMu = all2Mu;
                allPermutations = all2ClonePermutations;
            elseif sampleSubcloneNum(sample) == 2
                validMuInd = valid3MuInd;
                allMu = all3Mu;
                allPermutations = all3ClonePermutations;
            end
            
            [score, c, mu, alleles] = run(sampleMeasurements{sample}, sampleSegmentation{sample}, allPermutations, allMu, validMuInd, alleleMap, allAlleleDistanceMatrix, 'individual');
            
            samplePhasings{sample} = alleles;
            samplePloidy{sample} = c;
            samplePurity{sample} = mu;
            sampleScores{sample} = score;
            %Compute accuracy for mu
            relativePurityDistance = sum(abs(cell2mat(mu)/100 - sampleFrequencies{sample})) / size(mu{1},2);
            purityAccuracy = (1-relativePurityDistance);
        end
        
        %single-sample accuracy
        genomes = cellfun(@(x) x{2}, sampleGenomes, 'UniformOutput', 0);
        
        allGenomes = zeros(numSegments, length(sampleGenomes));
        for sample = 1:length(genomes)
           
            allGenomes(:,sample) = genomes{sample}';
            
        end
        
        ssSolutions = cell2mat(cellfun(@(x) x{1}{1}, samplePloidy, 'UniformOutput', 0));
        ssAlleleSolutionsPerSample = cellfun(@(x) x{1}{1}, samplePhasings, 'UniformOutput', 0);
        
        ssAlleleSolutionsPerSegment = cell.empty;
        
        for ssSegment = 1:length(ssAlleleSolutionsPerSample{1})
            tmp = cell.empty;
            for ssSample = 1:length(ssAlleleSolutionsPerSample)
                tmp = [tmp ssAlleleSolutionsPerSample{ssSample}{ssSegment}];
                
            end
            ssAlleleSolutionsPerSegment{ssSegment} = tmp;
        end
        
        
        groundTruthMatches = 0;
        ploidyDistance = abs(allGenomes - ssSolutions);
        ploidyDistance = ploidyDistance ~= 0;
        ploidySSAccuracy = (1 - sum(ploidyDistance(:)) / (numSegments * length(sampleGenomes)));
        alleleMismatches = 0;
        phasingMatches = 0;
        %how well do the alleles match?
        for alleleSegment = 1:length(ssAlleleSolutionsPerSegment)
            
            originalAlleles = cellfun(@(x) x{2}(alleleSegment), sampleGenomesAlleles, 'UniformOutput', 0);
            originalAlleleMatrix = cellfun(@(x) x, originalAlleles);
            %compute the distance between the two
            
            AMatchesTrue = strfind(originalAlleleMatrix, 'A');
            BMatchesTrue = strfind(originalAlleleMatrix, 'B');
            
            ALengthsTrue = cellfun(@length, AMatchesTrue);
            BLengthsTrue = cellfun(@length, BMatchesTrue);

            trueLAF = min(ALengthsTrue, BLengthsTrue) ./ (ALengthsTrue + BLengthsTrue);
            trueLAF(isnan(trueLAF)) = 0;
            
            AMatchesPred = strfind(ssAlleleSolutionsPerSegment{alleleSegment}, 'A');
            BMatchesPred = strfind(ssAlleleSolutionsPerSegment{alleleSegment}, 'B');
            
            ALengthsPred = cellfun(@length, AMatchesPred);
            BLengthsPred = cellfun(@length, BMatchesPred);
            
            predLAF = min(ALengthsPred, BLengthsPred) ./ (ALengthsPred + BLengthsPred);
            predLAF(isnan(predLAF)) = 0;
            
            mismatches = sum(trueLAF ~= predLAF);
            lafMatches = trueLAF == predLAF;
            alleleMismatches = alleleMismatches + mismatches;
            
            %also compute the mismatches for the ploidy, where the ploidy
            %does not match but the LAF does match increases the ambiguity
            %score
            solutionGenome = ssSolutions(alleleSegment,:);
            realGenome = allGenomes(alleleSegment,:);
            matches = solutionGenome ~= realGenome;
            ambiguities = lafMatches == 1 & matches == 1;
            groundTruthMatches = groundTruthMatches + sum(ambiguities);
            
            % compute the allele accuracy. If ploidy is same & LAF, then true
            matches = solutionGenome == realGenome;
            ambiguities = lafMatches == 1 & matches == 1;
            phasingMatches = phasingMatches + sum(ambiguities);

        end
        ssAlleleAccuracy = 1 - (mismatches / (numSegments * length(sampleGenomes)));
        ssGroundTruthScore = 1 - (groundTruthMatches / (numSegments * length(sampleGenomes)));
        ssPhasingAccuracy = (phasingMatches / (numSegments * length(sampleGenomes)));
        noiseSSPloidyAccuracy{i} = ploidySSAccuracy;
        noiseSSAlleleAccuracy{i} = ssAlleleAccuracy;
        noiseSSPhasingAccuracy{i} = ssPhasingAccuracy;

        noiseSSPloidySolutions{i} = ssSolutions;
        noiseSSAlleleSolutions{i} = samplePhasings;
        noiseSSGroundTruth{i} = ssGroundTruthScore;
        
        %% Include information across samples and re-compute the accuracy
        
        %Make all combinations, works fine without blocks with 5 samples
        %and 35 segments
        solutions = zeros(numSegments, length(samplePloidy));
        alleleSolutions = cell(numSegments, 1);
        for segmentInd = 1:numSegments

            combvecString = 'combvec(';

            for sample = 1:length(samplePloidy)

                %obtain the 1,2 best for this segment
                positionBest = samplePloidy{sample}{1};

                combvecString = [combvecString, '[', int2str(positionBest{1}(segmentInd)), ', ', int2str(positionBest{2}(segmentInd)), ']'];

                if(sample < length(samplePloidy))

                    combvecString = [combvecString, ', '];
                end

            end
            combvecString = [combvecString, ')'];

            combinations = eval(combvecString);

             normalClone = 2;

            %score for this segment: associate the segment with the correct
            %combination in C
            combinationAlleles = cell(1, length(combinations));
            combinationScores = zeros(length(combinations),1 );
            parfor combination = 1:length(combinations)
                LAFPs = zeros(1, length(samplePloidy));
                allAlleles = [];
                %obtain the LAFP for this combination and sample
                for sample = 1:length(samplePloidy) %ignore the first

                    %what is the score for this sample?
                    cn = combinations(sample,combination);
                    segmentC = cellfun(@(x) x(segmentInd), samplePloidy{sample}{1}, 'UniformOutput', 0);
                    correctInd = cell2mat(segmentC) == cn;
                    segmentLAFP = sampleScores{sample}{1}(segmentInd, correctInd);
                    LAFPs(sample) = segmentLAFP;
                    correspondingAlleles = samplePhasings{sample}{1}{correctInd}(segmentInd);
                    allAlleles = [allAlleles correspondingAlleles];
                end
                fullCombination = [normalClone; combinations(:, combination)];
                fullAlleleCombination = ['AB' allAlleles];
                combinationAlleles{combination} = allAlleles;
                %compute the total distance for the alleles

                %compute the overall distance
                distance = 1 / (sum(pdist(fullCombination, 'cityblock')) +1);
                %compute the score of the solution for a sample
                combinationScore = LAFPs * distance;
                combinationScores(combination) = prod(combinationScore);
            end
            %find the best CN of this segment for each sample
            [~, ind] = max(combinationScores);
            bestCombination = combinations(:,ind);
            bestCombinationAlleles = combinationAlleles{ind};
            %these are the solutions for all samples.
            solutions(segmentInd, :) = bestCombination;
            alleleSolutions{segmentInd} = bestCombinationAlleles;
        end

        %compute the accuracy
        %for all samples, how well do the copy numbers match?
        
        %collect all sample genomes in a matrix
        genomes = cellfun(@(x) x{2}, sampleGenomes, 'UniformOutput', 0);
        
        allGenomes = zeros(numSegments, length(sampleGenomes));
        for sample = 1:length(genomes)
           
            allGenomes(:,sample) = genomes{sample}';
            
        end
        groundTruthMatches = 0;
        ploidyDistance = abs(allGenomes - solutions);
        ploidyDistance = ploidyDistance ~= 0;
        ploidyAccuracy = (1 - sum(ploidyDistance(:)) / (numSegments * length(sampleGenomes)));
        alleleMismatches = 0;
        phasingMatches = 0;
        %how well do the alleles match?
        for alleleSegment = 1:length(alleleSolutions)
            
            originalAlleles = cellfun(@(x) x{2}(alleleSegment), sampleGenomesAlleles, 'UniformOutput', 0);
            originalAlleleMatrix = cellfun(@(x) x, originalAlleles);
            %compute the distance between the two
            
            AMatchesTrue = strfind(originalAlleleMatrix, 'A');
            BMatchesTrue = strfind(originalAlleleMatrix, 'B');
            
            ALengthsTrue = cellfun(@length, AMatchesTrue);
            BLengthsTrue = cellfun(@length, BMatchesTrue);

            trueLAF = min(ALengthsTrue, BLengthsTrue) ./ (ALengthsTrue + BLengthsTrue);
            trueLAF(isnan(trueLAF)) = 0;
            
            AMatchesPred = strfind(alleleSolutions{alleleSegment}, 'A');
            BMatchesPred = strfind(alleleSolutions{alleleSegment}, 'B');
            
            ALengthsPred = cellfun(@length, AMatchesPred);
            BLengthsPred = cellfun(@length, BMatchesPred);
            
            predLAF = min(ALengthsPred, BLengthsPred) ./ (ALengthsPred + BLengthsPred);
            predLAF(isnan(predLAF)) = 0;
            
            mismatches = sum(trueLAF ~= predLAF);
            lafMatches = trueLAF == predLAF;
            alleleMismatches = alleleMismatches + mismatches;
            
            %also compute the mismatches for the ploidy
            solutionGenome = solutions(alleleSegment,:);
            realGenome = allGenomes(alleleSegment,:);
            matches = solutionGenome ~= realGenome;
            ambiguities = lafMatches == 1 & matches == 1;
            groundTruthMatches = groundTruthMatches + sum(ambiguities);
            
            % compute the allele accuracy. If ploidy is same & LAF, then true
            matches = solutionGenome == realGenome;
            ambiguities = lafMatches == 1 & matches == 1;
            phasingMatches = phasingMatches + sum(ambiguities);
            
            
        end
        alleleAccuracy = 1 - (mismatches / (numSegments * length(sampleGenomes)));
        groundTruthScore = 1 - (groundTruthMatches / (numSegments * length(sampleGenomes)));
        phasingAccuracy = (phasingMatches / (numSegments * length(sampleGenomes)));
        
        noisePloidyAccuracy{i} = ploidyAccuracy;
        noisePurityAccuracy{i} = purityAccuracy;
        noiseAlleleAccuracy{i} = alleleAccuracy;
        noiseGroundTruth{i} = groundTruthScore;
        noisePhasingAccuracy{i} = phasingAccuracy;
        
        noisePuritySolutions{i} = samplePurity;
        noisePloidySolutions{i} = solutions;
        noiseAlleleSolutions{i} = alleleSolutions;
        
    end

    finalPloidyAccuracy{noiseInd} = noisePloidyAccuracy;
    finalPurityAccuracy{noiseInd} = noisePurityAccuracy;
    finalAlleleAccuracy{noiseInd} = noiseAlleleAccuracy;
    finalPhasingAccuracy{noiseInd} = noisePhasingAccuracy;

    finalPloidy{noiseInd} = noisePloidySolutions;
    finalPurity{noiseInd} = noisePuritySolutions;
    finalPhasing{noiseInd} = noiseAlleleSolutions;
    
    finalPhylogenies{noiseInd} = noisePhylogeny;
    finalSomVar{noiseInd} = noiseSomVar;
    finalGroundTruth{noiseInd} = noiseGroundTruth;

    finalSSPloidyAccuracy{noiseInd} = noiseSSPloidyAccuracy;
    finalSSAlleleAccuracy{noiseInd} = noiseSSAlleleAccuracy;
    finalSSPhasingAccuracy{noiseInd} = noiseSSPhasingAccuracy;

    finalSSPloidy{noiseInd} = noiseSSPloidySolutions;
    finalSSPhasing{noiseInd} = noiseSSAlleleSolutions;

    finalSSGroundTruth{noiseInd} = noiseSSGroundTruth;
end

%% process the results and make graphs for the results when the noise increases

% removed the updating because it does not seem necessary anymore?
updatedPloidyAccuracy = cellfun(@(x)(cell2mat(x)), finalPloidyAccuracy, 'UniformOutput', 0);

alleleScoreMeans = cellfun(@(x)mean(cell2mat(x)), finalAlleleAccuracy, 'UniformOutput', 0);
alleleScoreSd = cellfun(@(x)std(cell2mat(x)), finalAlleleAccuracy, 'UniformOutput', 0);
purityScoreMeans = cellfun(@(x)mean(cell2mat(x)), finalPurityAccuracy, 'UniformOutput', 0);
purityScoreSd = cellfun(@(x)std(cell2mat(x)), finalPurityAccuracy, 'UniformOutput', 0);
ploidyScoreMeans = cellfun(@(x)mean(x), updatedPloidyAccuracy, 'UniformOutput', 0);
ploidyScoreSd = cellfun(@(x)std(x), updatedPloidyAccuracy, 'UniformOutput', 0);
phasingScoreMeans = cellfun(@(x)mean(cell2mat(x)), finalPhasingAccuracy, 'UniformOutput', 0);
phasingScoreSd = cellfun(@(x)std(cell2mat(x)), finalPhasingAccuracy, 'UniformOutput', 0);


alleleScores = cellfun(@(x)(cell2mat(x)), finalPhasingAccuracy, 'UniformOutput', 0);
alleleScoresMat = zeros(length(alleleScores{1}), length(alleleScores)*2);

counter = 1;
for i = 1:2:(length(alleleScores)*2)
   alleleScoresMat(:,i) = cell2mat(finalGroundTruth{counter});
   alleleScoresMat(:,i+1) = alleleScores{counter};
   counter = counter + 1;
end

figure('color',[1,1,1]);
boxplot(alleleScoresMat,'factorgap',10,'color','bk')

set(gca,'xtick',1.8:4.3:50)
set(gca,'xticklabel',noiseLevels)

xlabel('Noise (?)');
ylabel('A accuracy');
figure;

finalGroundTruthMeans = cellfun(@(x)mean(cell2mat(x)), finalGroundTruth, 'UniformOutput', 0);
finalGroundTruthMat = cellfun(@(x)cell2mat(x), finalGroundTruth, 'UniformOutput', 0);

groundTruthScoresMat = zeros(length(alleleScores{1}), length(alleleScores));

for i = 1:length(alleleScores)
   groundTruthScoresMat(:,i) = finalGroundTruthMat{i};
end

alleleScores = cellfun(@(x)(cell2mat(x)), finalAlleleAccuracy, 'UniformOutput', 0);
alleleScoresMat = zeros(length(alleleScores{1}), length(alleleScores));

for i = 1:length(alleleScores)
   alleleScoresMat(:,i) = alleleScores{i};
end

boxplot(alleleScoresMat, noiseLevels, 'colors', [0 0 0])

hold on
plot(noiseLevels, cell2mat(finalGroundTruthMeans), 'b*')
xlabel('Noise (?)');
ylabel('LAF accuracy');
figure;

alleleScores = cellfun(@(x)(cell2mat(x)), finalPurityAccuracy, 'UniformOutput', 0);
alleleScoresMat = zeros(length(alleleScores{1}), length(alleleScores));

for i = 1:length(alleleScores)
   alleleScoresMat(:,i) = alleleScores{i};
end

boxplot(alleleScoresMat, noiseLevels, 'colors', [0 0 0])
xlabel('Noise (?)');
ylabel('µ accuracy');

hold on
plot(noiseLevels, cell2mat(finalGroundTruthMeans), 'b*')
figure;


alleleScoresMat = zeros(length(updatedPloidyAccuracy{1}), length(updatedPloidyAccuracy));
counter = 1;
for i = 1:2:(length(alleleScores)*2)
   alleleScoresMat(:,i) = cell2mat(finalGroundTruth{counter});
   alleleScoresMat(:,i+1) = updatedPloidyAccuracy{counter};
   counter = counter + 1;
end

figure('color',[1,1,1]);
boxplot(alleleScoresMat,'factorgap',10,'color','bk')

set(gca,'xtick',1.8:4.3:50)
set(gca,'xticklabel',noiseLevels)

xlabel('Noise (?)');
ylabel('C accuracy');
figure;

%%how often is comb 1 chosen and how often comb 2?

comb1Scores = cellfun(@(x)(cell2mat(x)), finalComb1Chosen, 'UniformOutput', 0);
comb1ScoresMat = zeros(length(finalComb1Chosen{1}), length(finalComb1Chosen));

for i = 1:length(finalComb1Chosen)
   comb1ScoresMat(:,i) = comb1Scores{i};
end

boxplot(comb1ScoresMat, noiseLevels,'positions',noiseLevels, 'colors', [0 0 0])
set(gca, 'YLim', [0,1]);

xlabel('Noise (?)');
ylabel('Ratio of choosing best copy number');
figure;

comb2Scores = cellfun(@(x)(cell2mat(x)), finalComb2Chosen, 'UniformOutput', 0);
comb2ScoresMat = zeros(length(finalComb2Chosen{1}), length(finalComb2Chosen));

for i = 1:length(finalComb2Chosen)
   comb2ScoresMat(:,i) = comb2Scores{i};
end

boxplot(comb2ScoresMat, noiseLevels,'positions',noiseLevels, 'colors', [0 0 0])
set(gca, 'YLim', [0,1]);


xlabel('Noise (?)');
ylabel('Ratio of choosing second best copy number');

%% repeat for SS results

alleleScoreMeans = cellfun(@(x)mean(cell2mat(x)), finalSSAlleleAccuracy, 'UniformOutput', 0);
alleleScoreSd = cellfun(@(x)std(cell2mat(x)), finalSSAlleleAccuracy, 'UniformOutput', 0);
ploidyScoreMeans = cellfun(@(x)mean(cell2mat(x)), finalSSPloidyAccuracy, 'UniformOutput', 0);
ploidyScoreSd = cellfun(@(x)std(cell2mat(x)), finalSSPloidyAccuracy, 'UniformOutput', 0);
phasingScoreMeans = cellfun(@(x)mean(cell2mat(x)), finalSSPhasingAccuracy, 'UniformOutput', 0);
phasingScoreSd = cellfun(@(x)std(cell2mat(x)), finalSSPhasingAccuracy, 'UniformOutput', 0);


finalSSGroundTruthMeans = cellfun(@(x)mean(cell2mat(x)), finalSSGroundTruth, 'UniformOutput', 0);

alleleScores = cellfun(@(x)cell2mat(x), finalSSPhasingAccuracy, 'UniformOutput', 0);
alleleScoresMat = zeros(length(alleleScores{1}), length(alleleScores)*2);

counter = 1;
for i = 1:2:(length(alleleScores)*2)
   alleleScoresMat(:,i) = cell2mat(finalGroundTruth{counter});
   alleleScoresMat(:,i+1) = alleleScores{counter};
   counter = counter + 1;
end

figure('color',[1,1,1]);
boxplot(alleleScoresMat,'factorgap',10,'color','bk')

set(gca,'xtick',1.8:4.3:50)
set(gca,'xticklabel',noiseLevels)

xlabel('Noise (?)');
ylabel('A accuracy');
figure;

alleleScores = cellfun(@(x)(cell2mat(x)), finalSSAlleleAccuracy, 'UniformOutput', 0);
alleleScoresMat = zeros(length(alleleScores{1}), length(alleleScores));

for i = 1:length(alleleScores)
   alleleScoresMat(:,i) = alleleScores{i};
end

boxplot(alleleScoresMat, noiseLevels, 'colors', [0 0 0])

hold on
plot(noiseLevels, cell2mat(finalSSGroundTruthMeans), 'b*')
xlabel('Noise (?)');
ylabel('LAF accuracy');
figure;

alleleScores = cellfun(@(x)(cell2mat(x)), finalSSPloidyAccuracy, 'UniformOutput', 0);
alleleScoresMat = zeros(length(alleleScores{1}), length(alleleScores)*2);

counter = 1;
for i = 1:2:(length(alleleScores)*2)
   alleleScoresMat(:,i) = cell2mat(finalGroundTruth{counter});
   alleleScoresMat(:,i+1) = alleleScores{counter};
   counter = counter + 1;
end

figure('color',[1,1,1]);
boxplot(alleleScoresMat,'factorgap',10,'color','bk')

set(gca,'xtick',1.8:4.3:50)
set(gca,'xticklabel',noiseLevels)

xlabel('Noise (?)');
ylabel('C accuracy');


%% validate the trees as well
validateTrees(finalPloidy, finalPhasing, finalPhylogenies, finalSomVar, finalSSPloidy, finalSSPhasing);
