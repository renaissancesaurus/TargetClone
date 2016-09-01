%% TargetClone
% TargetClone is a method to infer the most likely copy numbers, alleles
% and frequency of tumor subclones in a sample from lesser allele frequency
% measurements. It can also reconstruct the subclonal evolution tree of a
% tumor. 

% REQUIRED
% 1. At least one sample with one tumor subclone and potential normal cell
%    contamination
% 2. Lesser allele frequency measurements at SNP positions, tab-separated
% 3. Somatic variants (either binary or variant allele frequencies),
%    tab-separated
% 4. A segmentation of the genome that is consistent across all samples
% 5. A file with the corresponding sample names

% USAGE
% 1. Fill in the user parameters
% 2. Run the method from this script


%% user parameters

%required files
lafFile = 'Example/lafMeasurements.txt'; %make an example file!
vafFile = 'Example/somVarVaf.txt';
segmentationFile = 'Example/segmentation.txt';
sampleNamesFile = 'Example/sampleNames.txt';

cOutputFile = 'Output/c.txt';
aOutputFile = 'Output/a.txt';
muOutputFile = 'Output/mu.txt';
treeOutputSuffix = {'Output/', '_subclonalTree'}; %we will prepend the type of combination and automatically save as png

kmin = 0; %minimum allowed copy number
kmax = 6; %maximum allowed copy number
metric = 'median'; %either mean, median or individual
blockSize = 500000; %lower this value when less than 6 GB of RAM is available

%start the method
startTargetClone(lafFile, vafFile, segmentationFile, sampleNamesFile, cOutputFile, aOutputFile, ...
    muOutputFile, treeOutputSuffix, kmin, kmax, metric, blockSize);
