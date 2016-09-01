function [ normalizedDensity, locs ] = generateGaussianMixtureModel( mu, sigma, ratios )
%GENERATEGAUSSIANMIXTUREMODEL Generates the mixture model P(LAFi | Ci, mu)

	%this is the sigma that we estimate in the mixture model when building a mixture model with sigma 0.02. 
    sigma = 5.3619e-04; 

    fit = gmdistribution(mu', sigma, ratios);

    xx = 0:0.001:1;
    density = pdf(fit, xx');
    %scale to correct 0-1 range
    x = 0;
    y = max(ratios);
    minDensity = min(density);
     range = max(density) - minDensity;
     density = (density - minDensity) / range;

     range2 = y - x;
     normalized = (density*range2) + x;
    normalizedDensity = [xx', normalized];
    
    %find the start and end of each part of the density (valleys)
    %this to associate a LAF with the most likely phasing
   
    
    locs = find(diff(diff(-normalizedDensity(:,2))>0)<0);
end

