function run_all(outputDirectory)
if nargin < 1
    outputDirectory = fullfile(pwd, 'results');
end
experiments = {'reference', 'baseline', 'initial_data', 'basin_sections', ...
    'accumulation', 'allocation', 'mechanisms', 'sensitivity', ...
    'adaptation', 'matched_data', 'asymmetry'};
for k = 1:numel(experiments)
    run_experiment(experiments{k}, outputDirectory);
end
run_validation(outputDirectory);
end
