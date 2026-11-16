function result = run_experiment(experiment, outputDirectory)
if nargin < 1
    experiment = 'baseline';
end
if nargin < 2
    outputDirectory = fullfile(pwd, 'results');
end
tasks = psc.experiment_tasks(experiment);
outputs = cell(size(tasks));
for k = 1:numel(tasks)
    outputs{k} = psc.run_case(tasks(k), outputDirectory);
end
result = psc.export_results(tasks, outputs, outputDirectory, char(experiment));
end
