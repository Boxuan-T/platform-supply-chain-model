function report = run_validation(outputDirectory)
if nargin < 1
    outputDirectory = fullfile(pwd, 'results');
end
files = dir(fullfile(outputDirectory, '*.mat'));
cutoff = zeros(0, 6);
selected = struct('key', {}, 'p', {}, 'U', {}, 'group', {});
for k = 1:numel(files)
    data = load(fullfile(files(k).folder, files(k).name));
    if ~isfield(data, 'result') || ~isfield(data.result, 'tasks')
        continue;
    end
    result = data.result;
    for j = 1:numel(result.tasks)
        rows = result.outputs{j}.R;
        ids = find(rows(:, 1) <= 0);
        if isempty(ids)
            continue;
        end
        task = result.tasks(j);
        task.key = ['cutoff_', task.key];
        task.U = task.U(:, ids);
        checked = psc.run_case(task, outputDirectory, 'ode15s');
        error = max(abs(rows(ids, 3:6)-checked.R(:, 3:6)) ./ [1, 1, 1, 10], [], 2);
        cutoff = [cutoff; repmat([k, j], numel(ids), 1), ids, rows(ids, 1), checked.R(:, 1), error];
    end
    if ismember(result.experiment, {'reference', 'matched_data'})
        selected = [selected, result.tasks];
    elseif ismember(result.experiment, {'allocation', 'asymmetry'})
        for j = 1:numel(result.maxima)
            keys = [result.maxima(j).basinKeys(1), result.maxima(j).returnKeys(1)];
            if ~isempty(result.maxima(j).marginKeys)
                keys = [keys, result.maxima(j).marginKeys(1)];
            end
            selected = [selected, result.tasks(ismember({result.tasks.key}, keys))];
        end
    end
end
assert(~isempty(files));
assert(all(cutoff(:, 4) == cutoff(:, 5)));
checks = struct([]);
for k = 1:numel(selected)
    task = selected(k);
    coarse = psc.run_case(task, outputDirectory);
    stiff = psc.run_case(task, outputDirectory, 'ode15s');
    assert(all(coarse.R(:, 1) == stiff.R(:, 1)));
    returnError = max(abs(coarse.R(:, 8)-stiff.R(:, 8)));
    assert(returnError < 1e-5);
    fineStats = struct();
    if ~startsWith(task.key, 'matched_data')
        task.key = ['refined_', task.key];
        task.U = psc.strategy_grid(linspace(0.1, 0.9, 9), 4.5);
        fine = psc.run_case(task, outputDirectory);
        fineStats = psc.summarize_case(fine);
    end
    entry = struct('key', selected(k).key, 'returnError', returnError, ...
        'coarse', psc.summarize_case(coarse), 'refined', fineStats);
    checks = [checks, entry];
end
report = struct('files', {{files.name}}, 'cutoff', cutoff, 'checks', checks);
save(fullfile(outputDirectory, 'validation.mat'), 'report');
writetable(array2table(cutoff, 'VariableNames', ...
    {'file_index', 'configuration_index', 'trajectory_index', 'ode45_category', 'ode15s_category', 'terminal_error'}), ...
    fullfile(outputDirectory, 'cutoff_validation.csv'));
end
