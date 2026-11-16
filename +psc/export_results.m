function result = export_results(tasks, outputs, outputDirectory, experiment)
statistics = cellfun(@psc.summarize_case, outputs);
configurations = cell(numel(tasks), 1);
trajectories = cell(numel(tasks), 1);
for k = 1:numel(tasks)
    parameters = tasks(k).p;
    configuration = struct('key', string(tasks(k).key));
    parameterNames = fieldnames(parameters);
    for j = 1:numel(parameterNames)
        field = parameterNames{j};
        configuration.(field) = parameters.(field);
    end
    names = fieldnames(statistics(k));
    for j = 1:numel(names)
        configuration.(names{j}) = statistics(k).(names{j});
    end
    configurations{k} = struct2table(configuration);
    rows = [tasks(k).U', outputs{k}.R(:, 1:14)];
    trajectories{k} = array2table(rows, 'VariableNames', ...
        {'x0', 'y0', 'z0', 'D0', 'category', 'Tc', 'x', 'y', 'z', 'D', ...
        'alpha', 'meanW', 'meanP', 'meanM', 'meanR', 'residual', 'horizon', 'accounting_error'});
    trajectories{k}.key = repmat(string(tasks(k).key), size(rows, 1), 1);
    trace = outputs{k}.trace;
    if ~isempty(trace.t)
        path = array2table([trace.t, reshape(trace.u, 4, [])'], ...
            'VariableNames', {'time', 'x', 'y', 'z', 'D'});
        writetable(path, fullfile(outputDirectory, [tasks(k).key, '_trace.csv']));
    end
end
configurationTable = vertcat(configurations{:});
trajectoryTable = vertcat(trajectories{:});
maxima = struct([]);
if ismember(experiment, {'allocation', 'asymmetry'})
    groups = unique([tasks.group]);
    for g = groups
        ids = find([tasks.group] == g);
        shares = [[tasks(ids).p]];
        shares = [[shares.sP]', [shares.sM]', [shares.sR]'];
        basin = [statistics(ids).basin]';
        returns = [statistics(ids).meanW]';
        margins = cellfun(@(o) o.cooperation.margin, outputs(ids))';
        stable = cellfun(@(o) o.cooperation.stable, outputs(ids))';
        equal = find(max(abs(shares-1/3), [], 2) < 1e-12);
        basinIds = find(abs(basin-max(basin)) < 1e-12);
        returnIds = find(abs(returns-max(returns)) < 1e-9);
        marginIds = [];
        bestMargin = NaN;
        if any(stable)
            bestMargin = max(margins(stable));
            marginIds = find(stable & abs(margins-bestMargin) < 1e-10);
        end
        entry = struct('group', g, 'maxBasin', max(basin), ...
            'basinShares', shares(basinIds, :), 'basinKeys', {{tasks(ids(basinIds)).key}}, ...
            'maxReturn', max(returns), 'returnShares', shares(returnIds, :), ...
            'returnKeys', {{tasks(ids(returnIds)).key}}, 'maxMargin', bestMargin, ...
            'marginShares', shares(marginIds, :), 'marginKeys', {{tasks(ids(marginIds)).key}}, ...
            'equalBasin', basin(equal), 'equalReturn', returns(equal), ...
            'basinGap', max(basin)-basin(equal), 'returnGap', max(returns)-returns(equal));
        maxima = [maxima, entry];
    end
end
result = struct('experiment', experiment, 'tasks', tasks, 'outputs', {outputs}, ...
    'statistics', statistics, 'maxima', maxima);
save(fullfile(outputDirectory, [experiment, '.mat']), 'result', '-v7.3');
writetable(configurationTable, fullfile(outputDirectory, [experiment, '_configurations.csv']));
writetable(trajectoryTable, fullfile(outputDirectory, [experiment, '_trajectories.csv']));
end
