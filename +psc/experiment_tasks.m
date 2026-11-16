function tasks = experiment_tasks(experiment)
experiment = validatestring(experiment, {'reference', 'baseline', 'initial_data', ...
    'basin_sections', 'accumulation', 'allocation', 'mechanisms', ...
    'sensitivity', 'adaptation', 'matched_data', 'asymmetry'});
p = psc.reference_parameters();
U = psc.strategy_grid([0.1, 0.3, 0.5, 0.7, 0.9], 4.5);
tasks = struct('key', {}, 'p', {}, 'U', {}, 'group', {});
switch experiment
    case 'reference'
        add(p, U, 1);
    case 'baseline'
        for k = 1:3
            q = p;
            if k == 2
                q.delta = 1.5;
            elseif k == 3
                q.sP = 0.9;
                q.sM = 0.05;
                q.sR = 0.05;
            end
            add(q, [0.5; 0.5; 0.5; 4.5], k);
        end
    case 'initial_data'
        for d = [0, 0.25, 0.5, 1, 2]
            add(p, psc.strategy_grid(0.05 + 0.09 * (0:10), 9 * d), d);
        end
    case 'basin_sections'
        for z = [0.3, 0.7]
            for d = [0, 0.5, 1]
                [x, y] = ndgrid(linspace(0.05, 0.95, 51));
                states = [x(:), y(:), repmat(z, numel(x), 1), ...
                    repmat(9 * d, numel(x), 1)]';
                add(p, states, [z, d]);
            end
        end
    case 'accumulation'
        for a = linspace(0.7, 2, 21)
            for d = linspace(0.2, 1.6, 21)
                q = p;
                q.aH = a;
                q.delta = d;
                add(q, U, 1);
            end
        end
    case {'allocation', 'asymmetry'}
        scenarios = 1;
        if strcmp(experiment, 'asymmetry')
            scenarios = 1:3;
        end
        for scenario = scenarios
            q = p;
            if scenario == 2
                q.cM = 1.5;
            elseif scenario == 3
                q.qR = 1.5;
            end
            for i = 1:30
                for j = 0:30-i
                    q.sP = i / 30;
                    q.sM = j / 30;
                    q.sR = (30-i-j) / 30;
                    add(q, U, scenario);
                end
            end
        end
    case 'mechanisms'
        for k = 1:6
            q = p;
            switch k
                case 2
                    q.mode = 'fixed';
                case 3
                    q.mode = 'instant';
                case 4
                    q.qM = 1.5;
                    q.qR = 1.5;
                    q.qC = 0;
                case 5
                    q.revenue = 'exp';
                case 6
                    q.revenue = 'hill';
            end
            add(q, U, k);
        end
    case 'sensitivity'
        fields = {'qM', 'qR', 'qC', 'v', 'K', 'kH', 'CH', 'cM', 'cR'};
        ranges = [0.5, 2; 0.5, 2; 0, 2; 8, 16; 2, 8; 1.2, 2; 1.5, 3; 0.5, 2; 0.5, 2];
        for k = 1:numel(fields)
            for value = linspace(ranges(k, 1), ranges(k, 2), 9)
                q = p;
                q.(fields{k}) = value;
                add(q, U, k);
            end
        end
    case 'adaptation'
        fields = {'nuP', 'nuM', 'nuR'};
        for k = 1:3
            for value = [0.5, 2]
                q = p;
                q.(fields{k}) = value;
                add(q, U, k);
            end
        end
    case 'matched_data'
        U(4, :) = (p.aL + (p.aH-p.aL)*U(1, :)) .* ...
            (p.qM*U(2, :) + p.qR*U(3, :) + p.qC*U(2, :).*U(3, :)) / p.delta;
        add(p, U, 1);
        p.mode = 'instant';
        add(p, U, 2);
end

    function add(parameters, states, group)
        key = sprintf('%s_%04d', experiment, numel(tasks) + 1);
        tasks(end + 1) = struct('key', key, 'p', parameters, 'U', states, 'group', group);
    end
end
