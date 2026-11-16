function result = minimal_example()
p = psc.reference_parameters();
[rows, trace] = psc.simulate_batch([0.5; 0.5; 0.5; 4.5], p);
result = struct('parameters', p, 'outcome', rows(1, 1:14), ...
    'trace', trace, 'cooperation', psc.cooperation_metrics(p));
end
