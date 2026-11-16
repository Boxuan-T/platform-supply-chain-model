function U = strategy_grid(g, D)
[x, y, z] = ndgrid(g, g, g);
U = [x(:), y(:), z(:), zeros(numel(x), 1)+D]';
end
