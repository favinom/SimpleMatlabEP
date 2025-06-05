[X,Y,Z] = pg.getCoo;

% Reshape to 2D (X and Y are 2D arrays of size [nx, ny])
% Plot the grid (blue lines)
figure;
hold on;

% Plot horizontal grid lines
for j = 1:size(Y,2)
    plot(X(:,j), Y(:,j), 'b');  % horizontal lines
end

% Plot vertical grid lines
for i = 1:size(X,1)
    plot(X(i,:), Y(i,:), 'b');  % vertical lines
end

% Plot boundary nodes in red
plot(X(boundary_nodes), Y(boundary_nodes), 'ro', 'MarkerSize', 6, 'LineWidth', 1.5);




for i = 1:nk
    plot(el{i}(1), el{i}(2), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    text(el{i}(1), el{i}(2), sprintf('E%d', i), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end
