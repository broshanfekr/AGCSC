load('stl10_results.mat');


[~, idx] = max(results(:, 4));
max_row = results(idx, :);
acc = max_row(3);
nmi = max_row(4);
ari = max_row(5);
fprintf('\nacc = %.2f, nmi = %.2f, ari = %.2f\n', acc*100, nmi*100, ari*100);