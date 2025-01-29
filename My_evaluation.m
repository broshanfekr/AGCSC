%% Image clustering Test
clear, clc

%% Experiment Settings
addpath(genpath(pwd));
% load data
% database = {'ORL','YALEB','COIL20','COIL40','Umist','MNIST'};
database = 'stl10';

numdatas = length(database);

ProjectionType = 0;          % data projection type    
NormalizationType = 0;%2;         % data normalization type               


DataName = database;
%load(database{dataindex});
load(DataName)

%---------------------------------------------------------
FEA = zeros(size(fea, 1), size(fea, 2)*size(fea, 3));
for i = 1:size(fea, 1)
   tmp_img = fea(i, :, :, :);
   tmp_img = squeeze(tmp_img);
   tmp_img = rgb2gray(tmp_img);
   tmp_img = reshape(double(tmp_img), 1, size(fea, 2)*size(fea, 3));
   FEA(i, :) = tmp_img(1,:);
end
fea = FEA;
%---------------------------------------------------------

%X = double(fea);
L = double(gnd);

if min(unique(L)) == 0
    L = L + 1;
end
nbcluster = max(unique(L));


alpha = [1e-5,1e-4,1e-3,5e-3,0.01,0.05,0.1,0.5];
beta  = [1e-5,1e-4,1e-3,5e-3,0.01,0.05,0.1,0.5];
% algorithm
Timecell = cell(length(alpha), length(beta));


breakpoint_path = database + "_results.mat";
if exist(breakpoint_path)
    load(breakpoint_path);
else
    results = [];
end

for a=1:length(alpha)
    for b=1:length(beta)
        % save features and coefficient matrices
        tdir ="./Results/" + DataName;
        tfilename = num2str(a) + "_" + num2str(b) ;

        filename = tdir + "/" + tfilename + ".mat";

        disp(filename + "  " + 'File exists.');
        load(filename)

        C1 = C;
        [idx,~] = clu_ncut(C1,nbcluster);

        [acc, nmi, ari, f1, p, r] = compute_metrics(L, idx);
        results = [results; a, b, acc, nmi, ari];

        save(breakpoint_path, 'results');
    end
end

fprintf('\nacc = %.2f, nmi = %.2f, ari = %.2f\n', acc*100, nmi*100, ari*100)