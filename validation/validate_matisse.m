% VALIDATE_MATISSE  Export MATISSE results for the matisseR validation.
%
% Run in MATLAB with the MATISSE toolbox on the path, from matisseR/validation:
%     >> validate_matisse
% Writes validation/matisse/*.csv (no headers, same layout as validate_r.R); then, in R, from the
% package root: Rscript validation/compare.R
%
% Same inputs as the R side: MATISSE example rates (ratesdata, 1 x 1), zero cells filled (fill),
% New APC Model by WLS (apc defaults), CAPRICORN on filled rates.

out = fullfile(pwd, 'matisse');
if ~exist(out, 'dir'), mkdir(out); end
comps = {'lac', 'cac', 'ftt', 'fcp', 'ld'};

%%% APC: estimable functions (log scale: lot + X * B) and key parameters, 16 examples
EF = []; KEY = [];
for i = 1:16
    M = apc(fill(ratesdata(i)));
    for c = 1:numel(comps)
        F = ef(M, 'Comp', comps{c});
        v = F.data.data(:);
        EF = [EF; repmat([i c], numel(v), 1) (1:numel(v))' v]; %#ok<AGROW>
    end
    K = ef(M, 'Comp', 'key');                % Intercept, LAT, NetDrift, CAT, THETAa, THETAp, THETAc
    v = K.data.data(:);
    KEY = [KEY; repmat(i, numel(v), 1) (1:numel(v))' v]; %#ok<AGROW>
end
writematrix(EF, fullfile(out, 'apc_ef.csv'));
writematrix(KEY, fullfile(out, 'apc_key.csv'));

%%% CAPRICORN: {example numbers, cell size}; 5 = chunk to 5 x 5 (LTF)
cmp = {{[9 10], 1}, {[1 2], 1}, {7:10, 5}};
G = []; H = []; C = []; A = []; E = [];
for j = 1:numel(cmp)
    ex = cmp{j}{1}; d = cmp{j}{2};
    R = cell(1, numel(ex));
    for g = 1:numel(ex)
        Rg = ratesdata(ex(g));
        if d > 1, Rg = chunk(Rg, 'AgeBlock', d, 'PerBlock', d); end
        R{g} = fill(Rg);
    end
    S = capricorn(R{:});
    t = S.Tests.PH_Global.data;   G = [G; repmat(j, size(t, 1), 1) (1:size(t, 1))' t];        %#ok<AGROW> Test, df, PVAL
    t = S.Tests.Homogeneity.data; H = [H; repmat(j, size(t, 1), 1) (1:size(t, 1))' t];        %#ok<AGROW> DEV, df, PVAL
    t = S.Tests.Composite.data;   C = [C; repmat(j, size(t, 1), 1) (1:size(t, 1))' t(:, 3)];   %#ok<AGROW> PVAL
    t = S.AIC.data;               A = [A; repmat(j, size(t, 1), 1) (1:size(t, 1))' t];        %#ok<AGROW> mdf, bcdf, DEV, AICc, DELTA, rank, model
    for c = 1:numel(comps)
        [~, D] = ef(S, 'Comp', comps{c});        % differences, strata 1..K-1 vs K
        v = D.data.data;
        for g = 1:size(v, 2)
            E = [E; repmat([j c g], size(v, 1), 1) (1:size(v, 1))' v(:, g)]; %#ok<AGROW>
        end
    end
end
writematrix(G, fullfile(out, 'cap_global.csv'));
writematrix(H, fullfile(out, 'cap_homogeneity.csv'));
writematrix(C, fullfile(out, 'cap_composite.csv'));
writematrix(A, fullfile(out, 'cap_aic.csv'));
writematrix(E, fullfile(out, 'cap_ef.csv'));
disp(['wrote MATISSE results to ' out])
