close all;
clear all;

list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

N = [64, 128, 256, 512];

L2SBP42 = [5.2524*10^-5, 6.4869*10^-6, 8.1810*10^-7, 1.0359*10^-7];
L2SBP21 = [6.2767*10^-4, 1.5514*10^-4, 3.8586*10^-5, 9.6211*10^-6];
CSBP21 = [3.0098*10^-3, 7.3571*10^-4, 1.743*10^-4, 4.3337*10^-5];
CSBP42 = [2.2119*10^-4, 2.5079*10^-5, 3.0418*10^-6, 3.706*10^-7];
T2 = [L2SBP21(1) + 10^-4, (L2SBP21(1) + 10^-4)/4, (L2SBP21(1) + 10^-4)/4^2, (L2SBP21(1) + 10^-4)/4^3];
T3 = [L2SBP42(1) + 2*10^-4, (L2SBP42(1) + 2*10^-4)/8, (L2SBP42(1) + 2*10^-4)/8^2, (L2SBP42(1) + 2*10^-4)/8^3];
T3 = T3 / 4;

figure(1)

loglog(N,L2SBP21,'-bo', 'MarkerSize', 4, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'LineWidth', 1.2)
hold on
loglog(N,L2SBP42,'-ro', 'MarkerSize', 4, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'LineWidth', 1.2)

T2 = T2
T3 = T3
loglog(N,T2,'--k')
loglog(N,T3,'-.k')

xticks(N)

legend('sbp21', 'sbp42', '2nd order', '3rd order'); xlabel('$N_c$'); title("$l_2$ error")

set(gca, 'FontSize', 20)

exportgraphics(gcf,'linf_convergence.pdf')