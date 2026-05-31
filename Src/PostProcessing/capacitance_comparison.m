% Codice per creare grafico con confronto capacita, con problemi di
% dimensioni diverse.

eps_0 = 8.85e-12;
eps_r = 1;
x = [0.1,0.5,1,1.5];
C_ideal = eps_0*eps_r./x;
C_charge = [];
C_energy = [];

figure;

plot(x, C_ideal, '-o', 'LineWidth', 2, 'MarkerSize', 7); hold on;
plot(x, C_charge, '-s', 'LineWidth', 2, 'MarkerSize', 7);
plot(x, C_energy, '-^', 'LineWidth', 2, 'MarkerSize', 7);

grid on;

xlabel('L_y/L_x');
ylabel('Capacitance (F/m)');
title('Capacitance vs L_y');

legend('C\_ideal','C\_charge','C\_energy','Location','best');