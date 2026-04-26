% --- Errore L2 vs margine ---
    figure('Color','w');
    semilogy([3,20,50], [45.76,3.965,3.444], 'bd-', 'LineWidth', 2, 'MarkerSize', 8)
    xlabel('margin / Ly');
    ylabel('Relative L^2 error on u');
    title('Pure FEM accuracy vs domain size (reference: FEM-BEM)');
    grid on