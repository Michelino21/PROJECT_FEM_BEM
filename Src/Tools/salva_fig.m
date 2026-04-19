function salva_fig(nome, folder_name, SALVATAGGIO)
    if SALVATAGGIO == 1
        saveas(gcf, fullfile(folder_name, 'png', [nome '.png']));
        saveas(gcf, fullfile(folder_name, 'fig', [nome '.fig']));
    end
end