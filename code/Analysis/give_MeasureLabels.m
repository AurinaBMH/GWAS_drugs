function Mlabels = give_MeasureLabels(Mnames)

Mlabels = cell(length(Mnames),1);

for tt=1:length(Mnames)
    
    Mnames{tt} = strrep(Mnames{tt}, '_', ' ');
    
    if contains(Mnames{tt}, 'Allen') && contains(Mnames{tt}, 'Mapped')
        Mlabels{tt} = 'AHBA position  ';
        
    elseif contains(Mnames{tt}, 'Allen') && contains(Mnames{tt}, 'eQTL')
        Mlabels{tt} = 'AHBA eQTL  ';
        
    elseif contains(Mnames{tt}, 'default')
        Mlabels{tt} = 'SNP position  ';
        
    elseif contains(Mnames{tt}, 'eQTL') && contains(Mnames{tt}, '.P') && ~contains(Mnames{tt}, 'Allen')
        
        if contains(Mnames{tt}, 'Blood')
            Mlabels{tt} = 'eQTL Blood  ';
        elseif contains(Mnames{tt}, 'brain')
            Mlabels{tt} = 'eQTL Brain  ';
        elseif contains(Mnames{tt}, 'Terminal')
            Mlabels{tt} = 'eQTL SITI  ';
        elseif contains(Mnames{tt}, 'Omentum')
            Mlabels{tt} = 'eQTL AVO  ';
        elseif contains(Mnames{tt}, 'Subcutaneous')
            Mlabels{tt} = 'eQTL AS  ';
        else
            V = strsplit(Mnames{tt}, '.');
            Vnew = strrep(V{1,1}, 'eQTL', 'eQTL ');
            Mlabels{tt} = strcat(Vnew, "  ");
        end
        
    elseif contains(Mnames{tt}, 'PPI') && contains(Mnames{tt}, 'mapped')
        
        V = strsplit(Mnames{tt}, '.');
        Vnew = strrep(V{1,1},'mapped','position');
        Mlabels{tt} = strcat(Vnew, "  ");
        
    elseif contains(Mnames{tt}, 'PPI') && contains(Mnames{tt}, 'eQTL')
        
        V = strsplit(Mnames{tt}, '.');
        Mlabels{tt} = strcat(V{1,1}, "  ");
    elseif contains(Mnames{tt}, 'Combined')
        Mlabels{tt} = strcat(Mnames{tt}, "  ");
        
    else
        V = strsplit(Mnames{tt}, '.');
        Vnew = insertBefore(V{1,1}, V{1,1},"Chrom ");
        Mlabels{tt} = strcat(Vnew, "  ");
        
    end
end


end
