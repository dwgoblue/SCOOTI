function filtered = filter_gene_list(gene_table, model_genes)
  k = detect_gene_column(gene_table);
  genes = upper(string(gene_table{:,k}));
  filtered = intersect(genes, upper(model_genes));
  filtered = filtered(:);  % Column format
end


function gene_col = detect_gene_column(gene_table)
% Automatically detects whether gene names are in column 1 or 2
% Safely handles single-column tables

    nCols = width(gene_table);

    if nCols == 1
        % Only one column, assume it's the gene column
        gene_col = 1;

    elseif nCols >= 2
        % Check if the column contains mostly text (strings, char, or categorical)
        col1_is_text = iscellstr(gene_table{:,1}) || isstring(gene_table{:,1}) || iscategorical(gene_table{:,1});
        col2_is_text = iscellstr(gene_table{:,2}) || isstring(gene_table{:,2}) || iscategorical(gene_table{:,2});

        if col1_is_text && ~col2_is_text
            gene_col = 1;
        elseif col2_is_text && ~col1_is_text
            gene_col = 2;
        elseif col1_is_text && col2_is_text
            % If both are text, assume col1 is an index
            gene_col = 2;
        else
            error('Could not determine which column contains gene names.');
        end

    else
        error('Table must have at least one column.');
    end
end
