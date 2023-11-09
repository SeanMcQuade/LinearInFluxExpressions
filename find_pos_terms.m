function [positive_idx,negative_idx] = find_pos_terms(Current_column, All_variables)
%UNTITLED Summary of this function goes here
positive_idx = [];
negative_idx = [];
Current_col_double = double(subs(Current_column,All_variables,ones(1,53)));
for i = 1:length(Current_column)
        if Current_col_double(i) > 0 
            spec_row_idx = i;
            spec_col_idx = j;
            positive_idx = [positive_idx; i];
        elseif Current_col_double(i) < 0 
            spec_row_idx = i;
            spec_col_idx = j;
            negative_idx = [negative_idx; i];
        end
end
end


