function [pos_rows_in_col] = find_pos_rows_in_column(Current_col)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
positive_idx = [];
negative_idx = [];
for j = 1:length(Columns_to_change)
    for i = 1:length(C_0)
        if S1_double_trans(i,Columns_to_change(j)) > 0 
            spec_row_idx = i;
            spec_col_idx = j;
            positive_idx = [positive_idx; i,Columns_to_change(j)];
        elseif S1_double_trans(i,Columns_to_change(j)) < 0 
            spec_row_idx = i;
            spec_col_idx = j;
            negative_idx = [negative_idx; i,Columns_to_change(j)];
        end
    end
end

positive_idx
negative_idx
end

