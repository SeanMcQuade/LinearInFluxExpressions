function [positive_idx,negative_idx] = find_pos_neg_rows_in_column(S,S1_double,Columns_to_change)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
num_x = size(S1_double,1);
num_f = size(S1_double,2);
V_0 = eye(num_f);
C_0 = V_0*transpose(S);
S1_double_trans = transpose(S1_double);
[num_row_C_0, num_col_C_0] = size(C_0);
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
end

