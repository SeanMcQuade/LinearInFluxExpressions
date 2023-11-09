function [C_0shifted,V_1shifted,new_positive_terms_idx,new_negative_terms_idx] = shift_down_rows(C_0,V_1,positive_terms_idx,negative_terms_idx,new_rows_per_positive_row)
%shift down rows for extra space needed for the algorithm, must shift V
%rows as well
[num_rows,num_cols] = size(C_0);
[num_V_rows,num_V_cols] = size(V_1);
i=1;
C_0shifted = C_0;
V_1shifted = V_1;
while i < length(positive_terms_idx)+1
    shift_amount = new_rows_per_positive_row-1;
    start_idx = positive_terms_idx(i)+1;
    end_idx = length(C_0);
    C_0shifted((start_idx+shift_amount:end_idx+shift_amount),:) = C_0((start_idx:end_idx),:);
    V_1shifted((start_idx+shift_amount:end_idx+shift_amount),:) = V_1shifted((start_idx:end_idx),:);
    %add new rows of zeros
    C_0shifted((start_idx:start_idx+shift_amount-1),:) = zeros(new_rows_per_positive_row-1,num_cols);
    V_1shifted((start_idx:start_idx+shift_amount-1),:) = zeros(new_rows_per_positive_row-1,num_V_cols);
    C_0 = C_0shifted;
    V_1 = V_1shifted;
    %number_new_rows = length(C_0shifted);
    positive_terms_idx(i+1:end) = positive_terms_idx(i+1:end)+shift_amount;
    %which negative indices to shift?
    neg_indices_to_move = positive_terms_idx(i)<negative_terms_idx;
    negative_terms_idx(neg_indices_to_move) = negative_terms_idx(neg_indices_to_move)+shift_amount;
    i=i+1;
    %keep track of shifting of negative_idices
end
new_positive_terms_idx = positive_terms_idx;
new_negative_terms_idx = negative_terms_idx;
end
