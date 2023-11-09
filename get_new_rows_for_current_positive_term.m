function [C_new,V_new] = get_new_rows_for_current_positive_term(i,C_0_shifted,new_positive_terms_idx,new_rows_per_positive_row,V_0shifted,)
%UNTITLED Summary of this function goes here
C_new = C_0_shifted;
V_new = V_0shifted;
for k = 1:new_rows_per_positive_row
%     idx_new_row = new_positive_terms_idx(i)+k-1;
    new_positive_row = C_0_shifted(new_positive_terms_idx(i),:).*-negative_factors(k)./positive_factors(i)+C_new(negative_terms_idx(k),:)
    
    C_new(new_positive_terms_idx(i)+k-1,:) = new_positive_row;
    
    
    new_V_row = V_1(new_positive_terms_idx(i),:).*-negative_factors(k)./positive_factors(i)+V_1(negative_terms_idx(k),:);
    V_new(new_positive_terms_idx(i)+k-1,:) = new_V_row;
end
end

