function [positive_factors,negative_factors] = find_current_factors(Current_column, positive_idx,negative_idx)
%factor needed for combos in column j. Divide by positive term and multiply by (-)negative
positive_factors = Current_column(positive_idx);
negative_factors = Current_column(negative_idx);

end


