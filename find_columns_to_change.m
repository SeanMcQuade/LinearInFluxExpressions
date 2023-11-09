function [Columns_to_change,src_snk_idx_for_C_0] = find_columns_to_change(C_0,S,S1_double)
%find which entries of S are sources and sinks so that I can focus on the
%non source and sink columns of eye(num_f)*(S_transpose) = (C_0)
%   Detailed explanation goes here


[num_row_C_0, num_col_C_0] = size(C_0);

%find columns of S with sources and sinks using S1_no_parameters
src_snk_idx = [];
for i=1:num_row_C_0
    if nnz(S1_double(:,i))<2
        src_snk_idx = [src_snk_idx, i];
    end
end

%Rows corresponding to non_src_snk_idx should be sources and sinks
%find source and sink elements
for i = 1:length(src_snk_idx)
    col_of_src_snk = find(C_0(src_snk_idx(i),:));
    src_snk_element(i,:) = [src_snk_idx(i),col_of_src_snk];
end

src_snk_idx_for_C_0 = src_snk_element(:,2);
% search for columns of C_0 without src_snk_element. these are columns of
% C_0 not in src_snk_element(:,2)
Columns_to_change = setdiff(1:num_col_C_0,src_snk_element(:,2));

end

