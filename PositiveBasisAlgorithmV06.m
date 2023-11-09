%% function implementation to perform positive basis algorithm
%Currently working on this version
% even though there is a version V04.
clear all,clc
syms Ka F BP CL CL_R fup ful r_f Kmax estimated_emax scal drug_MIC estimated_exponent estimated_ec50       %general params
syms Vgu Vki Vli Vlu Vlim Vsp Vve Var Vot Vlc                                                                 %volume
syms Kpki Kpgu Kpsp Kpli Kplu kdiff_im kout_im Kpsp Kpot kdiff_c kout_c k_consumption                         %constants
syms Qgu Qki Qha Qgu Qsp Qlu Qot Qh  

Rate_constants = [Ka F BP CL CL_R fup ful r_f Kmax estimated_emax scal drug_MIC estimated_exponent estimated_ec50];
V_constants = [Vgu Vki Vli Vlu Vlim Vsp Vve Var Vot Vlc];
K_constants = [Kpki Kpgu Kpsp Kpli Kplu kdiff_im kout_im Kpsp Kpot kdiff_c kout_c k_consumption];
f_constants = [Qgu Qki Qha Qgu Qsp Qlu Qot Qh];
x_variables = sym('x',[1 9]);
All_variables = [Rate_constants,V_constants,K_constants,f_constants,x_variables];
%% load the stoichiometric matrix

Test_Stoich_COSBI %COSBI system

%Test_Stoich_RCT   %RCT system
S1_double = double(S1_no_parameters);
num_x = size(S1_double,1);
num_f = size(S1_double,2);
V_0 = eye(num_f);
C_0 = V_0*transpose(S);
C_1 = C_0
size_original=size(C_1)
C_new = C_1;
V_1 = sym(eye(num_f));
V_new = V_1;

%% Algorithm
[Columns_to_change,src_snk_idx_for_C_0] = find_columns_to_change(C_0,S,S1_double);
for j=1:length(Columns_to_change)
    % find positive and negative terms
    if j == 2
        pause = 1;
    end
    Current_column = C_new(:,Columns_to_change(j));
    [positive_terms_idx,negative_terms_idx] = find_pos_terms(Current_column, ...
        All_variables);
    %find factors we will use to change rows
    [positive_factors,negative_factors] = find_current_factors(Current_column, ...
        positive_terms_idx,negative_terms_idx);
    %total changes in this column
    net_change_number_of_rows = length(positive_terms_idx)*length(negative_terms_idx) ...
        -(length(positive_terms_idx)+length(negative_terms_idx))
    new_rows_per_positive_row = length(negative_terms_idx);
    [C_0_shifted,V_1shifted,new_positive_terms_idx,new_negative_terms_idx] = ...
        shift_down_rows(C_1,V_1,positive_terms_idx,negative_terms_idx,new_rows_per_positive_row);
    shift_expected_change = net_change_number_of_rows;
%     shifted_size = size(C_0_shifted)
%     negative_terms_idx
%     new_negative_terms_idx
    C_new = C_0_shifted;
    V_new = V_1shifted;

    %the changes must be different
    
    for i=1:length(new_positive_terms_idx)
        new_positive_row = sym([]);
        new_V_row = sym([]);
        for k = 1:new_rows_per_positive_row
            new_positive_row(k,:) = C_new(new_positive_terms_idx(i),:).*-negative_factors(k)./positive_factors(i)+C_new(new_negative_terms_idx(k),:);
            new_V_row(k,:) = V_new(new_positive_terms_idx(i),:).*-negative_factors(k)./positive_factors(i)+V_new(new_negative_terms_idx(k),:);
        end
        C_new((new_positive_terms_idx(i):new_positive_terms_idx(i)+k-1),:) = new_positive_row;
        V_new((new_positive_terms_idx(i):new_positive_terms_idx(i)+k-1),:) = new_V_row;
    end
    
    C_new(new_negative_terms_idx,:) = [];
    V_new(new_negative_terms_idx,:) = [];

    C_1 = C_new
    net_change_num_rows = net_change_number_of_rows;
    new_size = size(C_1)
%     V_1
%     size(V_1)
    V_1 = V_new;
    new_V_size = size(V_1);
end
%% Now do alg for the sources and sink columns; same process starting at Right most col

for j=flip(1:length(src_snk_idx_for_C_0)) %reverse the order to go right to left along columns
    Current_column = C_new(:,src_snk_idx_for_C_0(j));
    [positive_terms_idx,negative_terms_idx] = find_pos_terms(Current_column, ...
        All_variables);
    %find factors we will use to change rows
    [positive_factors,negative_factors] = find_current_factors(Current_column, ...
        positive_terms_idx,negative_terms_idx);
    %total changes in this column
    net_change_number_of_rows = length(positive_terms_idx)*length(negative_terms_idx) ...
        -(length(positive_terms_idx)+length(negative_terms_idx))
    new_rows_per_positive_row = length(negative_terms_idx);
    [C_0_shifted,V_1shifted,new_positive_terms_idx,new_negative_terms_idx] = ...
        shift_down_rows(C_1,V_1,positive_terms_idx,negative_terms_idx,new_rows_per_positive_row);
    shift_expected_change = net_change_number_of_rows;
%     shifted_size = size(C_0_shifted)
%     negative_terms_idx
%     new_negative_terms_idx
    C_new = C_0_shifted;
    V_new = V_1shifted;
    
    for i=1:length(new_positive_terms_idx)
        new_positive_row = sym([]);
        new_V_row = sym([]);
        for k = 1:new_rows_per_positive_row
            new_positive_row(k,:) = C_new(new_positive_terms_idx(i),:).*-negative_factors(k)./positive_factors(i)+C_new(new_negative_terms_idx(k),:);
            new_V_row(k,:) = V_new(new_positive_terms_idx(i),:).*-negative_factors(k)./positive_factors(i)+V_new(new_negative_terms_idx(k),:);
        end
        C_new((new_positive_terms_idx(i):new_positive_terms_idx(i)+k-1),:) = new_positive_row;
        V_new((new_positive_terms_idx(i):new_positive_terms_idx(i)+k-1),:) = new_V_row;
    end
    
    C_new(new_negative_terms_idx,:) = [];
    V_new(new_negative_terms_idx,:) = [];

    C_1 = C_new
    new_size = size(C_1)
    V_1 = V_new;
    new_V_size = size(V_1)  ;  
end
    



