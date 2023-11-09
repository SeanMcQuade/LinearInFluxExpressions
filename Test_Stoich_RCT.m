
%COSBI V03
x = sym('x',[1 6]);                                                                         %rates

S_RCT_transpose =[  1,   0,   0,   0,   0,   0
  0,   1,   0,   0,   0,   0
  0,   0,   1,   0,   0,   0
-x(1),   0,   0,  x(1),   0,   0
  0, -x(2),   0,  x(2),   0,   0
  0,   0, -x(3),  x(3),   0,   0
  0,   0,   0, -x(4),  x(4),   0
  0,   0,   0, -x(4),   0,  x(4)
  0,   0,   0,   0, -x(5),  x(5)
  0,   0,   0,   0,   0, -x(6)];

S = transpose(S_RCT_transpose);


 S1 = subs(S,x,ones(size(x)));
 S1_no_parameters = S1;
 %S2 = subs(S1,[general_params volume_params reaction_constants reaction_rates],[general_params1 volume_params1 reaction_constants1 reaction_rates1]);