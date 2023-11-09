
%COSBI V03
syms Ka F BP CL CL_R fup ful r_f Kmax estimated_emax scal drug_MIC estimated_exponent estimated_ec50       %general params
syms Vgu Vki Vli Vlu Vlim Vsp Vve Var Vot Vlc                                                                 %volume
syms Kpki Kpgu Kpsp Kpli Kplu kdiff_im kout_im Kpsp Kpot kdiff_c kout_c k_consumption                         %constants
syms Qgu Qki Qha Qgu Qsp Qlu Qot Qh                                                                           %rates

% general_params1 =
% 
%    1.0e+03 *
% 
%   Columns 1 through 9
% 
%     0.0001    0.0006    0.0010    2.2947         0    0.0000    0.0000    0.0000    0.0071
% 
%   Columns 10 through 14
% 
%     0.0001    0.0000    0.0001    0.0020    0.3765

general_params = [Ka F BP CL CL_R fup ful r_f Kmax estimated_emax scal drug_MIC estimated_exponent estimated_ec50];
volume_params = [Vgu Vki Vli Vlu Vlim Vsp Vve Var Vot Vlc];
reaction_constants = [Kpki Kpgu Kpsp Kpli Kplu kdiff_im kout_im Kpsp Kpot kdiff_c kout_c k_consumption];
reaction_rates = [Qgu Qki Qha Qgu Qsp Qlu Qot Qh];
x = sym('x', [1 10]); 
%Kpsp
 S(1,1)   = sym(1);                                                             %Source to x(1),f1            source rate changes with oral dose(x(9))
 f(1) = sym(1);
 S(1,2)   = (1/Vgu)*(-1)*(x(1)/Kpgu) * BP;                                   %out x(1) to x(3), f2
 S(3,2)   = (1/Vli)*(x(1)/Kpgu) *BP;                                         %x(1) into x(3),f2
 f(2) = Qha*Qgu;
 S(2,3)   = (1/Vki)*(-1)*(x(2)/Kpki) * BP;                                   %out x(2) to x(6), f3 
 f(3) = Qki;
 S(2,4)   = (1/Vki)*(-1)*x(2) * fup;                                         %out x(2) to excrete, f4
 f(4) = CL*CL_R;
 S(6,3)   = (1/Vve) * (x(2)/Kpki) * BP ;                                     %x(2) into x(6),f3
 %f(3) = Qki;
 S(3,5)   = (1/Vli)*(-1)*(x(3)/Kpli) * BP;                                   %out x(3) to x(6), f5
 f(5) = Qh;
 S(6,5)   = (1/Vve) * ((x(3)/Kpli) * BP);                                    %x(3) into x(6), f5
 %f(5) = Qh;
 S(4,6)   = (1/Vlu)*(-1)*(ful*x(4)/Kplu) *BP;                                %out x(4) to x(7), f6;
 f(6) = Qlu;
 S(7,6)   = (1/Var) * (ful*(x(4)/Kplu) * BP);                                %x(4) into x(7), f6;
 %f(6) = Qlu;
 S(4,7)   = (1/Vlu) * (-1)*x(4);                                             %out x(4) to x(9),f7
 f(7) = kdiff_im;
 S(9,7)  = (Vlu/Vlim) * x(4);                                               %x(4) into x(9),f7
 %f(7) = kdiffim;
 S(5,8)   = (1/Vsp)*(-1)*(x(5)/Kpsp) * BP;                                   %out x(5) to x(3),f8
 f(8) = Qsp;
 S(3,8)   = (1/Vli) * (x(5)/Kpsp) * BP;                                  %x(5) into x(3), f8
 %f(8) = Qsp;
 S(6,9)   = (1/Vve) * (-1)* x(6);                                            %out x(6) to x(4),f9
 f(9) = Qlu;
 S(4,9)   = (1/Vlu) * x(6);                                              %x(6) into x(4),f9
 %f(9) = Qlu;
 S(7,10)  = (1/Var) * (-1)* x(7);                              %out x(7) to x(1),f10 maybe multiply by 1/5 for all out x(7) ...
 %f(10) = ;
 S(7,11)  = (1/Var) *  (-1)* x(7);                            %out x(7) to x(2),f11
 %f(11) = Qlu/Qki;
 S(7,12)  = (1/Var) * (-1)* x(7);                              %out x(7) to x(3),f12
 %f(12) = Qlu;
 S(7,13)  = (1/Var) * (-1)* x(7);                                  %out x(7) to x(5),f13
 %f(13) = Qlu;
 S(7,14)  = (1/Var) * (-1)* x(7);                                  %out x(7) to x(8),f14
 %f(14) = Qlu;
 S(1,10)  = (1/Vgu) * x(7);                                                  %x(7) into x(1),f10 
 f(10) = Qgu;
 S(2,11)  = (1/Vki) * x(7);                                                  %x(7) into x(2),f11 
 f(11) = Qki;
 S(3,12)  = (1/Vli) * x(7);                                                  %x(7) into x(3),f12 
 f(12) = Qha;
 S(5,13)  = (1/Vsp) * x(7);                                                  %x(7) into x(5),f13 
 f(13) = Qsp;
 S(8,14)  = (1/Vot) * x(7);                                              %x(7) into x(8),f14 
 f(14) = Qot;
 S(8,15)  = (1/Vot) * (-1)*(x(8)/Kpot) * BP;                             %out x(8) to x(6), f15
 f(15) = Qot;
 S(6,15)  = (1/Vve) * (x(8)/Kpot) * BP;                                  %x(8) into x(6), f15
 %f(15) = Qot;
 S(9,16) = (-1)* x(9);                                 %out x(9) to x(4), f16 previously S(9,16) = (Vlim/Vlu)*(-1)* x(9)
 f(16) = kout_im;
 S(4,16)  = (Vlim/Vlu)*x(9);                                        %x(9) into x(4), f16
 %f(16) = kout_im;
 %S(10,17) = (Vlc/Vlim) * kout_c * x(1)1;                                     %caseum into x(10), f17          source (since mean(x(1)1)= 3.5643e-14 we ignore
 f(17) = kout_c;
 %S(10,17)  = 0; 
 f = transpose(f);
 N_S = null(S);
 
 % read in parameter values
 
general_params1 =   1.0e+03 * [0.0001    0.0006    0.0010    2.2947         0    0.0000    0.0000    0.0000    0.0071    0.0001    0.0000    0.0001    0.0020    0.3765];
    
volume_params1 =    [0.0004    0.0001    0.0005    0.0002    0.0000    0.0001    0.0013    0.0006    0.0218         0];    
% volume_params(5) is 0 which causes divide by zero error.
volume_params1(5) = 1e-4;
reaction_constants1 = [0.0513    0.0455    0.0472    0.0503    0.0488         0         0    0.0472    0.0400 0         0         0]; 
reaction_rates1 =    [0.1519    0.1971    0.0536    0.1519    0.0179    1.0374    0.6168    0.2234];
 
 S1 = subs(S,x,ones(size(x)));
 S1_no_parameters = subs(S1,[general_params volume_params reaction_constants reaction_rates],ones(size([general_params volume_params reaction_constants reaction_rates])));
 S2 = subs(S1,[general_params volume_params reaction_constants reaction_rates],[general_params1 volume_params1 reaction_constants1 reaction_rates1]);