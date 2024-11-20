clear all
close all

load map_20010701_msradj_1519_2124.mat
load 

% Nominal Sensor Alignment
%%%%%%%%%%%%%%%%%%%%%%%%%%
T_b_g0   = -eye(3);             % from nominal gyro, IRUA1 coordinate to body coordinate
T_g0_b   = T_b_g0';
T_b_dss0 = [  9.999996032678291e-001  8.117507288903545e-004   3.500425151614639e-004
             -8.116613295096767e-004  9.999996454193570e-001  -2.560902541433511e-004
             -3.502595848992596e-004  2.557916040517020e-004   9.999999021731627e-001 ]';          % from nominal DSS1 coordinate to body coordinate
T_dss0_b = T_b_dss0';
T_b_st0  = [  0  1  0           % from nominal AST1 coordinate to body coordinate
              0  0  1
              1  0  0 ];
T_st0_b  = T_b_st0';

