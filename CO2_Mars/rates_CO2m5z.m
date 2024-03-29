function [X,R,dTe,dS_ape,dS_el,dS_epro,dS_rea,dS_wall,K0,Rv,Kv0,Ralldlt]=rates_CO2m5z(Te,n,area_volume,Tg,me,e,S_abs,dt,Vibration,K,nsp,number_of_reactions,gas_heat,mass,aux,auxv,Kv,number_of_reaction,number_of_vreaction,nv,sim)

% R1:  e + CO2 => CO2^+ + 2e          k= 查表                                                ER=13.8eV
% R2:  e + CO2 => CO + O + e          k= 查表                   (激发解离)                   ER=11.46eV
% R3:  e + CO2 => CO + O^-            k= 查表                                            （附着反应无ER）
% R4:  e + O3  => O + O2 + e          k= 9.0e-16                (激发解离)           (应该存在ER，但未找到数据)
% R5:  e + O2  => 2O + e              k= 查表            (生成的O一般而言属于激发态)            ER=6eV
% R6:  e + O2  => O + O^-             k= 查表                                             (附着反应无ER)
% R7:  O^- + CO => CO2 + e            k= 5.5e-16
% R8:  O^- + O3 => 2O2 + e            k= 3.0e-16
% R9:  e + CO2^+ => CO + O            k= 2.0e-11*Te^(-0.5)*(Tg*1.16e4)^(-1)
% R10: O2^- + CO2^+ => CO + O2 + O    k= 6.0e-13
% R11: O + O2 + CO2 => O3 + CO2       k= 1.7e-42*(Tg*1.16e4)^(-1.2)                       (O+O2+M,M为CO2)
% R12: O^- + O2 => O + e + O2         k= 6.9e-16                                           (O- + M ,M为O2)
% R13: e + CO2 => e + CO2             k= 查表                                            （弹性碰撞反应ER=0）
% R14: O + O2^- => O2 + O^-           k= 3.31e-16
% R15: O^- + CO2 => O + CO2 + e       k= 4.0e-18                                          （O- + M ,M为CO2）
% R16: O^- + O => O2 + e              k= 2.3e-16
% R17: e + CO => C + O^-              k= 查表                                                 ER=9.2eV
% R18: e + CO => e + C + O            k= 查表                                                 ER=11.1eV
% R19: e + CO2^+ => C + O2            k= 3.94e-13*Te^(-0.4)
% R20: CO2 + C => 2CO                 k= 1.0e-21
% R21: C + O2 => O + CO               k= 3.0e-17
% R22: e + O2 => 2e + O2^+            k= 查表                                               ER=12.06eV
% R23: e + O2^+ => 2O                 k= 6.0e-13*(1/Te)^0.5*(1/(Tg*1.16e4))^0.5
% R24: O^- + O2^+ => O2 + O           k= 2.6e-14*(300/(Tg*1.16e4))^0.44
% R25: O^- + O2^+ => 3O               k= 4.2e-13*(300/(Tg*1.16e4))^0.44
% R26: O2^+ + O2^- => 2O2             k= 2.01e-13*(300/(Tg*1.16e4))^0.5
% R27: O2^+ + O2^- => O2 + 2O         k= 4.2e-13
% R28: O^- + CO2 + CO2 => CO3^- + CO2 k= 9.0e-41                                       （O- + CO2 + M,M为CO2）
% R29: CO3^- + CO => 2CO2 + e         k= 5.0e-19
% R30: CO3^- + CO2^+ => 2CO2 + O      k= 5.0e-13                                       （生成的CO2为vb振动态)
% R31: O + CO3^- => CO2 + O2^-        k= 8.0e-17
% R32: CO2 + e => O2^+ + C + 2e       k= K(1)*1/13                               ER未知,无法对应给出CO2e的电离阈值
% R33: O3 + e => O2^+ + O + 2e        k= 3.2e-17*(Te*1.16e4)^(0.5)*(1+0.15*Te)*exp(-12.93/Te)
% R34: O3 + e => O2 + O^-             k= 查表                                             （附着反应无ER）
% R35: CO2^+ + O => O2^+ + CO         k= 0.63*2.6e-16
% R36: CO2^+ + O2 => O2^+ + CO2       k= 6.4e-17
% R37: O3 + e => O2^- + O             k= 查表                                             （附着反应无ER）
% R38: CO3^- + O2^+ => CO2 + O2 + O   k= 3.0e-13
% R39: O2^- + O => O3 + e             k= 3.3e-16
% R40: CO2 + O^- + CO => CO3^- + CO   k= 1.5e-40                                       (O- + CO2 + M,M为CO)
% R41: CO2 + O^- + O2 => CO3^- + O2   k= 3.1e-40                                       （O- + CO2 + M,M为O2）
% R42: O + O + CO2 => O2 + CO2        k= 3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4))（O+O+M，其中M为CO2）
% R43: O^- + CO2^+ => O + CO2         k= 1.0e-13
% R44: e + CO2 => 2e + O + CO^+       k= 查表                                                ER=19.5eV
% R45: e + CO => 2e + CO^+            k= 查表                                                ER=14.01eV
% R46: CO2 + CO^+ => CO2^+ + CO       k= 1.0e-15
% R47: O2^+ + C => CO^+ + O           k= 5.2e-17
% R48: e + O3 + O2 => O3^- + O2       k= 4.6e-40                                          (e+O3+M,M为O2)
% R49: O3^- + CO2 => CO3^- + O2       k= 5.5e-16
% R50: O^- + O3 => O3^- + O           k= 8.0e-16
% R51: O2^- + O3 => O3^- + O2         k= 4.0e-16
% R52: O3^- + O => O2^- + O2          k= 2.5e-16
% R53: e + CO2 => 2e + CO + O^+       k= 查表                                                ER=19.1eV
% R54: e + O => 2e + O^+              k= 查表                                                ER=13.6eV
% R55: O^+ + CO2 => O2^+ + CO         k= 9.4e-16
% R56: O^+ + CO2 => O + CO2^+         k= 4.5e-16
% R57: CO2^+ + O => CO2 + O^+         k= 9.62e-17
% R58: e + CO => 2e + C^+ + O         k= 查表                                                ER=22eV
% R59: e + C => 2e + C^+              k= 查表                                               ER=11.2eV
% R60: C^+ + CO2 => CO^+ + CO         k= 1.1e-15
% R61: O2^+ + C => C^+ + O2           k= 5.2e-17
% R62: O2^- + CO2 + O2 => CO4^- + O2  k= 4.7e-41                                      (O2- + CO2 + M,M为O2)
% R63: CO4^- + CO2^+ => 2CO2 + O2     k= 5.0e-13
% R64: O2^+ + CO4^- => CO2 + 2O2      k= 3.0e-13
% R65: CO4^- + O => CO3^- + O2        k= 1.1e-16
% R66: CO4^- + O => CO2 + O2 + O^-    k= 1.4e-17
% R67: CO4^- + O => CO2 + O3^-        k= 1.4e-16
% R68: O^- + C => CO + e              k= 5.0e-16
% R69: e + CO4^+ => CO2 + O2          k= 1.61e-13*Te^(-0.5)
% R70: O2^+ + CO2 + CO2 => CO4^+ + CO2   k= 2.3e-41
% R71: e + O2^+ + CO2 => O2 + CO2     k=1.0e-38                                       (e + O2+ + M,M为CO2)
% R72: O2^- + CO2 + CO2 => CO4^- + CO2   k= 1.0e-41                         (O2- + CO2 + M.M为CO2)  ------1.339606
% R73: e + O2 + CO2 => O2^- + CO2        k= 2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4))
% R74: e + C2O2^+ => 2CO              k= 4.0e-13*Te^(-0.34) 
% R75: e + C2O3^+ => CO2 + CO         k= 5.4e-14*Te^(-0.7)
% R76: e + C2O4^+ => 2CO2             k= 2.0e-11*(Te^(-0.5))*(Tg*1.16e4)^(-1)
% R77: CO2^+ + CO2 + CO => C2O4^+ + CO  k= 3.0e-40
% R78: C2O3^+ + CO => CO2 + C2O2^+    k= 1.1e-15
% R79: C2O4^+ + CO => CO2 + C2O3^+    k= 9.0e-16
% R80: C2O3^+ + CO + CO2 => CO2 + C2O2^+ + CO2  k= 2.6e-38
% R81: C2O4^+ + CO + CO2 => CO2 + C2O3^+ + CO2  k= 4.2e-38
% R82: C2O2^+ + O2 => 2CO + O2^+      k= 5.0e-18
% R83: C2O2^+ + CO2 => CO^+ + CO + CO2    k= 1.0e-18
% R84: C2O2^+ + CO3^- => CO2 + 2CO + O       k= 5.0e-13
% R85: C2O2^+ + CO4^- => CO2 + 2CO + O2      k= 5.0e-13
% R86: C2O2^+ + O2^- => 2CO + O2      k= 6.0e-13
% R87: C2O3^+ + CO3^- => 2CO2 + CO + O       k= 5.0e-13
% R88: C2O3^+ + CO4^- => 2CO2 + CO + O2      k= 5.0e-13
% R89: C2O3^+ + O2^- => CO2 + CO + O2        k= 6.0e-13
% R90: C2O4^+ + CO2 => CO2^+ + CO2 + CO2  k= 1.0e-20
% R91: C2O4^+ + CO3^- => 3CO2 + O     k= 5.0e-13
% R92: C2O4^+ + CO4^- => 3CO2 + O2    k= 5.0e-13
% R93: C2O4^+ + O2^- => 2CO2 + O2     k= 6.0e-13
% R94: e + N2 => e + N2               k= 查表
% R95: e + N2 => 2e + N2^+            k= 查表                                          ER=15.58eV
% R96: e + N2 => e + 2N               k= 查表                                          ER=9.753eV
% R97: e + NO => N + O^-              k= 查表    
% R98: e + NO => 2e + NO^+            k= 查表                                          ER=9.26eV
% R99: e + N2O => N2 + O^-            k= 查表   
% R100:e + N2O => 2e + N2O^+          k= 查表                                          ER=12.89eV
% R101:e + NO2 => NO2^+ + 2e          k= 2.6e-15*Te^1/2*exp(-10/Te)                     ER=10eV
% R102:e + NO2 => NO^+ + O + 2e       k= 8.1e-15*Te^1/2*exp(-12.9/Te)                  ER=12.9eV
% R103:e + NO => N + O + e            k= 7.4e-15*exp(-6.5/Te)                           ER=6.5eV
% R104:e + NO2 => NO + O + e          k= 5.6e-15*exp(-3.11/Te)                         ER=3.11eV
% R105:e + N2O => N2 + O + e          k= 1.4e-15*exp(-1.67/Te)                         ER=1.67eV
% R106:e + NO2 + N2 => NO2^- + N2     k= 1.6e-43
% R107:e + NO2 + CO2 => NO2^- + CO2   k= 1.6e-43
% R108:O^- + N => NO + e              k= 2.2e-16
% R109:O^- + N2 => N2O + e            k= 1.0e-18
% R110:O^- + NO => NO2 + e            k= 2.5e-16
% R111:O2^- + N => NO2 + e            k= 4.0e-16
% R112:NO2^- + O => NO3 + e           k= 1.0e-18
% R113:NO2^- + N => N2 + O2 + e       k= 1.0e-18
% R114:NO3^- + O => NO2 + O2 + e      k= 1.0e-18
% R115:NO3^- + N => N2 + O3 + e       k= 1.0e-18
% R116:NO^- + N2 => NO + N2 + e       k= 5.0e-16
% R117:e + NO^+ => N + O              k= 7.0e-14*Te^(-0.5)
% R118:e + N2^+ => 2N                 k= 4.6e-14*Te^(-0.39)
% R119:e + NO2^+ => N + O2            k= 7.0e-14*Te^(-0.5)
% R120:e + N2O^+ => N2 + O            k= 7.0e-14*Te^(-0.5)
% R121:CO2^+ + NO => NO^+ + CO2       k= 1.2e-16
% R122:N2^+ + O2 => O2^+ + N2         k= 7.0e-17
% R123:N2^+ + NO => NO^+ + N2         k= 4.9e-16
% R124:N2^+ + CO => CO^+ + N2         k= 7.0e-17
% R125:N2^+ + CO2 => CO2^+ + N2       k= 9.0e-16
% R126:O2^+ + N => NO^+ + O           k= 1.9e-17
% R127:O2^+ + NO => NO^+ + O2         k= 7.0e-16
% R128:O2^+ + NO2 => NO2^+ + O2       k= 6.6e-16
% R129:NO^+ + O3 => NO2^+ + O2        k= 1.0e-20
% R130:NO2^+ + NO => NO^+ + NO2       k= 2.9e-16
% R131:N2O^+ + NO => NO^+ + N2O       k= 2.9e-16
% R132:O^- + NO2 => O2^- + NO         k= 1.0e-15
% R133:O^- + NO2 => NO2^- + O         k= 1.0e-15
% R134:O^- + N2O => NO^- + NO         k= 2.0e-16
% R135:O2^- + NO2 => NO2^- + O2       k= 2.0e-15
% R136:O2^- + N2O => O3^- + N2        k= 1.0e-17
% R137:O2^- + NO3 => NO3^- + O2       k= 5.0e-16
% R138:O3^- + NO => NO3^- + O         k= 1.0e-17
% R139:O3^- + NO => NO2^- + O2        k= 2.6e-18
% R140:O3^- + NO2 => NO2^- + O3       k= 7.0e-16
% R141:O3^- + NO2 => NO3^- + O2       k= 2.8e-16
% R142:O3^- + NO3 => NO3^- + O3       k= 5.0e-16
% R143:CO3^- + NO => NO2^- + CO2      k= 9.0e-18
% R144:CO3^- + NO2 => NO3^- + CO2     k= 1.0e-16
% R145:CO4^- + NO => NO3^- + CO2      k= 5.0e-17
% R146:NO^- + O2 => O2^- + NO         k= 9.0e-16
% R147:NO2^- + O3 => NO3^- + O2       k= 1.8e-17
% R148:NO2^- + NO2 => NO3^- + NO      k= 4.0e-18
% R149:NO2^- + NO3 => NO3^- + NO2     k= 3.0e-16
% R150:NO3^- + NO => NO2^- + NO2      k= 2.0e-17
% R151:O^- + NO + N2 => NO2^- + N2    k= 3.0e-40*(Tg*11600)^(-1)
% R152:CO2^+ + NO2^- => CO + NO2 + O  k= 6.0e-13
% R153:CO2^+ + NO3^- => CO + NO3 + O  k= 5.0e-13
% R154:NO^+ + CO3^- => CO2 + NO + O   k= 6.0e-13
% R155:NO^+ + CO4^- => CO2 + NO + O2  k= 6.0e-13
% R156:NO^+ + NO2^- => NO2 + N + O    k= 5.1e-13
% R157:NO^+ + NO3^- => NO3 + N + O    k= 8.1e-13
% R158:NO^+ + O2^- => O2 + N + O      k= 5.8e-13
% R159:NO^+ + O^- => NO + O           k= 4.9e-13
% R160:O2^+ + NO2^- => NO2 + 2O       k= 4.1e-13
% R161:O2^+ + NO3^- => NO3 + 2O       k= 1.3e-13
% R162:N2^+ + CO3^- => CO2 + N2 + O   k= 3.0e-13
% R163:N2^+ + CO4^- => CO2 + N2 + O2  k= 3.0e-13
% R164:N2^+ + NO2^- => NO2 + 2N       k= 4.1e-13
% R165:N2^+ + NO3^- => NO3 + 2N       k= 1.3e-13
% R166:N2^+ + O2^- => N2 + 2O         k= 4.2e-13
% R167:N2^+ + O^- => N2 + O           k= 1.0e-13
% R168:C2O4^+ + NO2^- => 2CO2 + NO2   k= 6.0e-13
% R169:C2O4^+ + NO3^- => 2CO2 + NO2 + O  k= 5.0e-13
% R170:C2O3^+ + NO2^- => CO2 + CO + NO2  k= 6.0e-13
% R171:C2O3^+ + NO3^- => CO2 + CO + NO2 + O  k=5.0e-13
% R172:C2O2^+ + NO2^- => 2CO + NO2    k= 6.0e-13
% R173:C2O2^+ + NO3^- => 2CO + NO2 + O   k=5.0e-13
% R174:N + CO2 => NO + CO             k= 1.0e-25
% R175:CO + NO2 => CO2 + NO           k= 1.48e-16*exp(-16967/(Tg*1.16e4))
% R176:NO + O + CO2 => CO2 + NO2      k= 1.13e-37*(Tg*1.16e4)^(-2.16)*exp(-529/(Tg*1.16e4))
% R177:NO^- + CO => NO + CO + e       k= 5.0e-19
% R178:NO^- + CO2 => NO + CO2 + e     k= 8.3e-18
% 振动激发态部分
% Rv1: e + CO2 => e + CO2e1           k= 查表                                                EE=7eV  
% Rv2: e + CO2 => e + CO2e2           k= 查表                                              EE=10.5eV
% Rv3: e + CO2e1 => CO2^+ + 2e        k= K(1)*0.944                                    EI=6.8eV(未必准确)
% Rv4: e + CO2e2 => CO2^+ + 2e        k= K(1)*0.099                                    EI=3.3eV(未必准确)
% Rv5: e + CO2e1 => 2e + O + CO^+     k= K(44)*3.189                                   EI=12.5eV(未必准确)
% Rv6: e + CO2e2 => 2e + O + CO^+     k= K(44)*0.735                                    EI=9eV(未必准确)
% Rv7: CO2e1 + CO^+ => CO2^+ + CO     k= K(46)*0.916
% Rv8: CO2e2 + CO^+ => CO2^+ + CO     k= K(46)*0.093
% Rv9: e + CO2e1 => 2e + CO + O^+     k= K(53)*2.988                                   EI=12.1eV(未必准确)
% Rv10: e + CO2e2 => 2e + CO + O^+     k= K(53)*0.671                                  EI=8.6eV((未必准确))
% Rv11: O^+ + CO2e1 => O + CO2^+       k= K(56)*0.916
% Rv12: O^+ + CO2e2 => O + CO2^+       k= K(56)*0.093
% Rv13: C^+ + CO2e1 => CO^+ + CO       k= K(60)*3.189
% Rv14: C^+ + CO2e2 => CO^+ + CO       k= K(60)*0.735
% Rv15: e + CO2 => e + CO2va           k= 查表                                              EV=0.083
% Rv16: e + CO2 => e + CO2vb           k= 查表                                              EV=0.167
% Rv17: e + CO2 => e + CO2vc           k= 查表                                              EV=0.252
% Rv18: e + CO2 => e + CO2vd           k= 查表                                              EV=0.339
% Rv19: CO2va + CO2 => 2CO2            k= 7.14e-14*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv20: CO2va + CO => CO2 + CO         k= 5.00e-14*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv21: CO2va + O2 => CO2 + O2         k= 5.00e-14*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv22: CO2vb + CO2 => 2CO2            k= 1.071e-15*exp(-137*(Tg*1.16e4)^(-1/3))
% Rv23: CO2vb + CO => CO2 + CO         k= 3.213e-15*exp(-137*(Tg*1.16e4)^(-1/3))
% Rv24: CO2vb + O2 => CO2 + O2         k= 3.213e-15*exp(-137*(Tg*1.16e4)^(-1/3))
% Rv25: CO2vb + CO2 => CO2va + CO2     k= 1.438e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv26: CO2vb + CO => CO2va + CO       k= 1.007e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv27: CO2vb + O2 => CO2va + O2       k= 1.007e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv28: CO2vc + CO2 => CO2va + CO2     k= 1.071e-15*exp(-137*(Tg*1.16e4)^(-1/3))
% Rv29: CO2vc + CO => CO2va + CO       k= 3.213e-15*exp(-137*(Tg*1.16e4)^(-1/3))
% Rv30: CO2vc + O2 => CO2va + O2       k= 3.213e-15*exp(-137*(Tg*1.16e4)^(-1/3))
% Rv31: e + CO2 => e + CO2v1           k= 查表                                            EV=0.291eV
% Rv32: e + CO2 => e + CO2v2           k= 查表                                             EV=0.58eV
% Rv33: e + CO2 => e + CO2v3           k= 查表                                             EV=0.86eV
% Rv34: e + CO2 => e + CO2v4           k= 查表                                             EV=1.14eV
% Rv35: e + CO2v1 => CO2^+ + 2e        k= K(1)                                             ER=13.8eV
% Rv36: e + CO2v2 => CO2^+ + 2e        k= K(1)                                             ER=13.8eV
% Rv37: e + CO2v3 => CO2^+ + 2e        k= K(1)                                             ER=13.8eV
% Rv38: e + CO2v4 => CO2^+ + 2e        k= K(1)                                             ER=13.8eV
% Rv39: e + CO2v1 => CO + O + e        k= 查表                                            ER=11.169eV
% Rv40: e + CO2v2 => CO + O + e        k= 查表                                            ER=10.88eV
% Rv41: e + CO2v3 => CO + O + e        k= 查表                                             ER=10.6eV
% Rv42: e + CO2v4 => CO + O + e        k= 查表                                            ER=10.32eV
% Rv43: e + CO2v1 => CO + O-           k= 查表                                        附着过程无ER，下同
% Rv44: e + CO2v2 => CO + O-           k= 查表
% Rv45: e + CO2v3 => CO + O-           k= 查表                    
% Rv46: e + CO2v4 => CO + O-           k= 查表
% Rv47: e + CO2v1 => O2^+ + C + 2e     k= K(1)*1/13                                   ER未知，k值未修正
% Rv48: e + CO2v2 => O2^+ + C + 2e     k= K(1)*1/13                                   ER未知，k值未修正
% Rv49: e + CO2v3 => O2^+ + C + 2e     k= K(1)*1/13                                   ER未知，k值未修正
% Rv50: e + CO2v4 => O2^+ + C + 2e     k= K(1)*1/13                                   ER未知，k值未修正
% Rv51: e + CO2v1 => 2e + O + CO^+     k= 查表                                            ER=19.21eV
% Rv52: e + CO2v2 => 2e + O + CO^+     k= 查表                                            ER=18.92eV
% Rv53: e + CO2v3 => 2e + O + CO^+     k= 查表                                            ER=18.64eV
% Rv54: e + CO2v4 => 2e + O + CO^+     k= 查表                                            ER=18.36eV
% Rv55: e + CO2v1 => 2e + CO + O^+     k= 查表                                            ER=18.81eV
% Rv56: e + CO2v2 => 2e + CO + O^+     k= 查表                                            ER=18.52eV
% Rv57: e + CO2v3 => 2e + CO + O^+     k= 查表                                            ER=18.24eV
% Rv58: e + CO2v4 => 2e + CO + O^+     k= 查表                                            ER=17.96eV
% Rv59: CO2v1 + CO2 => CO2va + CO2     k= 4.25e-7*exp(-407*(Tg*11600)^(-1/3)+824*(Tg*11600)^(-2/3))
% Rv60: CO2v1 + CO => CO2va + CO       k= Kv(59)*0.3
% Rv61: CO2v1 + O2 => CO2va + O2       k= Kv(59)*0.4
% Rv62: CO2v1 + CO2 => CO2vb + CO2     k= 8.57e-7*exp(-404*(Tg*11600)^(-1/3)+1096*(Tg*11600)^(-2/3))
% Rv63: CO2v1 + CO => CO2vb + CO       k= Kv(62)*0.3
% Rv64: CO2v1 + O2 => CO2vb + O2       k= Kv(62)*0.4
% Rv65: CO2v1 + CO2 => CO2vc + CO2     k= 1.43e-11*exp(-252*(Tg*11600)^(-1/3)+685*(Tg*11600)^(-2/3))
% Rv66: CO2v1 + CO => CO2vc + CO       k= Kv(65)*0.3;
% Rv67: CO2v1 + O2 => CO2vc + O2       k= Kv(65)*0.4;
% Rv68: CO2v2 + CO2 => CO2v1 + CO2     k= 见下文
% Rv69: CO2v2 + CO => CO2v1 + CO       k= 见下文
% Rv70: CO2v2 + O2 => CO2v1 + O2       k= 见下文
% Rv71: CO2v3 + CO2 => CO2v2 + CO2     k= 见下文
% Rv72: CO2v3 + CO => CO2v2 + CO       k= 见下文
% Rv73: CO2v3 + O2 => CO2v2 + O2       k= 见下文
% Rv74: CO2v4 + CO2 => CO2v3 + CO2     k= 见下文
% Rv75: CO2v4 + CO => CO2v3 + CO       k= 见下文
% Rv76: CO2v4 + O2 => CO2v3 + O2       k= 见下文
% Rv77: CO2v1 + CO2 => CO2va + CO2vb   k= 1.06e-11*exp(-242*(Tg*11600)^(-1/3)+633*(Tg*11600)^(-2/3))
% Rv78: CO2v2 + CO2 => CO2v1 + CO2va   k= 2.13e-11*exp(-242*(Tg*11600)^(-1/3)+633*(Tg*11600)^(-2/3))
% Rv79: CO2v2 + CO2 => CO2v1 + CO2vb   k= 2.13e-11*exp(-242*(Tg*11600)^(-1/3)+633*(Tg*11600)^(-2/3))
% Rv80: CO2v3 + CO2 => CO2v2 + CO2va   k= 2.13e-11*exp(-242*(Tg*11600)^(-1/3)+633*(Tg*11600)^(-2/3))
% Rv81: CO2v3 + CO2 => CO2v2 + CO2vb   k= 2.13e-11*exp(-242*(Tg*11600)^(-1/3)+633*(Tg*11600)^(-2/3))
% Rv82: CO2v4 + CO2 => CO2v3 + CO2va   k= 2.13e-11*exp(-242*(Tg*11600)^(-1/3)+633*(Tg*11600)^(-2/3))
% Rv83: CO2v4 + CO2 => CO2v3 + CO2vb   k= 2.13e-11*exp(-242*(Tg*11600)^(-1/3)+633*(Tg*11600)^(-2/3))
% Rv84: CO2vb + CO2 => CO2va + CO2va   k= 2.157e-15*exp(-88*(Tg*11600)^(-1/3)+233*(Tg*11600)^(-2/3))
% Rv85: CO2vc + CO2 => CO2va + CO2vb   k= 5.305e-15*exp(-88.2*(Tg*11600)^(-1/3)+233*(Tg*11600)^(-2/3))
% Rv86: CO2vc + CO2va => CO2vb + CO2vb k= 5.305e-15*exp(-88*(Tg*11600)^(-1/3)+233*(Tg*11600)^(-2/3))
% Rv87: CO2va + CO2 => CO + O + CO2    k= 3.91e-16*exp(-48688/(Tg*11600))
% Rv88: CO2vb + CO2 => CO + O + CO2    k= 3.91e-16*exp(-47852/(Tg*11600))
% Rv89: CO2vc + CO2 => CO + O + CO2    k= 3.91e-16*exp(-47110/(Tg*11600))
% Rv90: CO2vd + CO2 => CO + O + CO2    k= 3.91e-16*exp(-46368/(Tg*11600))
% Rv91: CO2v1 + CO2 => CO + O + CO2    k= 3.91e-16*exp(-46739/(Tg*11600))
% Rv92: CO2v2 + CO2 => CO + O + CO2    k= 3.91e-16*exp(-44048/(Tg*11600))
% Rv93: CO2v3 + CO2 => CO + O + CO2    k= 3.91e-16*exp(-41449/(Tg*11600))
% Rv94: CO2v4 + CO2 => CO + O + CO2    k= 3.91e-16*exp(-38851/(Tg*11600))
% Rv95: CO2va + O => CO + O2           k= 2.8e-17*exp(-26036/(Tg*11600))
% Rv96: CO2vb + O => CO + O2           k= 2.8e-17*exp(-25514/(Tg*11600))
% Rv97: CO2vc + O => CO + O2           k= 2.8e-17*exp(-25050/(Tg*11600))
% Rv98: CO2vd + O => CO + O2           k= 2.8e-17*exp(-24586/(Tg*11600))
% Rv99: CO2v1 + O => CO + O2           k= 2.8e-17*exp(-24818/(Tg*11600))
% Rv100:CO2v2 + O => CO + O2           k= 2.8e-17*exp(-23136/(Tg*11600))
% Rv101:CO2v3 + O => CO + O2           k= 2.8e-17*exp(-21512/(Tg*11600))
% Rv102:CO2v4 + O => CO + O2           k= 2.8e-17*exp(-19888/(Tg*11600))
% Rv103:e + CO2va => CO + O + e        k= 查表                                            ER=11.37eV
% Rv104:e + CO2vb => CO + O + e        k= 查表                                            ER=11.293eV
% Rv105:e + CO2vc => CO + O + e        k= 查表                                            ER=11.208eV
% Rv106:e + CO2vd => CO + O + e        k= 查表                                            ER=11.121eV
% Rv107:e + CO2va => CO + O-           k= 查表
% Rv108:e + CO2vb => CO + O-           k= 查表
% Rv109:e + CO2vc => CO + O-           k= 查表
% Rv110:e + CO2vd => CO + O-           k= 查表
% Rv111:e + CO2va => CO2^+ + 2e        k= K(1)                                             ER=13.8eV
% Rv112:e + CO2vb => CO2^+ + 2e        k= K(1)                                             ER=13.8eV
% Rv113:e + CO2vc => CO2^+ + 2e        k= K(1)                                             ER=13.8eV
% Rv114:e + CO2vd => CO2^+ + 2e        k= K(1)                                             ER=13.8eV
% Rv115:e + CO2va => O2^+ + C + 2e     k= K(1)*1/13                                    ER未知，k值未修正
% Rv116:e + CO2vb => O2^+ + C + 2e     k= K(1)*1/13                                    ER未知，k值未修正
% Rv117:e + CO2vc => O2^+ + C + 2e     k= K(1)*1/13                                    ER未知，k值未修正
% Rv118:e + CO2vd => O2^+ + C + 2e     k= K(1)*1/13                                    ER未知，k值未修正
% Rv119:e + CO2va => 2e + O + CO^+     k= 查表                                            ER=19.417eV                        
% Rv120:e + CO2vb => 2e + O + CO^+     k= 查表                                            ER=19.333eV                   
% Rv121:e + CO2vc => 2e + O + CO^+     k= 查表                                            ER=19.248eV                                 
% Rv122:e + CO2vd => 2e + O + CO^+     k= 查表                                            ER=19.161eV
% Rv123:e + CO2va => 2e + CO + O^+     k= 查表                                            ER=19.017eV
% Rv124:e + CO2vb => 2e + CO + O^+     k= 查表                                            ER=18.933eV                                   
% Rv125:e + CO2vc => 2e + CO + O^+     k= 查表                                            ER=18.848eV                                 
% Rv126:e + CO2vd => 2e + CO + O^+     k= 查表                                            ER=18.761eV                                       
% Rv127:CO2e1 + C => 2CO               k= 1.0e-21
% Rv128:CO2e2 + C => 2CO               k= 1.0e-21
% Rv129:CO2va + C => 2CO               k= 1.0e-21
% Rv130:CO2vb + C => 2CO               k= 1.0e-21
% Rv131:CO2vc + C => 2CO               k= 1.0e-21
% Rv132:CO2vd + C => 2CO               k= 1.0e-21
% Rv133:CO2v1 + C => 2CO               k= 1.0e-21
% Rv134:CO2v2 + C => 2CO               k= 1.0e-21
% Rv135:CO2v3 + C => 2CO               k= 1.0e-21
% Rv136:CO2v4 + C => 2CO               k= 1.0e-21
% Rv137:O^- + CO2e1 + CO2 => CO3^- + CO2 k= 9.0e-41
% Rv138:O^- + CO2e2 + CO2 => CO3^- + CO2 k= 9.0e-41
% Rv139:O^- + CO2va + CO2 => CO3^- + CO2 k= 9.0e-41
% Rv140:O^- + CO2vb + CO2 => CO3^- + CO2 k= 9.0e-41
% Rv141:O^- + CO2vc + CO2 => CO3^- + CO2 k= 9.0e-41
% Rv142:O^- + CO2vd + CO2 => CO3^- + CO2 k= 9.0e-41
% Rv143:O^- + CO2v1 + CO2 => CO3^- + CO2 k= 9.0e-41
% Rv144:O^- + CO2v2 + CO2 => CO3^- + CO2 k= 9.0e-41
% Rv145:O^- + CO2v3 + CO2 => CO3^- + CO2 k= 9.0e-41
% Rv146:O^- + CO2v4 + CO2 => CO3^- + CO2 k= 9.0e-41
% Rv147:CO2e1 + O^- + CO => CO3^- + CO   k= 1.5e-40
% Rv148:CO2e2 + O^- + CO => CO3^- + CO   k= 1.5e-40
% Rv149:CO2va + O^- + CO => CO3^- + CO   k= 1.5e-40
% Rv150:CO2vb + O^- + CO => CO3^- + CO   k= 1.5e-40
% Rv151:CO2vc + O^- + CO => CO3^- + CO   k= 1.5e-40
% Rv152:CO2vd + O^- + CO => CO3^- + CO   k= 1.5e-40
% Rv153:CO2v1 + O^- + CO => CO3^- + CO   k= 1.5e-40
% Rv154:CO2v2 + O^- + CO => CO3^- + CO   k= 1.5e-40
% Rv155:CO2v3 + O^- + CO => CO3^- + CO   k= 1.5e-40
% Rv156:CO2v4 + O^- + CO => CO3^- + CO   k= 1.5e-40
% Rv157:CO2e1 + O^- + O2 => CO3^- + O2   k= 3.1e-40
% Rv158:CO2e2 + O^- + O2 => CO3^- + O2   k= 3.1e-40
% Rv159:CO2va + O^- + O2 => CO3^- + O2   k= 3.1e-40
% Rv160:CO2vb + O^- + O2 => CO3^- + O2   k= 3.1e-40
% Rv161:CO2vc + O^- + O2 => CO3^- + O2   k= 3.1e-40
% Rv162:CO2vd + O^- + O2 => CO3^- + O2   k= 3.1e-40
% Rv163:CO2v1 + O^- + O2 => CO3^- + O2   k= 3.1e-40
% Rv164:CO2v2 + O^- + O2 => CO3^- + O2   k= 3.1e-40
% Rv165:CO2v3 + O^- + O2 => CO3^- + O2   k= 3.1e-40
% Rv166:CO2v4 + O^- + O2 => CO3^- + O2   k= 3.1e-40
% Rv167:CO2va + CO^+ => CO2^+ + CO       k= 1.0e-15
% Rv168:CO2vb + CO^+ => CO2^+ + CO       k= 1.0e-15
% Rv169:CO2vc + CO^+ => CO2^+ + CO       k= 1.0e-15
% Rv170:CO2vd + CO^+ => CO2^+ + CO       k= 1.0e-15
% Rv171:CO2v1 + CO^+ => CO2^+ + CO       k= 1.0e-15
% Rv172:CO2v2 + CO^+ => CO2^+ + CO       k= 1.0e-15
% Rv173:CO2v3 + CO^+ => CO2^+ + CO       k= 1.0e-15
% Rv174:CO2v4 + CO^+ => CO2^+ + CO       k= 1.0e-15
% Rv175:O3^- + CO2e1 => CO3^- + O2       k= 5.5e-16
% Rv176:O3^- + CO2e2 => CO3^- + O2       k= 5.5e-16
% Rv177:O3^- + CO2va => CO3^- + O2       k= 5.5e-16
% Rv178:O3^- + CO2vb => CO3^- + O2       k= 5.5e-16
% Rv179:O3^- + CO2vc => CO3^- + O2       k= 5.5e-16
% Rv180:O3^- + CO2vd => CO3^- + O2       k= 5.5e-16
% Rv181:O3^- + CO2v1 => CO3^- + O2       k= 5.5e-16
% Rv182:O3^- + CO2v2 => CO3^- + O2       k= 5.5e-16
% Rv183:O3^- + CO2v3 => CO3^- + O2       k= 5.5e-16
% Rv184:O3^- + CO2v4 => CO3^- + O2       k= 5.5e-16
% Rv185:O^+ + CO2e1 => O2^+ + CO         k= 9.4e-16
% Rv186:O^+ + CO2e2 => O2^+ + CO         k= 9.4e-16
% Rv187:O^+ + CO2va => O2^+ + CO         k= 9.4e-16
% Rv188:O^+ + CO2vb => O2^+ + CO         k= 9.4e-16
% Rv189:O^+ + CO2vc => O2^+ + CO         k= 9.4e-16
% Rv190:O^+ + CO2vd => O2^+ + CO         k= 9.4e-16
% Rv191:O^+ + CO2v1 => O2^+ + CO         k= 9.4e-16
% Rv192:O^+ + CO2v2 => O2^+ + CO         k= 9.4e-16
% Rv193:O^+ + CO2v3 => O2^+ + CO         k= 9.4e-16
% Rv194:O^+ + CO2v4 => O2^+ + CO         k= 9.4e-16
% Rv195:O^+ + CO2va => O + CO2^+         k= 4.5e-16
% Rv196:O^+ + CO2vb => O + CO2^+         k= 4.5e-16
% Rv197:O^+ + CO2vc => O + CO2^+         k= 4.5e-16
% Rv198:O^+ + CO2vd => O + CO2^+         k= 4.5e-16
% Rv199:O^+ + CO2v1 => O + CO2^+         k= 4.5e-16
% Rv200:O^+ + CO2v2 => O + CO2^+         k= 4.5e-16
% Rv201:O^+ + CO2v3 => O + CO2^+         k= 4.5e-16
% Rv202:O^+ + CO2v4 => O + CO2^+         k= 4.5e-16
% Rv203:C^+ + CO2va => CO^+ + CO         k= 1.1e-15
% Rv204:C^+ + CO2vb => CO^+ + CO         k= 1.1e-15
% Rv205:C^+ + CO2vc => CO^+ + CO         k= 1.1e-15
% Rv206:C^+ + CO2vd => CO^+ + CO         k= 1.1e-15
% Rv207:C^+ + CO2v1 => CO^+ + CO         k= 1.1e-15
% Rv208:C^+ + CO2v2 => CO^+ + CO         k= 1.1e-15
% Rv209:C^+ + CO2v3 => CO^+ + CO         k= 1.1e-15
% Rv210:C^+ + CO2v4 => CO^+ + CO         k= 1.1e-15
% Rv211:O2^- + CO2e1 + O2 => CO4^- + O2  k= 4.7e-41
% Rv212:O2^- + CO2e2 + O2 => CO4^- + O2  k= 4.7e-41
% Rv213:O2^- + CO2va + O2 => CO4^- + O2  k= 4.7e-41
% Rv214:O2^- + CO2vb + O2 => CO4^- + O2  k= 4.7e-41
% Rv215:O2^- + CO2vc + O2 => CO4^- + O2  k= 4.7e-41
% Rv216:O2^- + CO2vd + O2 => CO4^- + O2  k= 4.7e-41
% Rv217:O2^- + CO2v1 + O2 => CO4^- + O2  k= 4.7e-41
% Rv218:O2^- + CO2v2 + O2 => CO4^- + O2  k= 4.7e-41
% Rv219:O2^- + CO2v3 + O2 => CO4^- + O2  k= 4.7e-41
% Rv220:O2^- + CO2v4 + O2 => CO4^- + O2  k= 4.7e-41
% Rv221:e + CO2e1 => CO + O + e          k= 查表
% Rv222:e + CO2e2 => CO + O + e          k= 查表
% Rv223:e + CO2e1 => CO + O^-            k= 查表
% Rv224:e + CO2e2 => CO + O^-            k= 查表
% Rv225:O + O2 + CO2e1 => O3 + CO2e1     k= 1.7e-42*(Tg*1.16e4)^(-1.2)
% Rv226:O + O2 + CO2e2 => O3 + CO2e2     k= 1.7e-42*(Tg*1.16e4)^(-1.2)
% Rv227:O + O2 + CO2va => O3 + CO2va     k= 1.7e-42*(Tg*1.16e4)^(-1.2)
% Rv228:O + O2 + CO2vb => O3 + CO2vb     k= 1.7e-42*(Tg*1.16e4)^(-1.2)
% Rv229:O + O2 + CO2vc => O3 + CO2vc     k= 1.7e-42*(Tg*1.16e4)^(-1.2)
% Rv230:O + O2 + CO2vd => O3 + CO2vd     k= 1.7e-42*(Tg*1.16e4)^(-1.2)
% Rv231:O + O2 + CO2v1 => O3 + CO2v1     k= 1.7e-42*(Tg*1.16e4)^(-1.2)
% Rv232:O + O2 + CO2v2 => O3 + CO2v2     k= 1.7e-42*(Tg*1.16e4)^(-1.2)
% Rv233:O + O2 + CO2v3 => O3 + CO2v3     k= 1.7e-42*(Tg*1.16e4)^(-1.2)
% Rv234:O + O2 + CO2v4 => O3 + CO2v4     k= 1.7e-42*(Tg*1.16e4)^(-1.2)
% Rv235:O^- + CO2e1 => O + CO2e1 + e    k= 4.0e-18
% Rv236:O^- + CO2e2 => O + CO2e2 + e    k= 4.0e-18
% Rv237:O^- + CO2va => O + CO2va + e    k= 4.0e-18
% Rv238:O^- + CO2vb => O + CO2vb + e    k= 4.0e-18
% Rv239:O^- + CO2vc => O + CO2vc + e    k= 4.0e-18
% Rv240:O^- + CO2vd => O + CO2vd + e    k= 4.0e-18
% Rv241:O^- + CO2v1 => O + CO2v1 + e    k= 4.0e-18
% Rv242:O^- + CO2v2 => O + CO2v2 + e    k= 4.0e-18
% Rv243:O^- + CO2v3 => O + CO2v3 + e    k= 4.0e-18
% Rv244:O^- + CO2v4 => O + CO2v4 + e    k= 4.0e-18
% Rv245:O^- + CO2 + CO2e1 => CO3^- + CO2e1 k= 9.0e-41
% Rv246:O^- + CO2 + CO2e2 => CO3^- + CO2e2 k= 9.0e-41
% Rv247:O^- + CO2 + CO2va => CO3^- + CO2va k= 9.0e-41
% Rv248:O^- + CO2 + CO2vb => CO3^- + CO2vb k= 9.0e-41
% Rv249:O^- + CO2 + CO2vc => CO3^- + CO2vc k= 9.0e-41
% Rv250:O^- + CO2 + CO2vd => CO3^- + CO2vd k= 9.0e-41
% Rv251:O^- + CO2 + CO2v1 => CO3^- + CO2v1 k= 9.0e-41
% Rv252:O^- + CO2 + CO2v2 => CO3^- + CO2v2 k= 9.0e-41
% Rv253:O^- + CO2 + CO2v3 => CO3^- + CO2v3 k= 9.0e-41
% Rv254:O^- + CO2 + CO2v4 => CO3^- + CO2v4 k= 9.0e-41
% Rv255:O + O + CO2e1 => O2 + CO2e1     k= 3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4))
% Rv256:O + O + CO2e2 => O2 + CO2e2     k= 3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4))
% Rv257:O + O + CO2va => O2 + CO2va     k= 3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4))
% Rv258:O + O + CO2vb => O2 + CO2vb     k= 3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4))
% Rv259:O + O + CO2vc => O2 + CO2vc     k= 3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4))
% Rv260:O + O + CO2vd => O2 + CO2vd     k= 3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4))
% Rv261:O + O + CO2v1 => O2 + CO2v1     k= 3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4))
% Rv262:O + O + CO2v2 => O2 + CO2v2     k= 3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4))
% Rv263:O + O + CO2v3 => O2 + CO2v3     k= 3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4))
% Rv264:O + O + CO2v4 => O2 + CO2v4     k= 3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4))
% Rv265:O2^+ + CO2e1 + CO2 => CO4^+ + CO2   k= 2.3e-41
% Rv266:O2^+ + CO2e2 + CO2 => CO4^+ + CO2   k= 2.3e-41
% Rv267:O2^+ + CO2va + CO2 => CO4^+ + CO2   k= 2.3e-41
% Rv268:O2^+ + CO2vb + CO2 => CO4^+ + CO2   k= 2.3e-41
% Rv269:O2^+ + CO2vc + CO2 => CO4^+ + CO2   k= 2.3e-41
% Rv270:O2^+ + CO2vd + CO2 => CO4^+ + CO2   k= 2.3e-41
% Rv271:O2^+ + CO2v1 + CO2 => CO4^+ + CO2   k= 2.3e-41
% Rv272:O2^+ + CO2v2 + CO2 => CO4^+ + CO2   k= 2.3e-41
% Rv273:O2^+ + CO2v3 + CO2 => CO4^+ + CO2   k= 2.3e-41
% Rv274:O2^+ + CO2v4 + CO2 => CO4^+ + CO2   k= 2.3e-41
% Rv275:e + O2^+ + CO2e1 => O2 + CO2e1     k=1.0e-38
% Rv276:e + O2^+ + CO2e2 => O2 + CO2e2     k=1.0e-38
% Rv277:e + O2^+ + CO2va => O2 + CO2va     k=1.0e-38
% Rv278:e + O2^+ + CO2vb => O2 + CO2vb     k=1.0e-38
% Rv279:e + O2^+ + CO2vc => O2 + CO2vc     k=1.0e-38
% Rv280:e + O2^+ + CO2vd => O2 + CO2vd     k=1.0e-38
% Rv281:e + O2^+ + CO2v1 => O2 + CO2v1     k=1.0e-38
% Rv282:e + O2^+ + CO2v2 => O2 + CO2v2     k=1.0e-38
% Rv283:e + O2^+ + CO2v3 => O2 + CO2v3     k=1.0e-38
% Rv284:e + O2^+ + CO2v4 => O2 + CO2v4     k=1.0e-38
% Rv285:O2^- + CO2e1 + CO2 => CO4^- + CO2   k= 1.0e-41
% Rv286:O2^- + CO2e2 + CO2 => CO4^- + CO2   k= 1.0e-41
% Rv287:O2^- + CO2va + CO2 => CO4^- + CO2   k= 1.0e-41
% Rv288:O2^- + CO2vb + CO2 => CO4^- + CO2   k= 1.0e-41
% Rv289:O2^- + CO2vc + CO2 => CO4^- + CO2   k= 1.0e-41
% Rv290:O2^- + CO2vd + CO2 => CO4^- + CO2   k= 1.0e-41
% Rv291:O2^- + CO2v1 + CO2 => CO4^- + CO2   k= 1.0e-41
% Rv292:O2^- + CO2v2 + CO2 => CO4^- + CO2   k= 1.0e-41
% Rv293:O2^- + CO2v3 + CO2 => CO4^- + CO2   k= 1.0e-41
% Rv294:O2^- + CO2v4 + CO2 => CO4^- + CO2   k= 1.0e-41
% Rv295:O2^- + CO2 + CO2e1 => CO4^- + CO2e1   k= 1.0e-41
% Rv296:O2^- + CO2 + CO2e2 => CO4^- + CO2e2   k= 1.0e-41
% Rv297:O2^- + CO2 + CO2va => CO4^- + CO2va   k= 1.0e-41
% Rv298:O2^- + CO2 + CO2vb => CO4^- + CO2vb   k= 1.0e-41
% Rv299:O2^- + CO2 + CO2vc => CO4^- + CO2vc   k= 1.0e-41
% Rv300:O2^- + CO2 + CO2vd => CO4^- + CO2vd   k= 1.0e-41
% Rv301:O2^- + CO2 + CO2v1 => CO4^- + CO2v1   k= 1.0e-41
% Rv302:O2^- + CO2 + CO2v2 => CO4^- + CO2v2   k= 1.0e-41
% Rv303:O2^- + CO2 + CO2v3 => CO4^- + CO2v3   k= 1.0e-41
% Rv304:O2^- + CO2 + CO2v4 => CO4^- + CO2v4   k= 1.0e-41
% Rv305:CO2vc + CO2 => CO2vb + CO2   k= 2.897e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv306:CO2vc + CO => CO2vb + CO     k= 2.028e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv307:CO2vc + O2 => CO2vb + O2     k= 2.028e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv308:CO2vd + CO2 => CO2vb + CO2   k= 1.582e-11*exp(-272*(Tg*1.16e4)^(-1/3)+437*(Tg*1.16e4)^(-2/3))
% Rv309:CO2vd + CO => CO2vb + CO     k= 4.745e-11*exp(-272*(Tg*1.16e4)^(-1/3)+437*(Tg*1.16e4)^(-2/3))
% Rv310:CO2vd + O2 => CO2vb + O2     k= 4.745e-11*exp(-272*(Tg*1.16e4)^(-1/3)+437*(Tg*1.16e4)^(-2/3))
% Rv311:CO2vd + CO2 => CO2vc + CO2   k= 4.321e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv312:CO2vd + CO => CO2vc + CO     k= 3.025e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv313:CO2vd + O2 => CO2vc + O2     k= 3.025e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv314:CO2vd + CO2 => CO2v1 + CO2   k= 1.775e-17*exp(-108*(Tg*1.16e4)^(-1/3)+165*(Tg*1.16e4)^(-2/3))
% Rv315:CO2vd + CO => CO2v1 + CO     k= 2.13e-17*exp(-108*(Tg*1.16e4)^(-1/3)+165*(Tg*1.16e4)^(-2/3))
% Rv316:CO2vd + O2 => CO2v1 + O2     k= 2.13e-17*exp(-108*(Tg*1.16e4)^(-1/3)+165*(Tg*1.16e4)^(-2/3))
% Rv317:CO2vd + CO2 => CO2vb + CO2vb k= 2.384e-15*exp(-89*(Tg*1.16e4)^(-1/3)+234*(Tg*1.16e4)^(-2/3))
% Rv318:CO2vd + CO2 => CO2va + CO2vc k= 6.48e-15*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3))
% Rv319:CO2vd + CO2va => CO2vb + CO2vc   k= 1.442e-14*exp(-88.2*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3))
% Rv320:CO2vd + CO2vb => CO2vc + CO2vc   k= 2.628e-14*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3))
% Rv321:CO2v1 + CO2v1 => CO2 + CO2v2 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv322:CO2v1 + CO2v2 => CO2 + CO2v3 k= 2.927e-17*exp(21.4*(Tg*1.16e4)^(-1/3)-53*(Tg*1.16e4)^(-2/3))
% Rv323:CO2v1 + CO2v3 => CO2 + CO2v4 k= 4.825e-17*exp(20.1*(Tg*1.16e4)^(-1/3)-65*(Tg*1.16e4)^(-2/3))
% Rv324:CO2v2 + CO2v2 => CO2v1 + CO2v3   k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv325:CO2v2 + CO2v3 => CO2v1 + CO2v4   k= 2.927e-17*exp(21.4*(Tg*1.16e4)^(-1/3)-53*(Tg*1.16e4)^(-2/3))
% Rv326:CO2v3 + CO2v3 => CO2v2 + CO2v4   k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv327:e + O2^+ + CO2v1a => O2 + CO2v1a k=1.0e-38
% Rv328:e + O2^+ + CO2v1b => O2 + CO2v1b k=1.0e-38
% Rv329:e + O2^+ + CO2v1c => O2 + CO2v1c k=1.0e-38
% Rv330:CO2v1 + CO2v4 => CO2v2 + CO2v3   k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv331:CO2v1 + CO2v3 => CO2v2 + CO2v2   k= Kv(330)
% Rv332:CO2v2 + CO2v4 => CO2v3 + CO2v3   k= Kv(330)
% Rv333:e + CO2v1 => e + CO2v2           k= 查表                                                EV=0.29eV
% Rv334:e + CO2v2 => e + CO2v3           k= 查表                                                EV=0.28eV
% Rv335:e + CO2v3 => e + CO2v4           k= 查表                                                EV=0.28eV
% Rv336:e + CO2 => e + CO2v1a            k= 查表                                                EV=0.339eV
% Rv337:e + CO2 => e + CO2v1b            k= 查表                                                EV=0.422eV
% Rv338:e + CO2 => e + CO2v1c            k= 查表                                                EV=0.505eV
% Rv339:e + CO2v1a => CO2^+ + 2e         k= K(1)                                                ER=13.8eV
% Rv340:e + CO2v1b => CO2^+ + 2e         k= K(1)                                                ER=13.8eV
% Rv341:e + CO2v1c => CO2^+ + 2e         k= K(1)                                                ER=13.8eV
% Rv342:e + CO2v1a => O2^+ + C + 2e      k= ?
% Rv343:e + CO2v1b => O2^+ + C + 2e      k= ?
% Rv344:e + CO2v1c => O2^+ + C + 2e      k= ?
% Rv345:e + CO2v1a => CO^+ + O + 2e      k= 查表                                                ER=19.161eV
% Rv346:e + CO2v1b => CO^+ + O + 2e      k= 查表                                                ER=19.078eV
% Rv347:e + CO2v1c => CO^+ + O + 2e      k= 查表                                                ER=18.995eV
% Rv348:e + CO2v1a => O^+ + CO + 2e      k= 查表                                                ER=18.761eV
% Rv349:e + CO2v1b => O^+ + CO + 2e      k= 查表                                                ER=18.678eV
% Rv350:e + CO2v1c => O^+ + CO + 2e      k= 查表                                                ER=18.595eV
% Rv351:e + CO2v1a => CO + O + e         k= 查表                                                ER=11.121eV
% Rv352:e + CO2v1b => CO + O + e         k= 查表                                                ER=11.038eV
% Rv353:e + CO2v1c => CO + O + e         k= 查表                                                ER=10.955eV
% Rv354:e + CO2v1a => CO + O^-           k= 查表
% Rv355:e + CO2v1b => CO + O^-           k= 查表
% Rv356:e + CO2v1c => CO + O^-           k= 查表
% Rv357:CO2v1a + CO2 => CO2vc + CO2      k= 8.568e-7*exp(-404*(Tg*1.16e4)^(-1/3)+1096*(Tg*1.16e4)^(-2/3))
% Rv358:CO2v1a + CO => CO2vc + CO        k= 2.57e-7*exp(-404*(Tg*1.16e4)^(-1/3)+1096*(Tg*1.16e4)^(-2/3))
% Rv359:CO2v1a + O2 => CO2vc + O2        k= 3.427e-7*exp(-404*(Tg*1.16e4)^(-1/3)+1096*(Tg*1.16e4)^(-2/3))
% Rv360:CO2v1a + CO2 => CO2vd + CO2      k= 1.431e-11*exp(-252*(Tg*1.16e4)^(-1/3)+685*(Tg*1.16e4)^(-2/3))
% Rv361:CO2v1a + CO => CO2vd + CO        k= 4.293e-12*exp(-252*(Tg*1.16e4)^(-1/3)+685*(Tg*1.16e4)^(-2/3))
% Rv362:CO2v1a + O2 => CO2vd + O2        k= 5.724e-12*exp(-252*(Tg*1.16e4)^(-1/3)+685*(Tg*1.16e4)^(-2/3))
% Rv363:CO2v1a + CO2 => CO2v1 + CO2      k= 7.143e-14*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv364:CO2v1a + CO => CO2v1 + CO        k= 5e-14*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv365:CO2v1a + O2 => CO2v1 + O2        k= 5e-14*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv366:CO2v1b + CO2 => CO2vc + CO2      k= 3.218e-15*exp(-137*(Tg*1.16e4)^(-1/3))
% Rv367:CO2v1b + CO => CO2vc + CO        k= 9.976e-15*exp(-137*(Tg*1.16e4)^(-1/3))
% Rv368:CO2v1b + O2 => CO2vc + O2        k= 9.976e-15*exp(-137*(Tg*1.16e4)^(-1/3))
% Rv369:CO2v1b + CO2 => CO2vd + CO2      k= 6.447e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv370:CO2v1b + CO => CO2vd + CO        k= 4.513e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv371:CO2v1b + O2 => CO2vd + O2        k= 4.513e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv372:CO2v1b + CO2 => CO2v1 + CO2      k= 1.071e-15*exp(-137*(Tg*1.16e4)^(-1/3))
% Rv373:CO2v1b + CO => CO2v1 + CO        k= 3.32e-15*exp(-137*(Tg*1.16e4)^(-1/3))
% Rv374:CO2v1b + O2 => CO2v1 + O2        k= 3.32e-15*exp(-137*(Tg*1.16e4)^(-1/3))
% Rv375:CO2v1b + CO2 => CO2v1a + CO2     k= 6.447e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv376:CO2v1b + CO => CO2v1a + CO       k= 4.513e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv377:CO2v1b + O2 => CO2v1a + O2       k= 4.513e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv378:CO2v1c + CO2 => CO2v1a + CO2     k= 1.071e-15*exp(-137*(Tg*1.16e4)^(-1/3))
% Rv379:CO2v1c + CO => CO2v1a + CO       k= 3.32e-15*exp(-137*(Tg*1.16e4)^(-1/3))
% Rv380:CO2v1c + O2 => CO2v1a + O2       k= 3.32e-15*exp(-137*(Tg*1.16e4)^(-1/3))
% Rv381:CO2v1c + CO2 => CO2v1b + CO2     k= 2.897e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv382:CO2v1c + CO => CO2v1b + CO       k= 2.028e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv383:CO2v1c + O2 => CO2v1b + O2       k= 2.028e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv384:CO2v2 + CO2 => CO2v1a + CO2      k= 8.568e-7*exp(-406*(Tg*1.16e4)^(-1/3)+829*(Tg*1.16e4)^(-2/3))
% Rv385:CO2v2 + CO => CO2v1a + CO        k= 2.57e-7*exp(-406*(Tg*1.16e4)^(-1/3)+829*(Tg*1.16e4)^(-2/3))
% Rv386:CO2v2 + O2 => CO2v1a + O2        k= 3.427e-7*exp(-406*(Tg*1.16e4)^(-1/3)+829*(Tg*1.16e4)^(-2/3))
% Rv387:CO2v2 + CO2 => CO2v1b + CO2      k= 1.725e-6*exp(-404*(Tg*1.16e4)^(-1/3)+1098*(Tg*1.16e4)^(-2/3))
% Rv388:CO2v2 + CO => CO2v1b + CO        k= 5.175e-7*exp(-404*(Tg*1.16e4)^(-1/3)+1098*(Tg*1.16e4)^(-2/3))
% Rv389:CO2v2 + O2 => CO2v1b + O2        k= 6.9e-7*exp(-404*(Tg*1.16e4)^(-1/3)+1098*(Tg*1.16e4)^(-2/3))
% Rv390:CO2v2 + CO2 => CO2v1c + CO2      k= 2.882e-11*exp(-253*(Tg*1.16e4)^(-1/3)+683*(Tg*1.16e4)^(-2/3))
% Rv391:CO2v2 + CO => CO2v1c + CO        k= 8.646e-12*exp(-253*(Tg*1.16e4)^(-1/3)+683*(Tg*1.16e4)^(-2/3))
% Rv392:CO2v2 + O2 => CO2v1c + O2        k= 1.153e-11*exp(-253*(Tg*1.16e4)^(-1/3)+683*(Tg*1.16e4)^(-2/3))
% Rv393:CO2v2 + CO2 => CO2v1a + CO2vb    k= 4.299e-11*exp(-241*(Tg*1.16e4)^(-1/3)+637*(Tg*1.16e4)^(-2/3))
% Rv394:CO2v2 + CO2 => CO2v1b + CO2va    k= 4.299e-11*exp(-241*(Tg*1.16e4)^(-1/3)+635*(Tg*1.16e4)^(-2/3))
% Rv395:CO2v1a + CO2 => CO2v1 + CO2va    k= 1.071e-15*exp(-87.7*(Tg*1.16e4)^(-1/3)+230*(Tg*1.16e4)^(-2/3))
% Rv396:CO2v1a + CO2va => CO2v1 + CO2vb  k= 2.157e-15*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3))
% Rv397:CO2v1a + CO2vb => CO2v1 + CO2vc  k= 4.344e-15*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3))
% Rv398:CO2v1a + CO2vc => CO2v1 + CO2vd  k= 6.48e-15*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3))
% Rv399:CO2v1b + CO2 => CO2va + CO2vd    k= 9.667e-15*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3))
% Rv400:CO2v1b + CO2va => CO2vb + CO2vd  k= 9.667e-15*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3))
% Rv401:CO2v1b + CO2vb => CO2vc + CO2vd  k= 4.332e-14*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3))
% Rv402:CO2v1b + CO2vc => CO2vd + CO2vd  k= 5.848e-14*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3))
% Rv403:CO2v1b + CO2v1 => CO2v1a + CO2v1a  k= 2.157e-15*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3))
% Rv404:CO2v1c + CO2 => CO2va + CO2vd    k= 2.378e-14*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3))
% Rv405:CO2v1c + CO2va => CO2vb + CO2vd  k= 5.292e-14*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3))
% Rv406:CO2v1c + CO2vb => CO2vc + CO2vd  k= 1.066e-13*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3))
% Rv407:CO2v1c + CO2vc => CO2vd + CO2vd  k= 1.438e-13*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3))
% Rv408:CO2v1c + CO2v1 => CO2v1a + CO2v1b  k= 5.305e-15*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3))
% Rv409:CO2v1c + CO2v1a => CO2v1b + CO2v1b k= 5.305e-15*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3))
% Rv410:CO2 + CO2va => CO + O + CO2va    k= 3.91e-16*exp(-48688/(Tg*11600))
% Rv411:CO2 + CO2vb => CO + O + CO2vb    k= 3.91e-16*exp(-47852/(Tg*11600))
% Rv412:CO2 + CO2vc => CO + O + CO2vc    k= 3.91e-16*exp(-47110/(Tg*11600))
% Rv413:CO2 + CO2vd => CO + O + CO2vd    k= 3.91e-16*exp(-46368/(Tg*11600))
% Rv414:CO2 + CO2v1 => CO + O + CO2v1    k= 3.91e-16*exp(-46739/(Tg*11600))
% Rv415:CO2 + CO2v2 => CO + O + CO2v2    k= 3.91e-16*exp(-44048/(Tg*11600))
% Rv416:CO2 + CO2v3 => CO + O + CO2v3    k= 3.91e-16*exp(-41449/(Tg*11600))
% Rv417:CO2 + CO2v4 => CO + O + CO2v4    k= 3.91e-16*exp(-38851/(Tg*11600))
% Rv418:CO2v1a + CO2 => CO + O + CO2     k= 3.91e-16*exp(-46284/(Tg*11600))
% Rv419:CO2v1b + CO2 => CO + O + CO2     k= 3.91e-16*exp(-45514/(Tg*11600))
% Rv420:CO2v1c + CO2 => CO + O + CO2     k= 3.91e-16*exp(-44744/(Tg*11600))
% Rv421:CO2 + CO2v1a => CO + O + CO2v1a  k= 3.91e-16*exp(-46284/(Tg*11600))
% Rv422:CO2 + CO2v1b => CO + O + CO2v1b  k= 3.91e-16*exp(-45514/(Tg*11600))
% Rv423:CO2 + CO2v1c => CO + O + CO2v1c  k= 3.91e-16*exp(-44744/(Tg*11600))
% Rv424:CO2v1a + O => CO + O2            k= 2.8e-17*exp(-24534/(Tg*11600))
% Rv425:CO2v1b + O => CO + O2            k= 2.8e-17*exp(-24052/(Tg*11600))
% Rv426:CO2v1c + O => CO + O2            k= 2.8e-17*exp(-23571/(Tg*11600))
% Rv427:CO2v1a + C => 2CO                k= 1.0e-21
% Rv428:CO2v1b + C => 2CO                k= 1.0e-21
% Rv429:CO2v1c + C => 2CO                k= 1.0e-21
% Rv430:O^- + CO2v1a + CO2 => CO3^- + CO2 k= 9.0e-41
% Rv431:O^- + CO2v1b + CO2 => CO3^- + CO2 k= 9.0e-41
% Rv432:O^- + CO2v1c + CO2 => CO3^- + CO2 k= 9.0e-41
% Rv433:O^- + CO2v1a + CO => CO3^- + CO   k= 1.5e-40
% Rv434:O^- + CO2v1b + CO => CO3^- + CO   k= 1.5e-40
% Rv435:O^- + CO2v1c + CO => CO3^- + CO   k= 1.5e-40
% Rv436:O^- + CO2v1a + O2 => CO3^- + O2   k= 3.1e-40
% Rv437:O^- + CO2v1b + O2 => CO3^- + O2   k= 3.1e-40
% Rv438:O^- + CO2v1c + O2 => CO3^- + O2   k= 3.1e-40
% Rv439:O^- + CO2 + CO2v1a => CO3^- + CO2v1a k= 9.0e-41
% Rv440:O^- + CO2 + CO2v1b => CO3^- + CO2v1b k= 9.0e-41
% Rv441:O^- + CO2 + CO2v1c => CO3^- + CO2v1c k= 9.0e-41
% Rv442:CO2v1a + CO^+ => CO2^+ + CO       k= 1.0e-15
% Rv443:CO2v1b + CO^+ => CO2^+ + CO       k= 1.0e-15
% Rv444:CO2v1c + CO^+ => CO2^+ + CO       k= 1.0e-15
% Rv445:O3^- + CO2v1a => CO3^- + O2       k= 5.5e-16
% Rv446:O3^- + CO2v1b => CO3^- + O2       k= 5.5e-16
% Rv447:O3^- + CO2v1c => CO3^- + O2       k= 5.5e-16
% Rv448:O^+ + CO2v1a => O2^+ + CO         k= 9.4e-16
% Rv449:O^+ + CO2v1b => O2^+ + CO         k= 9.4e-16
% Rv450:O^+ + CO2v1c => O2^+ + CO         k= 9.4e-16
% Rv451:O^+ + CO2v1a => O + CO2^+         k= 4.5e-16
% Rv452:O^+ + CO2v1b => O + CO2^+         k= 4.5e-16
% Rv453:O^+ + CO2v1c => O + CO2^+         k= 4.5e-16
% Rv454:C^+ + CO2v1a => CO^+ + CO         k= 1.1e-15
% Rv455:C^+ + CO2v1b => CO^+ + CO         k= 1.1e-15
% Rv456:C^+ + CO2v1c => CO^+ + CO         k= 1.1e-15
% Rv457:O2^+ + CO2v1a + CO2 => CO4^+ + CO2   k= 2.3e-41
% Rv458:O2^+ + CO2v1b + CO2 => CO4^+ + CO2   k= 2.3e-41
% Rv459:O2^+ + CO2v1c + CO2 => CO4^+ + CO2   k= 2.3e-41
% Rv460:O2^- + CO2v1a + CO2 => CO4^- + CO2   k= 1.0e-41
% Rv461:O2^- + CO2v1b + CO2 => CO4^- + CO2   k= 1.0e-41
% Rv462:O2^- + CO2v1c + CO2 => CO4^- + CO2   k= 1.0e-41
% Rv463:O2^- + CO2v1a + O2 => CO4^- + O2  k= 4.7e-41
% Rv464:O2^- + CO2v1b + O2 => CO4^- + O2  k= 4.7e-41
% Rv465:O2^- + CO2v1c + O2 => CO4^- + O2  k= 4.7e-41
% Rv466:O2^- + CO2 + CO2v1a => CO4^- + CO2v1a   k= 1.0e-41
% Rv467:O2^- + CO2 + CO2v1b => CO4^- + CO2v1b   k= 1.0e-41
% Rv468:O2^- + CO2 + CO2v1c => CO4^- + CO2v1c   k= 1.0e-41
% Rv469:O + O2 + CO2v1a => O3 + CO2v1a     k= 1.7e-42*(Tg*1.16e4)^(-1.2)
% Rv470:O + O2 + CO2v1b => O3 + CO2v1b     k= 1.7e-42*(Tg*1.16e4)^(-1.2)
% Rv471:O + O2 + CO2v1c => O3 + CO2v1c     k= 1.7e-42*(Tg*1.16e4)^(-1.2)
% Rv472:O^- + CO2v1a => O + CO2v1a + e    k= 4.0e-18
% Rv473:O^- + CO2v1b => O + CO2v1b + e    k= 4.0e-18
% Rv474:O^- + CO2v1c => O + CO2v1c + e    k= 4.0e-18
% Rv475:O + O + CO2v1a => O2 + CO2v1a     k= 3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4))
% Rv476:O + O + CO2v1b => O2 + CO2v1b     k= 3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4))
% Rv477:O + O + CO2v1c => O2 + CO2v1c     k= 3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4))
% Rv478:e + O2 + CO2e1 => O2^- + CO2e1      k= 2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4))
% Rv479:e + O2 + CO2e2 => O2^- + CO2e2      k= 2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4))
% Rv480:e + O2 + CO2va => O2^- + CO2va      k= 2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4))
% Rv481:e + O2 + CO2vb => O2^- + CO2vb      k= 2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4))
% Rv482:e + O2 + CO2vc => O2^- + CO2vc      k= 2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4))
% Rv483:e + O2 + CO2vd => O2^- + CO2vd      k= 2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4))
% Rv484:e + O2 + CO2v1 => O2^- + CO2v1      k= 2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4))
% Rv485:e + O2 + CO2v2 => O2^- + CO2v2      k= 2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4))
% Rv486:e + O2 + CO2v3 => O2^- + CO2v3      k= 2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4))
% Rv487:e + O2 + CO2v4 => O2^- + CO2v4      k= 2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4))
% Rv488:e + O2 + CO2v1a => O2^- + CO2v1a    k= 2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4))
% Rv489:e + O2 + CO2v1b => O2^- + CO2v1b    k= 2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4))
% Rv490:e + O2 + CO2v1c => O2^- + CO2v1c    k= 2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4))
% Rv491:e + CO2 => e + CO2v5           k= 查表                                            EV=1.43eV
% Rv492:e + CO2 => e + CO2v6           k= 查表                                            EV=1.70eV
% Rv493:e + CO2 => e + CO2v7           k= 查表                                            EV=1.97eV
% Rv494:e + CO2 => e + CO2v8           k= 查表                                            EV=2.24eV
% Rv495:e + CO2v4 => e + CO2v5         k= 查表                                            EV=0.29eV
% Rv496:e + CO2v5 => e + CO2v6         k= 查表                                            EV=0.27eV
% Rv497:e + CO2v6 => e + CO2v7         k= 查表                                            EV=0.27eV
% Rv498:e + CO2v7 => e + CO2v8         k= 查表                                            EV=0.27eV
% Rv499:e + CO2v5 => CO2^+ + 2e        k= K(1)                                             ER=13.8eV
% Rv500:e + CO2v6 => CO2^+ + 2e        k= K(1)                                             ER=13.8eV
% Rv501:e + CO2v7 => CO2^+ + 2e        k= K(1)                                             ER=13.8eV
% Rv502:e + CO2v8 => CO2^+ + 2e        k= K(1)                                             ER=13.8eV
% Rv503:e + CO2v5 => O2^+ + C + 2e     k= ?                                    
% Rv504:e + CO2v6 => O2^+ + C + 2e     k= ?                                   
% Rv505:e + CO2v7 => O2^+ + C + 2e     k= ?                                    
% Rv506:e + CO2v8 => O2^+ + C + 2e     k= ?                                    
% Rv507:e + CO2v5 => 2e + O + CO^+     k= 查表                                            ER=18.07eV                        
% Rv508:e + CO2v6 => 2e + O + CO^+     k= 查表                                            ER=17.8eV                   
% Rv509:e + CO2v7 => 2e + O + CO^+     k= 查表                                            ER=17.53eV                                 
% Rv510:e + CO2v8 => 2e + O + CO^+     k= 查表                                            ER=17.26eV
% Rv511:e + CO2v5 => 2e + CO + O^+     k= 查表                                            ER=17.67eV
% Rv512:e + CO2v6 => 2e + CO + O^+     k= 查表                                            ER=17.4eV
% Rv513:e + CO2v7 => 2e + CO + O^+     k= 查表                                            ER=17.13eV
% Rv514:e + CO2v8 => 2e + CO + O^+     k= 查表                                            ER=16.86eV
% Rv515:e + CO2v5 => CO + O + e        k= 查表                                            ER=10.03eV
% Rv516:e + CO2v6 => CO + O + e        k= 查表                                            ER=9.76eV
% Rv517:e + CO2v7 => CO + O + e        k= 查表                                            ER=9.49eV
% Rv518:e + CO2v8 => CO + O + e        k= 查表                                            ER=9.22eV
% Rv519:e + CO2v5 => CO + O-           k= 查表                                        附着过程无ER，下同
% Rv520:e + CO2v6 => CO + O-           k= 查表
% Rv521:e + CO2v7 => CO + O-           k= 查表                    
% Rv522:e + CO2v8 => CO + O-           k= 查表
% Rv523:CO2v5 + CO2 => CO2v4 + CO2     k= 见下文
% Rv524:CO2v5 + CO => CO2v4 + CO       k= 见下文
% Rv525:CO2v5 + O2 => CO2v4 + O2       k= 见下文
% Rv526:CO2v6 + CO2 => CO2v5 + CO2     k= 见下文
% Rv527:CO2v6 + CO => CO2v5 + CO       k= 见下文
% Rv528:CO2v6 + O2 => CO2v5 + O2       k= 见下文
% Rv529:CO2v7 + CO2 => CO2v6 + CO2     k= 见下文
% Rv530:CO2v7 + CO => CO2v6 + CO       k= 见下文
% Rv531:CO2v7 + O2 => CO2v6 + O2       k= 见下文
% Rv532:CO2v8 + CO2 => CO2v7 + CO2     k= 见下文
% Rv533:CO2v8 + CO => CO2v7 + CO       k= 见下文
% Rv534:CO2v8 + O2 => CO2v7 + O2       k= 见下文
% Rv535:CO2v5 + CO2 => CO2v4 + CO2va   k= 2.13e-11*exp(-242*(Tg*11600)^(-1/3)+633*(Tg*11600)^(-2/3))
% Rv536:CO2v5 + CO2 => CO2v4 + CO2vb   k= 2.13e-11*exp(-242*(Tg*11600)^(-1/3)+633*(Tg*11600)^(-2/3))
% Rv537:CO2v6 + CO2 => CO2v5 + CO2va   k= 2.13e-11*exp(-242*(Tg*11600)^(-1/3)+633*(Tg*11600)^(-2/3))
% Rv538:CO2v6 + CO2 => CO2v5 + CO2vb   k= 2.13e-11*exp(-242*(Tg*11600)^(-1/3)+633*(Tg*11600)^(-2/3))
% Rv539:CO2v7 + CO2 => CO2v6 + CO2va   k= 2.13e-11*exp(-242*(Tg*11600)^(-1/3)+633*(Tg*11600)^(-2/3))
% Rv540:CO2v7 + CO2 => CO2v6 + CO2vb   k= 2.13e-11*exp(-242*(Tg*11600)^(-1/3)+633*(Tg*11600)^(-2/3))
% Rv541:CO2v8 + CO2 => CO2v7 + CO2va   k= 2.13e-11*exp(-242*(Tg*11600)^(-1/3)+633*(Tg*11600)^(-2/3))
% Rv542:CO2v8 + CO2 => CO2v7 + CO2vb   k= 2.13e-11*exp(-242*(Tg*11600)^(-1/3)+633*(Tg*11600)^(-2/3))
% Rv543:CO2v1 + CO2v4 => CO2 + CO2v5   k= 7.199e-17*exp(18.4*(Tg*1.16e4)^(-1/3)-76*(Tg*1.16e4)^(-2/3))
% Rv544:CO2v2 + CO2v4 => CO2v1 + CO2v5 k= 4.825e-17*exp(20.1*(Tg*1.16e4)^(-1/3)-65*(Tg*1.16e4)^(-2/3))
% Rv545:CO2v2 + CO2v5 => CO2v1 + CO2v6 k= 7.199e-17*exp(18.4*(Tg*1.16e4)^(-1/3)-76*(Tg*1.16e4)^(-2/3))
% Rv546:CO2v3 + CO2v4 => CO2v2 + CO2v5 k= 2.927e-17*exp(21.4*(Tg*1.16e4)^(-1/3)-53*(Tg*1.16e4)^(-2/3))
% Rv547:CO2v3 + CO2v5 => CO2v2 + CO2v6 k= 4.825e-17*exp(20.1*(Tg*1.16e4)^(-1/3)-65*(Tg*1.16e4)^(-2/3))
% Rv548:CO2v3 + CO2v6 => CO2v2 + CO2v7 k= 7.199e-17*exp(18.4*(Tg*1.16e4)^(-1/3)-76*(Tg*1.16e4)^(-2/3))
% Rv549:CO2v4 + CO2v4 => CO2v3 + CO2v5 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv550:CO2v4 + CO2v5 => CO2v3 + CO2v6 k= 2.927e-17*exp(21.4*(Tg*1.16e4)^(-1/3)-53*(Tg*1.16e4)^(-2/3))
% Rv551:CO2v4 + CO2v6 => CO2v3 + CO2v7 k= 4.825e-17*exp(20.1*(Tg*1.16e4)^(-1/3)-65*(Tg*1.16e4)^(-2/3))
% Rv552:CO2v4 + CO2v7 => CO2v3 + CO2v8 k= 7.199e-17*exp(18.4*(Tg*1.16e4)^(-1/3)-76*(Tg*1.16e4)^(-2/3))
% Rv553:CO2v5 + CO2v5 => CO2v4 + CO2v6 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv554:CO2v5 + CO2v6 => CO2v4 + CO2v7 k= 2.927e-17*exp(21.4*(Tg*1.16e4)^(-1/3)-53*(Tg*1.16e4)^(-2/3))
% Rv555:CO2v5 + CO2v7 => CO2v4 + CO2v8 k= 4.825e-17*exp(20.1*(Tg*1.16e4)^(-1/3)-65*(Tg*1.16e4)^(-2/3))
% Rv556:CO2v6 + CO2v6 => CO2v5 + CO2v7 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv557:CO2v6 + CO2v7 => CO2v5 + CO2v8 k= 2.927e-17*exp(21.4*(Tg*1.16e4)^(-1/3)-53*(Tg*1.16e4)^(-2/3))
% Rv558:CO2v7 + CO2v7 => CO2v6 + CO2v8 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv559:CO2v1 + CO2v5 => CO2v2 + CO2v4 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv560:CO2v1 + CO2v6 => CO2v2 + CO2v5 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv561:CO2v1 + CO2v7 => CO2v2 + CO2v6 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv562:CO2v1 + CO2v8 => CO2v2 + CO2v7 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv563:CO2v2 + CO2v5 => CO2v3 + CO2v4 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv564:CO2v2 + CO2v6 => CO2v3 + CO2v5 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv565:CO2v2 + CO2v7 => CO2v3 + CO2v6 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv566:CO2v2 + CO2v8 => CO2v3 + CO2v7 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv567:CO2v3 + CO2v5 => CO2v4 + CO2v4 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv568:CO2v3 + CO2v6 => CO2v4 + CO2v5 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv569:CO2v3 + CO2v7 => CO2v4 + CO2v6 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv570:CO2v3 + CO2v8 => CO2v4 + CO2v7 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv571:CO2v4 + CO2v6 => CO2v5 + CO2v5 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv572:CO2v4 + CO2v7 => CO2v5 + CO2v6 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv573:CO2v4 + CO2v8 => CO2v5 + CO2v7 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv574:CO2v5 + CO2v7 => CO2v6 + CO2v6 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv575:CO2v5 + CO2v8 => CO2v6 + CO2v7 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv576:CO2v6 + CO2v8 => CO2v7 + CO2v7 k= 1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3))
% Rv577:CO2 + CO2v5 => CO + O + CO2v5    k= 3.91e-16*exp(-36160/(Tg*11600))
% Rv578:CO2 + CO2v6 => CO + O + CO2v6    k= 3.91e-16*exp(-33654/(Tg*11600))
% Rv579:CO2 + CO2v7 => CO + O + CO2v7    k= 3.91e-16*exp(-31148/(Tg*11600))
% Rv580:CO2 + CO2v8 => CO + O + CO2v8    k= 3.91e-16*exp(-28643/(Tg*11600))
% Rv581:CO2 + CO2v5 => CO + O + CO2      k= 3.91e-16*exp(-36160/(Tg*11600))
% Rv582:CO2 + CO2v6 => CO + O + CO2      k= 3.91e-16*exp(-33654/(Tg*11600))
% Rv583:CO2 + CO2v7 => CO + O + CO2      k= 3.91e-16*exp(-31148/(Tg*11600))
% Rv584:CO2 + CO2v8 => CO + O + CO2      k= 3.91e-16*exp(-28643/(Tg*11600))
% Rv585:CO2v5 + O => CO + O2             k= 2.8e-17*exp(-18206/(Tg*11600))
% Rv586:CO2v6 + O => CO + O2             k= 2.8e-17*exp(-16640/(Tg*11600))
% Rv587:CO2v7 + O => CO + O2             k= 2.8e-17*exp(-15074/(Tg*11600))
% Rv588:CO2v8 + O => CO + O2             k= 2.8e-17*exp(-13508/(Tg*11600))
% Rv589:CO2v5 + C => 2CO                 k= 1.0e-21
% Rv590:CO2v6 + C => 2CO                 k= 1.0e-21
% Rv591:CO2v7 + C => 2CO                 k= 1.0e-21
% Rv592:CO2v8 + C => 2CO                 k= 1.0e-21
% Rv593:O^- + CO2v5 + CO2 => CO3^- + CO2 k= 9.0e-41
% Rv594:O^- + CO2v6 + CO2 => CO3^- + CO2 k= 9.0e-41
% Rv595:O^- + CO2v7 + CO2 => CO3^- + CO2 k= 9.0e-41
% Rv596:O^- + CO2v8 + CO2 => CO3^- + CO2 k= 9.0e-41
% Rv597:O^- + CO2v5 + CO => CO3^- + CO   k= 1.5e-40
% Rv598:O^- + CO2v6 + CO => CO3^- + CO   k= 1.5e-40
% Rv599:O^- + CO2v7 + CO => CO3^- + CO   k= 1.5e-40
% Rv600:O^- + CO2v8 + CO => CO3^- + CO   k= 1.5e-40
% Rv601:O^- + CO2v5 + O2 => CO3^- + O2   k= 3.1e-40
% Rv602:O^- + CO2v6 + O2 => CO3^- + O2   k= 3.1e-40
% Rv603:O^- + CO2v7 + O2 => CO3^- + O2   k= 3.1e-40
% Rv604:O^- + CO2v8 + O2 => CO3^- + O2   k= 3.1e-40
% Rv605:O^- + CO2 + CO2v5 => CO3^- + CO2v5 k= 9.0e-41
% Rv606:O^- + CO2 + CO2v6 => CO3^- + CO2v6 k= 9.0e-41
% Rv607:O^- + CO2 + CO2v7 => CO3^- + CO2v7 k= 9.0e-41
% Rv608:O^- + CO2 + CO2v8 => CO3^- + CO2v8 k= 9.0e-41
% Rv609:e + O2^+ + CO2v5 => O2 + CO2v5   k= 1.0e-38
% Rv610:e + O2^+ + CO2v6 => O2 + CO2v6   k= 1.0e-38
% Rv611:e + O2^+ + CO2v7 => O2 + CO2v7   k= 1.0e-38
% Rv612:e + O2^+ + CO2v8 => O2 + CO2v8   k= 1.0e-38
% Rv613:CO2v5 + CO^+ => CO2^+ + CO       k= 1.0e-15
% Rv614:CO2v6 + CO^+ => CO2^+ + CO       k= 1.0e-15
% Rv615:CO2v7 + CO^+ => CO2^+ + CO       k= 1.0e-15
% Rv616:CO2v8 + CO^+ => CO2^+ + CO       k= 1.0e-15
% Rv617:O3^- + CO2v5 => CO3^- + O2       k= 5.5e-16
% Rv618:O3^- + CO2v6 => CO3^- + O2       k= 5.5e-16
% Rv619:O3^- + CO2v7 => CO3^- + O2       k= 5.5e-16
% Rv620:O3^- + CO2v8 => CO3^- + O2       k= 5.5e-16
% Rv621:O^+ + CO2v5 => O2^+ + CO         k= 9.4e-16
% Rv622:O^+ + CO2v6 => O2^+ + CO         k= 9.4e-16
% Rv623:O^+ + CO2v7 => O2^+ + CO         k= 9.4e-16
% Rv624:O^+ + CO2v8 => O2^+ + CO         k= 9.4e-16
% Rv625:O^+ + CO2v5 => O + CO2^+         k= 4.5e-16
% Rv626:O^+ + CO2v6 => O + CO2^+         k= 4.5e-16
% Rv627:O^+ + CO2v7 => O + CO2^+         k= 4.5e-16
% Rv628:O^+ + CO2v8 => O + CO2^+         k= 4.5e-16
% Rv629:C^+ + CO2v5 => CO^+ + CO         k= 1.1e-15
% Rv630:C^+ + CO2v6 => CO^+ + CO         k= 1.1e-15
% Rv631:C^+ + CO2v7 => CO^+ + CO         k= 1.1e-15
% Rv632:C^+ + CO2v8 => CO^+ + CO         k= 1.1e-15
% Rv633:O2^+ + CO2v5 + CO2 => CO4^+ + CO2   k= 2.3e-41
% Rv634:O2^+ + CO2v6 + CO2 => CO4^+ + CO2   k= 2.3e-41
% Rv635:O2^+ + CO2v7 + CO2 => CO4^+ + CO2   k= 2.3e-41
% Rv636:O2^+ + CO2v8 + CO2 => CO4^+ + CO2   k= 2.3e-41
% Rv637:O2^- + CO2v5 + CO2 => CO4^- + CO2   k= 1.0e-41
% Rv638:O2^- + CO2v6 + CO2 => CO4^- + CO2   k= 1.0e-41
% Rv639:O2^- + CO2v7 + CO2 => CO4^- + CO2   k= 1.0e-41
% Rv640:O2^- + CO2v8 + CO2 => CO4^- + CO2   k= 1.0e-41
% Rv641:O2^- + CO2v5 + O2 => CO4^- + O2  k= 4.7e-41
% Rv642:O2^- + CO2v6 + O2 => CO4^- + O2  k= 4.7e-41
% Rv643:O2^- + CO2v7 + O2 => CO4^- + O2  k= 4.7e-41
% Rv644:O2^- + CO2v8 + O2 => CO4^- + O2  k= 4.7e-41
% Rv645:O2^- + CO2 + CO2v5 => CO4^- + CO2v5   k= 1.0e-41
% Rv646:O2^- + CO2 + CO2v6 => CO4^- + CO2v6   k= 1.0e-41
% Rv647:O2^- + CO2 + CO2v7 => CO4^- + CO2v7   k= 1.0e-41
% Rv648:O2^- + CO2 + CO2v8 => CO4^- + CO2v8   k= 1.0e-41
% Rv649:O + O2 + CO2v5 => O3 + CO2v5     k= 1.7e-42*(Tg*1.16e4)^(-1.2)
% Rv650:O + O2 + CO2v6 => O3 + CO2v6     k= 1.7e-42*(Tg*1.16e4)^(-1.2)
% Rv651:O + O2 + CO2v7 => O3 + CO2v7     k= 1.7e-42*(Tg*1.16e4)^(-1.2)
% Rv652:O + O2 + CO2v8 => O3 + CO2v8     k= 1.7e-42*(Tg*1.16e4)^(-1.2)
% Rv653:O^- + CO2v5 => O + CO2v5 + e     k= 4.0e-18
% Rv654:O^- + CO2v6 => O + CO2v6 + e     k= 4.0e-18
% Rv655:O^- + CO2v7 => O + CO2v7 + e     k= 4.0e-18
% Rv656:O^- + CO2v8 => O + CO2v8 + e     k= 4.0e-18
% Rv657:O + O + CO2v5 => O2 + CO2v5     k= 3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4))
% Rv658:O + O + CO2v6 => O2 + CO2v6     k= 3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4))
% Rv659:O + O + CO2v7 => O2 + CO2v7     k= 3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4))
% Rv660:O + O + CO2v8 => O2 + CO2v8     k= 3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4))
% Rv661:e + O2 + CO2v5 => O2^- + CO2v5      k= 2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4))
% Rv662:e + O2 + CO2v6 => O2^- + CO2v6      k= 2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4))
% Rv663:e + O2 + CO2v7 => O2^- + CO2v7      k= 2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4))
% Rv664:e + O2 + CO2v8 => O2^- + CO2v8      k= 2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4))
% Rv665:e + CO2 => e + CO2vk                k= 查表                                         EV=2.5eV
% Rv666:e + CO2vk => CO2^+ + 2e             k= K(1)                                         ER=13.8eV
% Rv667:e + CO2vk => O2^+ + C + 2e          k= ?
% Rv668:e + CO2vk => 2e + O + CO^+          k= 查表                                         ER=17.0eV
% Rv669:e + CO2vk => 2e + CO + O^+          k= 查表                                         ER=16.6eV
% Rv670:e + CO2vk => CO + O + e             k= 查表                                         ER=8.96eV
% Rv671:e + CO2vk => CO + O-                k= 查表
% Rv672:CO2v1 + CO2 => CO2vk + CO2      k= 1.43e-11*exp(-252*(Tg*11600)^(-1/3)+685*(Tg*11600)^(-2/3))
% Rv673:CO2v1 + CO => CO2vk + CO        k= 4.29e-12*exp(-252*(Tg*11600)^(-1/3)+685*(Tg*11600)^(-2/3))
% Rv674:CO2v1 + O2 => CO2vk + O2        k= 5.72e-12*exp(-252*(Tg*11600)^(-1/3)+685*(Tg*11600)^(-2/3))
% Rv675:CO2vk + CO2 => CO2va + CO2      k= 1.071e-15*exp(-137*(Tg*1.16e4)^(-1/3))
% Rv676:CO2vk + CO => CO2va + CO        k= 3.213e-15*exp(-137*(Tg*1.16e4)^(-1/3))
% Rv677:CO2vk + O2 => CO2va + O2        k= 3.213e-15*exp(-137*(Tg*1.16e4)^(-1/3))
% Rv678:CO2vk + CO2 => CO2vb + CO2      k= 2.897e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv679:CO2vk + CO => CO2vb + CO        k= 2.028e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv680:CO2vk + O2 => CO2vb + O2        k= 2.028e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3))
% Rv681:CO2 + CO2vk => CO + O + CO2vk   k= 3.91e-16*exp(-26230/(Tg*11600))
% Rv682:CO2 + CO2vk => CO + O + CO2     k= 3.91e-16*exp(-26230/(Tg*11600))
% Rv683:CO2vk + O => CO + O2            k= 2.8e-17*exp(-12000/(Tg*11600))
% Rv684:CO2vk + C => 2CO                k= 1.0e-21
% Rv685:O^- + CO2vk + CO2 => CO3^- + CO2 k= 9.0e-41
% Rv686:O^- + CO2vk + CO => CO3^- + CO   k= 1.5e-40
% Rv687:O^- + CO2vk + O2 => CO3^- + O2   k= 3.1e-40
% Rv688:O^- + CO2 + CO2vk => CO3^- + CO2vk k= 9.0e-41
% Rv689:e + O2^+ + CO2vk => O2 + CO2vk   k= 1.0e-38
% Rv690:CO2vk + CO^+ => CO2^+ + CO       k= 1.0e-15
% Rv691:O3^- + CO2vk => CO3^- + O2       k= 5.5e-16
% Rv692:O^+ + CO2vk => O2^+ + CO         k= 9.4e-16
% Rv693:O^+ + CO2vk => O + CO2^+         k= 4.5e-16
% Rv694:C^+ + CO2vk => CO^+ + CO         k= 1.1e-15
% Rv695:O2^+ + CO2vk + CO2 => CO4^+ + CO2   k= 2.3e-41
% Rv696:O2^- + CO2vk + CO2 => CO4^- + CO2   k= 1.0e-41
% Rv697:O2^- + CO2vk + O2 => CO4^- + O2  k= 4.7e-41
% Rv698:O2^- + CO2 + CO2vk => CO4^- + CO2vk   k= 1.0e-41
% Rv699:O + O2 + CO2vk => O3 + CO2vk     k= 1.7e-42*(Tg*1.16e4)^(-1.2)
% Rv700:O^- + CO2vk => O + CO2vk + e     k= 4.0e-18
% Rv701:O + O + CO2vk => O2 + CO2vk      k= 3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4))
% Rv702:e + O2 + CO2vk => O2^- + CO2vk   k= 2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4))
% Rv703: CO2vk + CO2 => CO2va + CO2vb    k= 5.305e-15*exp(-88*(Tg*11600)^(-1/3)+233*(Tg*11600)^(-2/3))
% Rv704:CO2^+ + CO2e1 + CO => C2O4^+ + CO  k= 3.0e-40
% Rv705:CO2^+ + CO2e2 + CO => C2O4^+ + CO  k= 3.0e-40
% Rv706:CO2^+ + CO2va + CO => C2O4^+ + CO  k= 3.0e-40
% Rv707:CO2^+ + CO2vb + CO => C2O4^+ + CO  k= 3.0e-40
% Rv708:CO2^+ + CO2vc + CO => C2O4^+ + CO  k= 3.0e-40
% Rv709:CO2^+ + CO2vd + CO => C2O4^+ + CO  k= 3.0e-40
% Rv710:CO2^+ + CO2v1 + CO => C2O4^+ + CO  k= 3.0e-40
% Rv711:CO2^+ + CO2v2 + CO => C2O4^+ + CO  k= 3.0e-40
% Rv712:CO2^+ + CO2v3 + CO => C2O4^+ + CO  k= 3.0e-40
% Rv713:CO2^+ + CO2v4 + CO => C2O4^+ + CO  k= 3.0e-40
% Rv714:CO2^+ + CO2v1a + CO => C2O4^+ + CO  k= 3.0e-40
% Rv715:CO2^+ + CO2v1b + CO => C2O4^+ + CO  k= 3.0e-40
% Rv716:CO2^+ + CO2v1c + CO => C2O4^+ + CO  k= 3.0e-40
% Rv717:CO2^+ + CO2v5 + CO => C2O4^+ + CO  k= 3.0e-40
% Rv718:CO2^+ + CO2v6 + CO => C2O4^+ + CO  k= 3.0e-40
% Rv719:CO2^+ + CO2v7 + CO => C2O4^+ + CO  k= 3.0e-40
% Rv720:CO2^+ + CO2v8 + CO => C2O4^+ + CO  k= 3.0e-40
% Rv721:CO2^+ + CO2vk + CO => C2O4^+ + CO  k= 3.0e-40
% Rv722:C2O3^+ + CO + CO2e1 => CO2 + C2O2^+ + CO2e1  k= 2.6e-38
% Rv723:C2O3^+ + CO + CO2e2 => CO2 + C2O2^+ + CO2e2  k= 2.6e-38
% Rv724:C2O3^+ + CO + CO2va => CO2 + C2O2^+ + CO2va  k= 2.6e-38
% Rv725:C2O3^+ + CO + CO2vb => CO2 + C2O2^+ + CO2vb  k= 2.6e-38
% Rv726:C2O3^+ + CO + CO2vc => CO2 + C2O2^+ + CO2vc  k= 2.6e-38
% Rv727:C2O3^+ + CO + CO2vd => CO2 + C2O2^+ + CO2vd  k= 2.6e-38
% Rv728:C2O3^+ + CO + CO2v1 => CO2 + C2O2^+ + CO2v1  k= 2.6e-38
% Rv729:C2O3^+ + CO + CO2v2 => CO2 + C2O2^+ + CO2v2  k= 2.6e-38
% Rv730:C2O3^+ + CO + CO2v3 => CO2 + C2O2^+ + CO2v3  k= 2.6e-38
% Rv731:C2O3^+ + CO + CO2v4 => CO2 + C2O2^+ + CO2v4  k= 2.6e-38
% Rv732:C2O3^+ + CO + CO2v1a => CO2 + C2O2^+ + CO2v1a  k= 2.6e-38
% Rv733:C2O3^+ + CO + CO2v1b => CO2 + C2O2^+ + CO2v1b  k= 2.6e-38
% Rv734:C2O3^+ + CO + CO2v1c => CO2 + C2O2^+ + CO2v1c  k= 2.6e-38
% Rv735:C2O3^+ + CO + CO2v5 => CO2 + C2O2^+ + CO2v5  k= 2.6e-38
% Rv736:C2O3^+ + CO + CO2v6 => CO2 + C2O2^+ + CO2v6  k= 2.6e-38
% Rv737:C2O3^+ + CO + CO2v7 => CO2 + C2O2^+ + CO2v7  k= 2.6e-38
% Rv738:C2O3^+ + CO + CO2v8 => CO2 + C2O2^+ + CO2v8  k= 2.6e-38
% Rv739:C2O3^+ + CO + CO2vk => CO2 + C2O2^+ + CO2vk  k= 2.6e-38
% Rv740:C2O4^+ + CO + CO2e1 => CO2 + C2O3^+ + CO2e1  k= 4.2e-38
% Rv741:C2O4^+ + CO + CO2e2 => CO2 + C2O3^+ + CO2e2  k= 4.2e-38
% Rv742:C2O4^+ + CO + CO2va => CO2 + C2O3^+ + CO2va  k= 4.2e-38
% Rv743:C2O4^+ + CO + CO2vb => CO2 + C2O3^+ + CO2vb  k= 4.2e-38
% Rv744:C2O4^+ + CO + CO2vc => CO2 + C2O3^+ + CO2vc  k= 4.2e-38
% Rv745:C2O4^+ + CO + CO2vd => CO2 + C2O3^+ + CO2vd  k= 4.2e-38
% Rv746:C2O4^+ + CO + CO2v1 => CO2 + C2O3^+ + CO2v1  k= 4.2e-38
% Rv747:C2O4^+ + CO + CO2v2 => CO2 + C2O3^+ + CO2v2  k= 4.2e-38
% Rv748:C2O4^+ + CO + CO2v3 => CO2 + C2O3^+ + CO2v3  k= 4.2e-38
% Rv749:C2O4^+ + CO + CO2v4 => CO2 + C2O3^+ + CO2v4  k= 4.2e-38
% Rv750:C2O4^+ + CO + CO2v1a => CO2 + C2O3^+ + CO2v1a  k= 4.2e-38
% Rv751:C2O4^+ + CO + CO2v1b => CO2 + C2O3^+ + CO2v1b  k= 4.2e-38
% Rv752:C2O4^+ + CO + CO2v1c => CO2 + C2O3^+ + CO2v1c  k= 4.2e-38
% Rv753:C2O4^+ + CO + CO2v5 => CO2 + C2O3^+ + CO2v5  k= 4.2e-38
% Rv754:C2O4^+ + CO + CO2v6 => CO2 + C2O3^+ + CO2v6  k= 4.2e-38
% Rv755:C2O4^+ + CO + CO2v7 => CO2 + C2O3^+ + CO2v7  k= 4.2e-38
% Rv756:C2O4^+ + CO + CO2v8 => CO2 + C2O3^+ + CO2v8  k= 4.2e-38
% Rv757:C2O4^+ + CO + CO2vk => CO2 + C2O3^+ + CO2vk  k= 4.2e-38
% Rv758:C2O2^+ + CO2e1 => CO^+ + CO + CO2e1    k= 1.0e-18
% Rv759:C2O2^+ + CO2e2 => CO^+ + CO + CO2e2    k= 1.0e-18
% Rv760:C2O2^+ + CO2va => CO^+ + CO + CO2va    k= 1.0e-18
% Rv761:C2O2^+ + CO2vb => CO^+ + CO + CO2vb    k= 1.0e-18
% Rv762:C2O2^+ + CO2vc => CO^+ + CO + CO2vc    k= 1.0e-18
% Rv763:C2O2^+ + CO2vd => CO^+ + CO + CO2vd    k= 1.0e-18
% Rv764:C2O2^+ + CO2v1 => CO^+ + CO + CO2v1    k= 1.0e-18
% Rv765:C2O2^+ + CO2v2 => CO^+ + CO + CO2v2    k= 1.0e-18
% Rv766:C2O2^+ + CO2v3 => CO^+ + CO + CO2v3    k= 1.0e-18
% Rv767:C2O2^+ + CO2v4 => CO^+ + CO + CO2v4    k= 1.0e-18
% Rv768:C2O2^+ + CO2v1a => CO^+ + CO + CO2v1a    k= 1.0e-18
% Rv769:C2O2^+ + CO2v1b => CO^+ + CO + CO2v1b    k= 1.0e-18
% Rv770:C2O2^+ + CO2v1c => CO^+ + CO + CO2v1c    k= 1.0e-18
% Rv771:C2O2^+ + CO2v5 => CO^+ + CO + CO2v5    k= 1.0e-18
% Rv772:C2O2^+ + CO2v6 => CO^+ + CO + CO2v6    k= 1.0e-18
% Rv773:C2O2^+ + CO2v7 => CO^+ + CO + CO2v7    k= 1.0e-18
% Rv774:C2O2^+ + CO2v8 => CO^+ + CO + CO2v8    k= 1.0e-18
% Rv775:C2O2^+ + CO2vk => CO^+ + CO + CO2vk    k= 1.0e-18
% Rv776:C2O4^+ + CO2e1 => CO2^+ + CO2 + CO2e1  k= 1.0e-20
% Rv777:C2O4^+ + CO2e2 => CO2^+ + CO2 + CO2e2  k= 1.0e-20
% Rv778:C2O4^+ + CO2va => CO2^+ + CO2 + CO2va  k= 1.0e-20
% Rv779:C2O4^+ + CO2vb => CO2^+ + CO2 + CO2vb  k= 1.0e-20
% Rv780:C2O4^+ + CO2vc => CO2^+ + CO2 + CO2vc  k= 1.0e-20
% Rv781:C2O4^+ + CO2vd => CO2^+ + CO2 + CO2vd  k= 1.0e-20
% Rv782:C2O4^+ + CO2v1 => CO2^+ + CO2 + CO2v1  k= 1.0e-20
% Rv783:C2O4^+ + CO2v2 => CO2^+ + CO2 + CO2v2  k= 1.0e-20
% Rv784:C2O4^+ + CO2v3 => CO2^+ + CO2 + CO2v3  k= 1.0e-20
% Rv785:C2O4^+ + CO2v4 => CO2^+ + CO2 + CO2v4  k= 1.0e-20
% Rv786:C2O4^+ + CO2v1a => CO2^+ + CO2 + CO2v1a  k= 1.0e-20
% Rv787:C2O4^+ + CO2v1b => CO2^+ + CO2 + CO2v1b  k= 1.0e-20
% Rv788:C2O4^+ + CO2v1c => CO2^+ + CO2 + CO2v1c  k= 1.0e-20
% Rv789:C2O4^+ + CO2v5 => CO2^+ + CO2 + CO2v5  k= 1.0e-20
% Rv790:C2O4^+ + CO2v6 => CO2^+ + CO2 + CO2v6  k= 1.0e-20
% Rv791:C2O4^+ + CO2v7 => CO2^+ + CO2 + CO2v7  k= 1.0e-20
% Rv792:C2O4^+ + CO2v8 => CO2^+ + CO2 + CO2v8  k= 1.0e-20
% Rv793:C2O4^+ + CO2vk => CO2^+ + CO2 + CO2vk  k= 1.0e-20
% Rv794:N2^+ + CO2e1 => CO2^+ + N2       k= K(125)*0.916
% Rv795:N2^+ + CO2e2 => CO2^+ + N2       k= K(125)*0.093
% Rv796:N2^+ + CO2va => CO2^+ + N2       k= 9.0e-16
% Rv797:N2^+ + CO2vb => CO2^+ + N2       k= 9.0e-16
% Rv798:N2^+ + CO2vc => CO2^+ + N2       k= 9.0e-16
% Rv799:N2^+ + CO2vd => CO2^+ + N2       k= 9.0e-16
% Rv800:N2^+ + CO2v1 => CO2^+ + N2       k= 9.0e-16
% Rv801:N2^+ + CO2v2 => CO2^+ + N2       k= 9.0e-16
% Rv802:N2^+ + CO2v3 => CO2^+ + N2       k= 9.0e-16
% Rv803:N2^+ + CO2v4 => CO2^+ + N2       k= 9.0e-16
% Rv804:N2^+ + CO2v1a => CO2^+ + N2       k= 9.0e-16
% Rv805:N2^+ + CO2v1b => CO2^+ + N2       k= 9.0e-16
% Rv806:N2^+ + CO2v1c => CO2^+ + N2       k= 9.0e-16
% Rv807:N2^+ + CO2v5 => CO2^+ + N2       k= 9.0e-16
% Rv808:N2^+ + CO2v6 => CO2^+ + N2       k= 9.0e-16
% Rv809:N2^+ + CO2v7 => CO2^+ + N2       k= 9.0e-16
% Rv810:N2^+ + CO2v8 => CO2^+ + N2       k= 9.0e-16
% Rv811:N2^+ + CO2vk => CO2^+ + N2       k= 9.0e-16
% Rv812:N + CO2e1 => NO + CO             k= 1.0e-25
% Rv813:N + CO2e2 => NO + CO             k= 1.0e-25
% Rv814:N + CO2va => NO + CO             k= 1.0e-25
% Rv815:N + CO2vb => NO + CO             k= 1.0e-25
% Rv816:N + CO2vc => NO + CO             k= 1.0e-25
% Rv817:N + CO2vd => NO + CO             k= 1.0e-25
% Rv818:N + CO2v1 => NO + CO             k= 1.0e-25
% Rv819:N + CO2v2 => NO + CO             k= 1.0e-25
% Rv820:N + CO2v3 => NO + CO             k= 1.0e-25
% Rv821:N + CO2v4 => NO + CO             k= 1.0e-25
% Rv822:N + CO2v1a => NO + CO             k= 1.0e-25
% Rv823:N + CO2v1b => NO + CO             k= 1.0e-25
% Rv824:N + CO2v1c => NO + CO             k= 1.0e-25
% Rv825:N + CO2v5 => NO + CO             k= 1.0e-25
% Rv826:N + CO2v6 => NO + CO             k= 1.0e-25
% Rv827:N + CO2v7 => NO + CO             k= 1.0e-25
% Rv828:N + CO2v8 => NO + CO             k= 1.0e-25
% Rv829:N + CO2vk => NO + CO             k= 1.0e-25
% Rv830:NO + O + CO2e1 => CO2e1 + NO2    k= 1.13e-37*(Tg*1.16e4)^(-2.16)*exp(-529/(Tg*1.16e4))
% Rv831:NO + O + CO2e2 => CO2e2 + NO2    k= 1.13e-37*(Tg*1.16e4)^(-2.16)*exp(-529/(Tg*1.16e4))
% Rv832:NO + O + CO2va => CO2va + NO2    k= 1.13e-37*(Tg*1.16e4)^(-2.16)*exp(-529/(Tg*1.16e4))
% Rv833:NO + O + CO2vb => CO2vb + NO2    k= 1.13e-37*(Tg*1.16e4)^(-2.16)*exp(-529/(Tg*1.16e4))
% Rv834:NO + O + CO2vc => CO2vc + NO2    k= 1.13e-37*(Tg*1.16e4)^(-2.16)*exp(-529/(Tg*1.16e4))
% Rv835:NO + O + CO2vd => CO2vd + NO2    k= 1.13e-37*(Tg*1.16e4)^(-2.16)*exp(-529/(Tg*1.16e4))
% Rv836:NO + O + CO2v1 => CO2v1 + NO2    k= 1.13e-37*(Tg*1.16e4)^(-2.16)*exp(-529/(Tg*1.16e4))
% Rv837:NO + O + CO2v2 => CO2v2 + NO2    k= 1.13e-37*(Tg*1.16e4)^(-2.16)*exp(-529/(Tg*1.16e4))
% Rv838:NO + O + CO2v3 => CO2v3 + NO2    k= 1.13e-37*(Tg*1.16e4)^(-2.16)*exp(-529/(Tg*1.16e4))
% Rv839:NO + O + CO2v4 => CO2v4 + NO2    k= 1.13e-37*(Tg*1.16e4)^(-2.16)*exp(-529/(Tg*1.16e4))
% Rv840:NO + O + CO2v1a => CO2v1a + NO2    k= 1.13e-37*(Tg*1.16e4)^(-2.16)*exp(-529/(Tg*1.16e4))
% Rv841:NO + O + CO2v1b => CO2v1b + NO2    k= 1.13e-37*(Tg*1.16e4)^(-2.16)*exp(-529/(Tg*1.16e4))
% Rv842:NO + O + CO2v1c => CO2v1c + NO2    k= 1.13e-37*(Tg*1.16e4)^(-2.16)*exp(-529/(Tg*1.16e4))
% Rv843:NO + O + CO2v5 => CO2v5 + NO2    k= 1.13e-37*(Tg*1.16e4)^(-2.16)*exp(-529/(Tg*1.16e4))
% Rv844:NO + O + CO2v6 => CO2v6 + NO2    k= 1.13e-37*(Tg*1.16e4)^(-2.16)*exp(-529/(Tg*1.16e4))
% Rv845:NO + O + CO2v7 => CO2v7 + NO2    k= 1.13e-37*(Tg*1.16e4)^(-2.16)*exp(-529/(Tg*1.16e4))
% Rv846:NO + O + CO2v8 => CO2v8 + NO2    k= 1.13e-37*(Tg*1.16e4)^(-2.16)*exp(-529/(Tg*1.16e4))
% Rv847:NO + O + CO2vk => CO2vk + NO2    k= 1.13e-37*(Tg*1.16e4)^(-2.16)*exp(-529/(Tg*1.16e4))
% Rv848:NO^- + CO2e1 => NO + CO2e1 + e   k= 8.3e-18
% Rv849:NO^- + CO2e2 => NO + CO2e2 + e   k= 8.3e-18
% Rv850:NO^- + CO2va => NO + CO2va + e   k= 8.3e-18
% Rv851:NO^- + CO2vb => NO + CO2vb + e   k= 8.3e-18
% Rv852:NO^- + CO2vc => NO + CO2vc + e   k= 8.3e-18
% Rv853:NO^- + CO2vd => NO + CO2vd + e   k= 8.3e-18
% Rv854:NO^- + CO2v1 => NO + CO2v1 + e   k= 8.3e-18
% Rv855:NO^- + CO2v2 => NO + CO2v2 + e   k= 8.3e-18
% Rv856:NO^- + CO2v3 => NO + CO2v3 + e   k= 8.3e-18
% Rv857:NO^- + CO2v4 => NO + CO2v4 + e   k= 8.3e-18
% Rv858:NO^- + CO2v1a => NO + CO2v1a + e   k= 8.3e-18
% Rv859:NO^- + CO2v1b => NO + CO2v1b + e   k= 8.3e-18
% Rv860:NO^- + CO2v1c => NO + CO2v1c + e   k= 8.3e-18
% Rv861:NO^- + CO2v5 => NO + CO2v5 + e   k= 8.3e-18
% Rv862:NO^- + CO2v6 => NO + CO2v6 + e   k= 8.3e-18
% Rv863:NO^- + CO2v7 => NO + CO2v7 + e   k= 8.3e-18
% Rv864:NO^- + CO2v8 => NO + CO2v8 + e   k= 8.3e-18
% Rv865:NO^- + CO2vk => NO + CO2vk + e   k= 8.3e-18

Ralldlt=0;
R(1:number_of_reaction)=0;
Rv=zeros(1,number_of_vreaction);
X=zeros(nsp,number_of_reactions+1);
Xgv=zeros(length(n),number_of_vreaction);
Xvv=zeros(length(nv),number_of_vreaction);
% Rate=zeros(1,8);
% RateK=zeros(1,8);

Ngg=n*n';
Ngv=n*nv';
Nvv=nv*nv';

%% Function of Te
%基态存在电子得失
% K(1)=lookups('reaction-rate-cross_section-e_CO2ionization.lut',Te);
% K(3)=lookups('reaction-rate-cross_section-e_CO2attachment.lut',Te);
% K(6)=lookups('reaction-rate-cross_section-e_O2attachment.lut',Te);
K(1)=interp1(aux{1}(:,1),aux{1}(:,2),Te);
K(3)=interp1(aux{3}(:,1),aux{3}(:,2),Te);
K(6)=interp1(aux{6}(:,1),aux{6}(:,2),Te);
K(9)=2.0e-11*Te^(-0.5)*(Tg*1.16e4)^(-1);
% K(17)=lookups('reaction-rate-cross_section-e_COattachment.lut',Te);
K(17)=interp1(aux{17}(:,1),aux{17}(:,2),Te);
K(19)=3.92e-13*Te^(-0.4);
% K(22)=lookups('reaction-rate-cross_section-e_O2ionization.lut',Te);
K(22)=interp1(aux{22}(:,1),aux{22}(:,2),Te);
K(23)=6.0e-13*(1/Te)^0.5*(1/(Tg*1.16e4))^0.5;
% K(32)=K(1)*1/13;
K(32)=0;
K(33)=3.2e-17*(Te*1.16e4)^(0.5)*(1+0.15*Te)*exp(-12.93/Te);
% K(34)=lookups('reaction-rate-cross_section-e_O3attachment_O-.lut',Te);
% K(37)=lookups('reaction-rate-cross_section-e_O3attachment_O2-.lut',Te);
% K(44)=lookups('reaction-rate-cross_section-e_CO2ionization_CO+.lut',Te);
% K(45)=lookups('reaction-rate-cross_section-e_COionization.lut',Te);
% K(53)=lookups('reaction-rate-cross_section-e_CO2ionization_O+.lut',Te);
% K(54)=lookups('reaction-rate-cross_section-e_Oionization.lut',Te);
% K(58)=lookups('reaction-rate-cross_section-e_COionization_C+.lut',Te);
% K(59)=lookups('reaction-rate-cross_section-e_Cionization.lut',Te);
K(34)=interp1(aux{34}(:,1),aux{34}(:,2),Te);
K(37)=interp1(aux{37}(:,1),aux{37}(:,2),Te);
K(44)=interp1(aux{44}(:,1),aux{44}(:,2),Te);
K(45)=interp1(aux{45}(:,1),aux{45}(:,2),Te);
K(53)=interp1(aux{53}(:,1),aux{53}(:,2),Te);
K(54)=interp1(aux{54}(:,1),aux{54}(:,2),Te);
K(58)=interp1(aux{58}(:,1),aux{58}(:,2),Te);
K(59)=interp1(aux{59}(:,1),aux{59}(:,2),Te);
%K(69)=1.61e-13*Te^(-0.5);
K(69)=0;
K(74)=4.0e-13*Te^(-0.34);
K(75)=5.4e-14*Te^(-0.7);
K(76)=2.0e-11*(Te^(-0.5))*(Tg*1.16e4)^(-1);
K(95)=interp1(aux{95}(:,1),aux{95}(:,2),Te);
K(98)=interp1(aux{98}(:,1),aux{98}(:,2),Te);
K(100)=interp1(aux{100}(:,1),aux{100}(:,2),Te);
K(101)=2.6e-15*Te^1/2*exp(-10/Te);
K(102)=8.1e-15*Te^1/2*exp(-12.9/Te);
K(117)=7.0e-14*Te^(-0.5);
K(118)=4.6e-14*Te^(-0.39);
K(119)=7.0e-14*Te^(-0.5);
K(120)=7.0e-14*Te^(-0.5);

%基态无电子得失
% K(2)=lookups('reaction-rate-cross_section-e_CO2dissociation.lut',Te);
% K(5)=lookups('reaction-rate-cross_section-e_O2dissociation.lut',Te);
% K(13)=lookups('reaction-rate-cross_section-e_CO2elastic.lut',Te);%e与CO2发生弹性碰撞的速率系数
% K(18)=lookups('reaction-rate-cross_section-e_COexcitation.lut',Te);
K(2)=interp1(aux{2}(:,1),aux{2}(:,2),Te);
K(5)=interp1(aux{5}(:,1),aux{5}(:,2),Te);
K(13)=interp1(aux{13}(:,1),aux{13}(:,2),Te);
K(18)=interp1(aux{18}(:,1),aux{18}(:,2),Te);
K(94)=interp1(aux{94}(:,1),aux{94}(:,2),Te);
K(96)=interp1(aux{96}(:,1),aux{96}(:,2),Te);
K(97)=interp1(aux{97}(:,1),aux{97}(:,2),Te);
K(99)=interp1(aux{99}(:,1),aux{99}(:,2),Te);
K(103)=7.4e-15*exp(-6.5/Te) ;
K(104)=5.6e-15*exp(-3.11/Te);
K(105)=1.4e-15*exp(-1.67/Te);

% 振动态
if Vibration==1
%     K(73)=lookups('reaction-rate-cross_section-CO2e1.lut',Te);
%     K(74)=lookups('reaction-rate-cross_section-CO2e2.lut',Te);
    Kv(1)=interp1(auxv{1}(:,1),auxv{1}(:,2),Te);
    Kv(2)=interp1(auxv{2}(:,1),auxv{2}(:,2),Te);
    Kv(3)=K(1)*0.944;
    Kv(4)=K(1)*0.099;
    Kv(5)=K(44)*3.191;
    Kv(6)=K(44)*0.735;
    Kv(9)=K(53)*2.988;
    Kv(10)=K(53)*0.671;
%     K(87)=lookups('reaction-rate-cross_section-CO2va.lut',Te);
%     K(88)=lookups('reaction-rate-cross_section-CO2vb.lut',Te);
%     K(89)=lookups('reaction-rate-cross_section-CO2vc.lut',Te);
%     K(90)=lookups('reaction-rate-cross_section-CO2vd.lut',Te);
%     K(103)=lookups('reaction-rate-cross_section-CO2v1.lut',Te);
%     K(104)=lookups('reaction-rate-cross_section-CO2v2.lut',Te);
%     K(105)=lookups('reaction-rate-cross_section-CO2v3.lut',Te);
%     K(106)=lookups('reaction-rate-cross_section-CO2v4.lut',Te);
    Kv(15)=interp1(auxv{15}(:,1),auxv{15}(:,2),Te);
    Kv(16)=interp1(auxv{16}(:,1),auxv{16}(:,2),Te);
    Kv(17)=interp1(auxv{17}(:,1),auxv{17}(:,2),Te);
    Kv(18)=interp1(auxv{18}(:,1),auxv{18}(:,2),Te);
    Kv(31)=interp1(auxv{31}(:,1),auxv{31}(:,2),Te);
    Kv(32)=interp1(auxv{32}(:,1),auxv{32}(:,2),Te);
    Kv(33)=interp1(auxv{33}(:,1),auxv{33}(:,2),Te);
    Kv(34)=interp1(auxv{34}(:,1),auxv{34}(:,2),Te);
    Kv(35)=K(1);
    Kv(36)=K(1);
    Kv(37)=K(1);
    Kv(38)=K(1);
%     K(111)=lookups('reaction-rate-cross_section-CO2v1dissociation.lut',Te);
%     K(112)=lookups('reaction-rate-cross_section-CO2v2dissociation.lut',Te);
%     K(113)=lookups('reaction-rate-cross_section-CO2v3dissociation.lut',Te);
%     K(114)=lookups('reaction-rate-cross_section-CO2v4dissociation.lut',Te);
%     K(115)=lookups('reaction-rate-cross_section-CO2v1attachment.lut',Te);
%     K(116)=lookups('reaction-rate-cross_section-CO2v2attachment.lut',Te);
    Kv(39)=interp1(auxv{39}(:,1),auxv{39}(:,2),Te);
    Kv(40)=interp1(auxv{40}(:,1),auxv{40}(:,2),Te);
    Kv(41)=interp1(auxv{41}(:,1),auxv{41}(:,2),Te);
    Kv(42)=interp1(auxv{42}(:,1),auxv{42}(:,2),Te);
    Kv(43)=interp1(auxv{43}(:,1),auxv{43}(:,2),Te);
    Kv(44)=interp1(auxv{44}(:,1),auxv{44}(:,2),Te);
    Kv(45)=interp1(auxv{45}(:,1),auxv{45}(:,2),Te);
    Kv(46)=interp1(auxv{46}(:,1),auxv{46}(:,2),Te);
%     K(119)=K(1)*1/13;
%     K(120)=K(1)*1/13;
%     K(121)=K(1)*1/13;
%     K(122)=K(1)*1/13;
    Kv(47)=0;
    Kv(48)=0;
    Kv(49)=0;
    Kv(50)=0;
%     K(123)=lookups('reaction-rate-cross_section-CO2v1ionization_CO+.lut',Te);
%     K(124)=lookups('reaction-rate-cross_section-CO2v2ionization_CO+.lut',Te);
%     K(125)=lookups('reaction-rate-cross_section-CO2v3ionization_CO+.lut',Te);
%     K(126)=lookups('reaction-rate-cross_section-CO2v4ionization_CO+.lut',Te);
%     K(127)=lookups('reaction-rate-cross_section-CO2v1ionization_O+.lut',Te);
%     K(128)=lookups('reaction-rate-cross_section-CO2v2ionization_O+.lut',Te);
%     K(129)=lookups('reaction-rate-cross_section-CO2v3ionization_O+.lut',Te);
%     K(130)=lookups('reaction-rate-cross_section-CO2v4ionization_O+.lut',Te);
    Kv(51)=interp1(auxv{51}(:,1),auxv{51}(:,2),Te);
    Kv(52)=interp1(auxv{52}(:,1),auxv{52}(:,2),Te);
    Kv(53)=interp1(auxv{53}(:,1),auxv{53}(:,2),Te);
    Kv(54)=interp1(auxv{54}(:,1),auxv{54}(:,2),Te);
    Kv(55)=interp1(auxv{55}(:,1),auxv{55}(:,2),Te);
    Kv(56)=interp1(auxv{56}(:,1),auxv{56}(:,2),Te);
    Kv(57)=interp1(auxv{57}(:,1),auxv{57}(:,2),Te);
    Kv(58)=interp1(auxv{58}(:,1),auxv{58}(:,2),Te);
    %对称拉伸及弯曲模式补充
%     K(175)=lookups('reaction-rate-cross_section-CO2vadissociation.lut',Te);
%     K(176)=lookups('reaction-rate-cross_section-CO2vbdissociation.lut',Te);
%     K(177)=lookups('reaction-rate-cross_section-CO2vcdissociation.lut',Te);
%     K(178)=lookups( 'reaction-rate-cross_section-CO2vddissociation.lut',Te);
%     K(179)=lookups('reaction-rate-cross_section-CO2vaattachment.lut',Te);
%     K(180)=lookups('reaction-rate-cross_section-CO2vbattachment.lut',Te);
%     K(181)=lookups('reaction-rate-cross_section-CO2vcattachment.lut',Te);
%     K(182)=lookups('reaction-rate-cross_section-CO2vdattachment.lut',Te);
    Kv(103)=interp1(auxv{103}(:,1),auxv{103}(:,2),Te);
    Kv(104)=interp1(auxv{104}(:,1),auxv{104}(:,2),Te);
    Kv(105)=interp1(auxv{105}(:,1),auxv{105}(:,2),Te);
    Kv(106)=interp1(auxv{106}(:,1),auxv{106}(:,2),Te);
    Kv(107)=interp1(auxv{107}(:,1),auxv{107}(:,2),Te);
    Kv(108)=interp1(auxv{108}(:,1),auxv{108}(:,2),Te);
    Kv(109)=interp1(auxv{109}(:,1),auxv{109}(:,2),Te);
    Kv(110)=interp1(auxv{110}(:,1),auxv{110}(:,2),Te);
    Kv(111)=K(1);
    Kv(112)=K(1);
    Kv(113)=K(1);
    Kv(114)=K(1);
%     K(187)=K(1)*1/13;
%     K(188)=K(1)*1/13;
%     K(189)=K(1)*1/13;
%     K(190)=K(1)*1/13; 
    Kv(115)=0;
    Kv(116)=0;
    Kv(117)=0;
    Kv(118)=0;
%     K(191)=lookups('reaction-rate-cross_section-CO2vaionization_CO+.lut',Te);
%     K(192)=lookups('reaction-rate-cross_section-CO2vbionization_CO+.lut',Te);
%     K(193)=lookups('reaction-rate-cross_section-CO2vcionization_CO+.lut',Te);
%     K(194)=lookups('reaction-rate-cross_section-CO2vdionization_CO+.lut',Te);
%     K(195)=lookups('reaction-rate-cross_section-CO2vaionization_O+.lut',Te);
%     K(196)=lookups('reaction-rate-cross_section-CO2vbionization_O+.lut',Te);
%     K(197)=lookups('reaction-rate-cross_section-CO2vcionization_O+.lut',Te);
%     K(198)=lookups('reaction-rate-cross_section-CO2vdionization_O+.lut',Te);
    Kv(119)=interp1(auxv{119}(:,1),auxv{119}(:,2),Te);
    Kv(120)=interp1(auxv{120}(:,1),auxv{120}(:,2),Te);
    Kv(121)=interp1(auxv{121}(:,1),auxv{121}(:,2),Te);
    Kv(122)=interp1(auxv{122}(:,1),auxv{122}(:,2),Te);
    Kv(123)=interp1(auxv{123}(:,1),auxv{123}(:,2),Te);
    Kv(124)=interp1(auxv{124}(:,1),auxv{124}(:,2),Te);
    Kv(125)=interp1(auxv{125}(:,1),auxv{125}(:,2),Te);
    Kv(126)=interp1(auxv{126}(:,1),auxv{126}(:,2),Te);
    Kv(221)=K(2);
    Kv(222)=K(2);
    Kv(223)=K(3);
    Kv(224)=K(3);
    Kv(333)=Kv(31);
    Kv(334)=interp1(auxv{334}(:,1),auxv{334}(:,2),Te);
    Kv(335)=Kv(334);
    Kv(336)=interp1(auxv{336}(:,1),auxv{336}(:,2),Te);
    Kv(337)=interp1(auxv{337}(:,1),auxv{337}(:,2),Te);
    Kv(338)=interp1(auxv{338}(:,1),auxv{338}(:,2),Te);
    Kv(339)=K(1);
    Kv(340)=K(1);
    Kv(341)=K(1);
    Kv(342)=0;
    Kv(343)=0;
    Kv(344)=0;
    Kv(345)=interp1(auxv{345}(:,1),auxv{345}(:,2),Te);
    Kv(346)=interp1(auxv{346}(:,1),auxv{346}(:,2),Te);
    Kv(347)=interp1(auxv{347}(:,1),auxv{347}(:,2),Te);
    Kv(348)=interp1(auxv{348}(:,1),auxv{348}(:,2),Te);
    Kv(349)=interp1(auxv{349}(:,1),auxv{349}(:,2),Te);
    Kv(350)=interp1(auxv{350}(:,1),auxv{350}(:,2),Te);
    Kv(351)=interp1(auxv{351}(:,1),auxv{351}(:,2),Te);
    Kv(352)=interp1(auxv{352}(:,1),auxv{352}(:,2),Te);
    Kv(353)=interp1(auxv{353}(:,1),auxv{353}(:,2),Te);
    Kv(354)=interp1(auxv{354}(:,1),auxv{354}(:,2),Te);
    Kv(355)=interp1(auxv{355}(:,1),auxv{355}(:,2),Te);
    Kv(356)=interp1(auxv{356}(:,1),auxv{356}(:,2),Te);
    Kv(491)=interp1(auxv{491}(:,1),auxv{491}(:,2),Te);
    Kv(492)=interp1(auxv{492}(:,1),auxv{492}(:,2),Te);
    Kv(493)=interp1(auxv{493}(:,1),auxv{493}(:,2),Te);
    Kv(494)=interp1(auxv{494}(:,1),auxv{494}(:,2),Te);
    Kv(495)=Kv(31);
    Kv(496)=interp1(auxv{496}(:,1),auxv{496}(:,2),Te);
    Kv(497:498)=Kv(496);
    Kv(499:502)=K(1);
    Kv(503:506)=0;
    Kv(507)=interp1(auxv{507}(:,1),auxv{507}(:,2),Te);
    Kv(508)=interp1(auxv{508}(:,1),auxv{508}(:,2),Te);
    Kv(509)=interp1(auxv{509}(:,1),auxv{509}(:,2),Te);
    Kv(510)=interp1(auxv{510}(:,1),auxv{510}(:,2),Te);
    Kv(511)=interp1(auxv{511}(:,1),auxv{511}(:,2),Te);
    Kv(512)=interp1(auxv{512}(:,1),auxv{512}(:,2),Te);
    Kv(513)=interp1(auxv{513}(:,1),auxv{513}(:,2),Te);
    Kv(514)=interp1(auxv{514}(:,1),auxv{514}(:,2),Te);
    Kv(515)=interp1(auxv{515}(:,1),auxv{515}(:,2),Te);
    Kv(516)=interp1(auxv{516}(:,1),auxv{516}(:,2),Te);
    Kv(517)=interp1(auxv{517}(:,1),auxv{517}(:,2),Te);
    Kv(518)=interp1(auxv{518}(:,1),auxv{518}(:,2),Te);
    Kv(519)=interp1(auxv{519}(:,1),auxv{519}(:,2),Te);
    Kv(520)=interp1(auxv{520}(:,1),auxv{520}(:,2),Te);
    Kv(521)=interp1(auxv{521}(:,1),auxv{521}(:,2),Te);
    Kv(522)=interp1(auxv{522}(:,1),auxv{522}(:,2),Te);
    Kv(665)=interp1(auxv{665}(:,1),auxv{665}(:,2),Te);
    Kv(666)=K(1);
    Kv(667)=0;
    Kv(668)=interp1(auxv{668}(:,1),auxv{668}(:,2),Te);
    Kv(669)=interp1(auxv{669}(:,1),auxv{669}(:,2),Te);
    Kv(670)=interp1(auxv{670}(:,1),auxv{670}(:,2),Te);
    Kv(671)=interp1(auxv{671}(:,1),auxv{671}(:,2),Te);
end 

%% Function without Te

if gas_heat==1
    % 基态存在电子得失
    K(7)=5.5e-16;
    K(8)=3.0e-16;
    K(12)=6.9e-16;
    K(15)=4.0e-18;
    K(16)=2.3e-16;
    K(29)=5.0e-19;
    K(39)=3.3e-16;
    K(48)=4.6e-40;
    K(68)=5.0e-16;
    K(71)=1.0e-38;
    K(73)=2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4));
    K(106)=1.6e-43;
    K(107)=1.6e-43;
    K(108)=2.2e-16;
    K(109)=1.0e-18;
    K(110)=2.5e-16;
    K(111)=4.0e-16;
    K(112)=1.0e-18;
    K(113)=1.0e-18;
    K(114)=1.0e-18;
    K(115)=1.0e-18;
    K(116)=5.0e-16;
    K(177)=5.0e-19;
    K(178)=8.3e-18;
    
    %基态无电子得失
    K(4)=9.0e-16;
    K(10)=6.0e-13;
    K(11)=1.7e-42*(Tg*1.16e4)^(-1.2);
    K(14)=3.31e-16;
    K(20)=1.0e-21;
    K(21)=3.0e-17;
    K(24)=2.6e-14*(300/(Tg*1.16e4))^0.44;
    K(25)=4.2e-13*(300/(Tg*1.16e4))^0.44;
    K(26)=2.01e-13*(300/(Tg*1.16e4))^0.5;
    K(27)=4.2e-13;
    K(28)=9.0e-41;
    K(30)=5.0e-13;
    K(31)=8.0e-17;
    K(35)=0.63*2.6e-16;
    K(36)=6.4e-17;
    K(38)=3.0e-13;
    K(40)=1.5e-40;
    K(41)=3.1e-40;
    K(42)=3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4));
    K(43)=1.0e-13;
    K(46)=1.0e-15;
    K(47)=5.2e-17;
    K(49)=5.5e-16;
    K(50)=8.0e-16;
    K(51)=4.0e-16;
    K(52)=2.5e-16;
    K(55)=9.4e-16;
    K(56)=4.5e-16;
    K(57)=9.62e-17;
    K(60)=1.1e-15;
    K(61)=5.2e-17;
    K(62)=4.7e-41;
    K(63)=5.0e-13;
    K(64)=3.0e-13;
    K(65)=1.1e-16;
    K(66)=1.4e-17;
    K(67)=1.4e-16;
    %K(70)=2.3e-41;
    K(70)=0;
    K(72)=1.0e-41;
    K(77)=3.0e-40;
    K(78)=1.1e-15;
    K(79)=9.0e-16;
    K(80)=2.6e-38;
    K(81)=4.2e-38;
    K(82)=5.0e-18;
    K(83)=1.0e-18;
    K(84)=5.0e-13;
    K(85)=5.0e-13;
    K(86)=6.0e-13;
    K(87)=5.0e-13;
    K(88)=5.0e-13;
    K(89)=6.0e-13;
    K(90)=1.0e-20;
    K(91)=5.0e-13;
    K(92)=5.0e-13;
    K(93)=6.0e-13;
    K(121)=1.2e-16;
    K(122)=7.0e-17;
    K(123)=4.9e-16;
    K(124)=7.0e-17;
    K(125)=9.0e-16;
    K(126)=1.9e-17;
    K(127)=7.0e-16;
    K(128)=6.6e-16;
    K(129)=1.0e-20;
    K(130)=2.9e-16;
    K(131)=2.9e-16;
    K(132)=1.0e-15;
    K(133)=1.0e-15;
    K(134)=2.0e-16;
    K(135)=2.0e-15;
    K(136)=1.0e-17;
    K(137)=5.0e-16;
    K(138)=1.0e-17;
    K(139)=2.6e-18;
    K(140)=7.0e-16;
    K(141)=2.8e-16;
    K(142)=5.0e-16;
    K(143)=9.0e-18;
    K(144)=1.0e-16;
    K(145)=5.0e-17;
    K(146)=9.0e-16;
    K(147)=1.8e-17;
    K(148)=4.0e-18;
    K(149)=3.0e-16;
    K(150)=2.0e-17;
    K(151)=3.0e-40*(Tg*11600)^(-1);
    K(152)=6.0e-13;
    K(153)=5.0e-13;
    K(154)=6.0e-13;
    K(155)=6.0e-13;
    K(156)=5.1e-13;
    K(157)=8.1e-13;
    K(158)=5.8e-13;
    K(159)=4.9e-13;
    K(160)=4.1e-13;
    K(161)=1.3e-13;
    K(162)=3.0e-13;
    K(163)=3.0e-13;
    K(164)=4.1e-13;
    K(165)=1.3e-13;
    K(166)=4.2e-13;
    K(167)=1.0e-13;
    K(168)=6.0e-13;
    K(169)=5.0e-13;
    K(170)=6.0e-13;
    K(171)=5.0e-13;
    K(172)=6.0e-13;
    K(173)=5.0e-13;
    K(174)=1.0e-25;
    K(175)=1.48e-16*exp(-16967/(Tg*1.16e4));
    K(176)=1.13e-37*(Tg*1.16e4)^(-2.16)*exp(-529/(Tg*1.16e4));
    
    %振动态
    if Vibration==1
        Kv(7)=K(46)*0.916;
        Kv(8)=K(46)*0.093;
        Kv(11)=K(56)*0.916;
        Kv(12)=K(56)*0.093;
        Kv(13)=K(60)*3.189;
        Kv(14)=K(60)*0.735;
        Kv(19)=7.14e-14*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
        Kv(20)=5.00e-14*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
        Kv(21)=5.00e-14*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
        Kv(22)=1.071e-15*exp(-137*(Tg*1.16e4)^(-1/3));
        Kv(23)=3.213e-15*exp(-137*(Tg*1.16e4)^(-1/3));
        Kv(24)=3.213e-15*exp(-137*(Tg*1.16e4)^(-1/3));
        Kv(25)=1.438e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
        Kv(26)=1.007e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
        Kv(27)=1.007e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
        Kv(28)=1.071e-15*exp(-137*(Tg*1.16e4)^(-1/3));
        Kv(29)=3.213e-15*exp(-137*(Tg*1.16e4)^(-1/3));
        Kv(30)=3.213e-15*exp(-137*(Tg*1.16e4)^(-1/3));
        Kv(59)=4.25e-7*exp(-407*(Tg*11600)^(-1/3)+824*(Tg*11600)^(-2/3));
        Kv(60)=Kv(59)*0.3;
        Kv(61)=Kv(59)*0.4;
        Kv(62)=8.57e-7*exp(-404*(Tg*11600)^(-1/3)+1096*(Tg*11600)^(-2/3));
        Kv(63)=Kv(62)*0.3;
        Kv(64)=Kv(62)*0.4;
        Kv(65)=1.43e-11*exp(-252*(Tg*11600)^(-1/3)+685*(Tg*11600)^(-2/3));
        Kv(66)=Kv(65)*0.3;
        Kv(67)=Kv(65)*0.4;
        %VT弛豫
        r(1,1)=10.555*Tg^(-1/2);%1→0,M=CO2
        r(2,1)=9.01*Tg^(-1/2);%1→0,M=CO
        r(3,1)=9.11*Tg^(-1/2);%1→0,M=O2
        r(1:3,2)=r(1:3,1);%2→1同1→0,均0.29
        r(1,3)=10.191*Tg^(-1/2);%3→2，0.28
        r(2,3)=8.7*Tg^(-1/2);
        r(3,3)=8.794*Tg^(-1/2);
        r(1:3,4)=r(1:3,3);
        r(1:3,5)=r(1:3,1);
        r(1,6)=9.827*Tg^(-1/2);
        r(2,6)=8.391*Tg^(-1/2);
        r(3,6)=8.48*Tg^(-1/2);
        r(1:3,7)=r(1:3,6);
        r(1:3,8)=r(1:3,6);
        Fr=0.5.*(3-exp(-2/3.*r)).*exp(-2/3.*r);
        Kv(68)=(2.0074*Kv(59)+2.002*Kv(62)+1.9697*Kv(65))*Fr(1,2)/Fr(1,1);
        Kv(69)=(2.0074*Kv(60)+2.002*Kv(63)+1.9697*Kv(66))*Fr(2,2)/Fr(2,1);
        Kv(70)=(2.0074*Kv(61)+2.002*Kv(64)+1.9697*Kv(67))*Fr(3,2)/Fr(3,1);
        Kv(71)=(3.0224*Kv(59)+3.006*Kv(62)+2.9106*Kv(65))*Fr(1,3)/Fr(1,1);
        Kv(72)=(3.0224*Kv(60)+3.006*Kv(63)+2.9106*Kv(66))*Fr(2,3)/Fr(2,1);
        Kv(73)=(3.0224*Kv(61)+3.006*Kv(64)+2.9106*Kv(67))*Fr(3,3)/Fr(3,1);
        Kv(74)=(4.0451*Kv(59)+4.012*Kv(62)+3.8238*Kv(65))*Fr(1,4)/Fr(1,1);
        Kv(75)=(4.0451*Kv(60)+4.012*Kv(63)+3.8238*Kv(66))*Fr(2,4)/Fr(2,1);
        Kv(76)=(4.0451*Kv(61)+4.012*Kv(64)+3.8238*Kv(67))*Fr(3,4)/Fr(3,1);
        Kv(77)=1.06e-11*exp(-242*(Tg*11600)^(-1/3)+633*(Tg*11600)^(-2/3));
        Kv(78:83)=2.13e-11*exp(-242*(Tg*11600)^(-1/3)+633*(Tg*11600)^(-2/3));
        Kv(84)=2.157e-15*exp(-88*(Tg*11600)^(-1/3)+233*(Tg*11600)^(-2/3));
        Kv(85)=5.305e-15*exp(-88.2*(Tg*11600)^(-1/3)+233*(Tg*11600)^(-2/3));
        Kv(86)=5.305e-15*exp(-88*(Tg*11600)^(-1/3)+233*(Tg*11600)^(-2/3));
        %V-N
        Kv(87)=3.91e-16*exp(-48688/(Tg*11600));
        Kv(88)=3.91e-16*exp(-47852/(Tg*11600));
        Kv(89)=3.91e-16*exp(-47110/(Tg*11600));
        Kv(90)=3.91e-16*exp(-46368/(Tg*11600));
        Kv(91)=3.91e-16*exp(-46739/(Tg*11600));
        Kv(92)=3.91e-16*exp(-44048/(Tg*11600));
        Kv(93)=3.91e-16*exp(-41449/(Tg*11600));
        Kv(94)=3.91e-16*exp(-38851/(Tg*11600));
        Kv(95)=2.8e-17*exp(-26036/(Tg*11600));
        Kv(96)=2.8e-17*exp(-25514/(Tg*11600));
        Kv(97)=2.8e-17*exp(-25050/(Tg*11600));
        Kv(98)=2.8e-17*exp(-24586/(Tg*11600));
        Kv(99)=2.8e-17*exp(-24818/(Tg*11600));
        Kv(100)=2.8e-17*exp(-23136/(Tg*11600));
        Kv(101)=2.8e-17*exp(-21512/(Tg*11600));
        Kv(102)=2.8e-17*exp(-19888/(Tg*11600));
        %离子反应补充
        Kv(127:136)=1.0e-21;
        Kv(137:146)=9.0e-41;
        Kv(147:156)=1.5e-40;
        Kv(157:166)=3.1e-40;
        Kv(167:174)=1.0e-15;
        Kv(175:184)=5.5e-16;
        Kv(185:194)=9.4e-16;
        Kv(195:202)=4.5e-16;
        Kv(203:210)=1.1e-15;
        Kv(211:220)=4.7e-41;
        %
        Kv(225:234)=1.7e-42*(Tg*1.16e4)^(-1.2);
        Kv(235:244)=4.0e-18;
        Kv(245:254)=9.0e-41;
        Kv(255:264)=3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4));
        %Kv(265:274)=2.3e-41;
        Kv(265:274)=0;
        Kv(275:284)=1.0e-38;
        Kv(285:304)=1.0e-41;
        Kv(305)=2.897e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
        Kv(306)=2.028e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
        Kv(307)=2.028e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
        Kv(308)=1.582e-11*exp(-272*(Tg*1.16e4)^(-1/3)+437*(Tg*1.16e4)^(-2/3));
        Kv(309)=4.745e-11*exp(-272*(Tg*1.16e4)^(-1/3)+437*(Tg*1.16e4)^(-2/3));
        Kv(310)=4.745e-11*exp(-272*(Tg*1.16e4)^(-1/3)+437*(Tg*1.16e4)^(-2/3));
        Kv(311)=4.321e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
        Kv(312)=3.025e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
        Kv(313)=3.025e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
        Kv(314)=1.775e-17*exp(-108*(Tg*1.16e4)^(-1/3)+165*(Tg*1.16e4)^(-2/3));
        Kv(315)=2.13e-17*exp(-108*(Tg*1.16e4)^(-1/3)+165*(Tg*1.16e4)^(-2/3));
        Kv(316)=2.13e-17*exp(-108*(Tg*1.16e4)^(-1/3)+165*(Tg*1.16e4)^(-2/3));
        Kv(317)=2.384e-15*exp(-89*(Tg*1.16e4)^(-1/3)+234*(Tg*1.16e4)^(-2/3));
        Kv(318)=6.48e-15*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3));
        Kv(319)=1.442e-14*exp(-88.2*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3));
        Kv(320)=2.628e-14*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3));
        Kv(321)=1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3));
        Kv(322)=2.927e-17*exp(21.4*(Tg*1.16e4)^(-1/3)-53*(Tg*1.16e4)^(-2/3));
        Kv(323)=4.825e-17*exp(20.1*(Tg*1.16e4)^(-1/3)-65*(Tg*1.16e4)^(-2/3));
        Kv(324)=1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3));
        Kv(325)=2.927e-17*exp(21.4*(Tg*1.16e4)^(-1/3)-53*(Tg*1.16e4)^(-2/3));
        Kv(326)=1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3));
        Kv(327:329)=1.0e-38;
        Kv(330:332)=1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3));
        Kv(357)=8.568e-7*exp(-404*(Tg*1.16e4)^(-1/3)+1096*(Tg*1.16e4)^(-2/3));
        Kv(358)=2.57e-7*exp(-404*(Tg*1.16e4)^(-1/3)+1096*(Tg*1.16e4)^(-2/3));
        Kv(359)=3.427e-7*exp(-404*(Tg*1.16e4)^(-1/3)+1096*(Tg*1.16e4)^(-2/3));
        Kv(360)=1.431e-11*exp(-252*(Tg*1.16e4)^(-1/3)+685*(Tg*1.16e4)^(-2/3));
        Kv(361)=4.293e-12*exp(-252*(Tg*1.16e4)^(-1/3)+685*(Tg*1.16e4)^(-2/3));
        Kv(362)=5.724e-12*exp(-252*(Tg*1.16e4)^(-1/3)+685*(Tg*1.16e4)^(-2/3));
        Kv(363)=7.143e-14*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
        Kv(364)=5e-14*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
        Kv(365)=5e-14*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
        Kv(366)=3.218e-15*exp(-137*(Tg*1.16e4)^(-1/3));
        Kv(367)=9.976e-15*exp(-137*(Tg*1.16e4)^(-1/3));
        Kv(368)=9.976e-15*exp(-137*(Tg*1.16e4)^(-1/3));
        Kv(369)=6.447e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
        Kv(370)=4.513e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
        Kv(371)=4.513e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
        Kv(372)=1.071e-15*exp(-137*(Tg*1.16e4)^(-1/3));
        Kv(373)=3.32e-15*exp(-137*(Tg*1.16e4)^(-1/3));
        Kv(374)=3.32e-15*exp(-137*(Tg*1.16e4)^(-1/3));
        Kv(375)=6.447e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
        Kv(376)=4.513e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
        Kv(377)=4.513e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
        Kv(378)=1.071e-15*exp(-137*(Tg*1.16e4)^(-1/3));
        Kv(379)=3.32e-15*exp(-137*(Tg*1.16e4)^(-1/3));
        Kv(380)=3.32e-15*exp(-137*(Tg*1.16e4)^(-1/3));
        Kv(381)=2.897e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
        Kv(382)=2.028e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
        Kv(383)=2.028e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
        Kv(384)=8.568e-7*exp(-406*(Tg*1.16e4)^(-1/3)+829*(Tg*1.16e4)^(-2/3));
        Kv(385)=2.57e-7*exp(-406*(Tg*1.16e4)^(-1/3)+829*(Tg*1.16e4)^(-2/3));
        Kv(386)=3.427e-7*exp(-406*(Tg*1.16e4)^(-1/3)+829*(Tg*1.16e4)^(-2/3));
        Kv(387)=1.725e-6*exp(-404*(Tg*1.16e4)^(-1/3)+1098*(Tg*1.16e4)^(-2/3));
        Kv(388)=5.175e-7*exp(-404*(Tg*1.16e4)^(-1/3)+1098*(Tg*1.16e4)^(-2/3));
        Kv(389)=6.9e-7*exp(-404*(Tg*1.16e4)^(-1/3)+1098*(Tg*1.16e4)^(-2/3));
        Kv(390)=2.882e-11*exp(-253*(Tg*1.16e4)^(-1/3)+683*(Tg*1.16e4)^(-2/3));
        Kv(391)=8.646e-12*exp(-253*(Tg*1.16e4)^(-1/3)+683*(Tg*1.16e4)^(-2/3));
        Kv(392)=1.153e-11*exp(-253*(Tg*1.16e4)^(-1/3)+683*(Tg*1.16e4)^(-2/3));
        Kv(393)=4.299e-11*exp(-241*(Tg*1.16e4)^(-1/3)+637*(Tg*1.16e4)^(-2/3));
        Kv(394)=4.299e-11*exp(-241*(Tg*1.16e4)^(-1/3)+635*(Tg*1.16e4)^(-2/3));
        Kv(395)=1.071e-15*exp(-87.7*(Tg*1.16e4)^(-1/3)+230*(Tg*1.16e4)^(-2/3));
        Kv(396)=2.157e-15*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3));
        Kv(397)=4.344e-15*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3));
        Kv(398)=6.48e-15*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3));
        Kv(399)=9.667e-15*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3));
        Kv(400)=9.667e-15*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3));
        Kv(401)=4.332e-14*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3));
        Kv(402)=5.848e-14*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3));
        Kv(403)=2.157e-15*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3));
        Kv(404)=2.378e-14*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3));
        Kv(405)=5.292e-14*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3));
        Kv(406)=1.066e-13*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3));
        Kv(407)=1.438e-13*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3));
        Kv(408)=5.305e-15*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3));
        Kv(409)=5.305e-15*exp(-88*(Tg*1.16e4)^(-1/3)+233*(Tg*1.16e4)^(-2/3));
        Kv(410:417)=Kv(87:94);
        Kv(418)=3.91e-16*exp(-46284/(Tg*11600));
        Kv(419)=3.91e-16*exp(-45514/(Tg*11600));
        Kv(420)=3.91e-16*exp(-44744/(Tg*11600));
        Kv(421:423)=Kv(418:420);
        Kv(424)=2.8e-17*exp(-24534/(Tg*11600));
        Kv(425)=2.8e-17*exp(-24052/(Tg*11600));
        Kv(426)=2.8e-17*exp(-23571/(Tg*11600));
        Kv(427:429)=1e-21;
        Kv(430:432)=9.0e-41;
        Kv(433:435)=1.5e-40;
        Kv(436:438)=3.1e-40;
        Kv(439:441)=9.0e-41;
        Kv(442:444)=1.0e-15;
        Kv(445:447)=5.5e-16;
        Kv(448:450)=9.4e-16;
        Kv(451:453)=4.5e-16;
        Kv(454:456)=1.1e-15;
        Kv(457:459)=0;
        Kv(460:462)=1.0e-41;
        Kv(463:465)=4.7e-41;
        Kv(466:468)=1.0e-41;
        Kv(469:471)=1.7e-42*(Tg*1.16e4)^(-1.2);
        Kv(472:474)=4.0e-18;
        Kv(475:477)=3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4));
        Kv(478:490)=2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4));
        Kv(523)=(5.075*Kv(59)+5.02*Kv(62)+4.711*Kv(65))*Fr(1,5)/Fr(1,1);
        Kv(524)=(5.075*Kv(60)+5.02*Kv(63)+4.711*Kv(66))*Fr(2,5)/Fr(2,1);
        Kv(525)=(5.075*Kv(61)+5.02*Kv(64)+4.711*Kv(67))*Fr(3,5)/Fr(3,1);
        Kv(526)=(6.114*Kv(59)+6.03*Kv(62)+5.572*Kv(65))*Fr(1,6)/Fr(1,1);
        Kv(527)=(6.114*Kv(60)+6.03*Kv(63)+5.572*Kv(66))*Fr(2,6)/Fr(2,1);
        Kv(528)=(6.114*Kv(61)+6.03*Kv(64)+5.572*Kv(67))*Fr(3,6)/Fr(3,1);
        Kv(529)=(7.16*Kv(59)+7.042*Kv(62)+6.409*Kv(65))*Fr(1,7)/Fr(1,1);
        Kv(530)=(7.16*Kv(60)+7.042*Kv(63)+6.409*Kv(66))*Fr(2,7)/Fr(2,1);
        Kv(531)=(7.16*Kv(61)+7.042*Kv(64)+6.409*Kv(67))*Fr(3,7)/Fr(3,1);
        Kv(532)=(8.214*Kv(59)+8.056*Kv(62)+7.223*Kv(65))*Fr(1,8)/Fr(1,1);
        Kv(533)=(8.214*Kv(60)+8.056*Kv(63)+7.223*Kv(66))*Fr(2,8)/Fr(2,1);
        Kv(534)=(8.214*Kv(61)+8.056*Kv(64)+7.223*Kv(67))*Fr(3,8)/Fr(3,1);
        Kv(535:542)=2.03e-11*exp(-242*(Tg*11600)^(-1/3)+633*(Tg*11600)^(-2/3));
        Kv(543)=7.199e-17*exp(18.4*(Tg*1.16e4)^(-1/3)-76*(Tg*1.16e4)^(-2/3));
        Kv(544)=4.825e-17*exp(20.1*(Tg*1.16e4)^(-1/3)-65*(Tg*1.16e4)^(-2/3));
        Kv(545)=7.199e-17*exp(18.4*(Tg*1.16e4)^(-1/3)-76*(Tg*1.16e4)^(-2/3));
        Kv(546)=2.927e-17*exp(21.4*(Tg*1.16e4)^(-1/3)-53*(Tg*1.16e4)^(-2/3));
        Kv(547)=4.825e-17*exp(20.1*(Tg*1.16e4)^(-1/3)-65*(Tg*1.16e4)^(-2/3));
        Kv(548)=7.199e-17*exp(18.4*(Tg*1.16e4)^(-1/3)-76*(Tg*1.16e4)^(-2/3));
        Kv(549)=1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3));
        Kv(550)=2.927e-17*exp(21.4*(Tg*1.16e4)^(-1/3)-53*(Tg*1.16e4)^(-2/3));
        Kv(551)=4.825e-17*exp(20.1*(Tg*1.16e4)^(-1/3)-65*(Tg*1.16e4)^(-2/3));
        Kv(552)=7.199e-17*exp(18.4*(Tg*1.16e4)^(-1/3)-76*(Tg*1.16e4)^(-2/3));
        Kv(553)=1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3));
        Kv(554)=2.927e-17*exp(21.4*(Tg*1.16e4)^(-1/3)-53*(Tg*1.16e4)^(-2/3));
        Kv(555)=4.825e-17*exp(20.1*(Tg*1.16e4)^(-1/3)-65*(Tg*1.16e4)^(-2/3));
        Kv(556)=1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3));
        Kv(557)=2.927e-17*exp(21.4*(Tg*1.16e4)^(-1/3)-53*(Tg*1.16e4)^(-2/3));
        Kv(558)=1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3));
        Kv(559:576)=1.453e-17*exp(22.1*(Tg*1.16e4)^(-1/3)-40.3*(Tg*1.16e4)^(-2/3));
        Kv(577)=3.91e-16*exp(-36160/(Tg*11600));
        Kv(578)=3.91e-16*exp(-33654/(Tg*11600));
        Kv(579)=3.91e-16*exp(-31148/(Tg*11600));
        Kv(580)=3.91e-16*exp(-28643/(Tg*11600));
        Kv(581:584)=Kv(577:580);
        Kv(585)=2.8e-17*exp(-18206/(Tg*11600));
        Kv(586)=2.8e-17*exp(-16640/(Tg*11600));
        Kv(587)=2.8e-17*exp(-15074/(Tg*11600));
        Kv(588)=2.8e-17*exp(-13508/(Tg*11600));
        Kv(589:592)=1e-21;
        Kv(593:596)=9.0e-41;
        Kv(597:600)=1.5e-40;
        Kv(601:604)=3.1e-40;
        Kv(605:608)=9.0e-41;
        Kv(609:612)=1.0e-38;
        Kv(613:616)=1.0e-15;
        Kv(617:620)=5.5e-16;
        Kv(621:624)=9.4e-16;
        Kv(625:628)=4.5e-16;
        Kv(629:632)=1.1e-15;
        Kv(633:636)=0;
        Kv(637:640)=1.0e-41;
        Kv(641:644)=4.7e-41;
        Kv(645:648)=1.0e-41;
        Kv(649:652)=1.7e-42*(Tg*1.16e4)^(-1.2);
        Kv(653:656)=4.0e-18;
        Kv(657:660)=3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4));
        Kv(661:664)=2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4));
        Kv(672)=1.43e-11*exp(-252*(Tg*11600)^(-1/3)+685*(Tg*11600)^(-2/3));
        Kv(673)=4.29e-12*exp(-252*(Tg*11600)^(-1/3)+685*(Tg*11600)^(-2/3));
        Kv(674)=5.72e-12*exp(-252*(Tg*11600)^(-1/3)+685*(Tg*11600)^(-2/3));
        Kv(675)=1.071e-15*exp(-137*(Tg*1.16e4)^(-1/3));
        Kv(676)=3.213e-15*exp(-137*(Tg*1.16e4)^(-1/3));
        Kv(677)=3.213e-15*exp(-137*(Tg*1.16e4)^(-1/3));
        Kv(678)=2.897e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
        Kv(679)=2.028e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
        Kv(680)=2.028e-13*exp(-177*(Tg*1.16e4)^(-1/3)+451*(Tg*1.16e4)^(-2/3));
        Kv(681:682)=3.91e-16*exp(-26230/(Tg*11600));
        Kv(683)=2.8e-17*exp(-12000/(Tg*11600));
        Kv(684)=1.0e-21;
        Kv(685)=9.0e-41;
        Kv(686)=1.5e-40;
        Kv(687)=3.1e-40;
        Kv(688)=9.0e-41;
        Kv(689)=1.0e-38;
        Kv(690)=1.0e-15;
        Kv(691)=5.5e-16;
        Kv(692)=9.4e-16;
        Kv(693)=4.5e-16;
        Kv(694)=1.1e-15;
        Kv(695)=0;
        Kv(696)=1.0e-41;
        Kv(697)=4.7e-41;
        Kv(698)=1.0e-41;
        Kv(699)=1.7e-42*(Tg*1.16e4)^(-1.2);
        Kv(700)=4.0e-18;
        Kv(701)=3.81e-42*(Tg*1.16e4)^(-1)*exp(-170/(Tg*1.16e4));
        Kv(702)=2.2e-41*(300/(Tg*1.16e4))^(1.5)*exp(-600/(Tg*1.16e4));
        Kv(703)=5.305e-15*exp(-88*(Tg*11600)^(-1/3)+233*(Tg*11600)^(-2/3));
        Kv(704:721)=3.0e-40;
        Kv(722:739)=2.6e-38;
        Kv(740:757)=4.2e-38;
        Kv(758:775)=1.0e-18;
        Kv(776:793)=1.0e-20;
        Kv(794)=K(125)*0.916;
        Kv(795)=K(125)*0.093;
        Kv(796:811)=9.0e-16;
        Kv(812:829)=1.0e-25;
        Kv(830:847)=1.13e-37*(Tg*1.16e4)^(-2.16)*exp(-529/(Tg*1.16e4));
        Kv(848:865)=8.3e-18;
        %K(1,number_of_reaction+1:number_of_reactions)=Kv(1,1:number_of_vreaction);
    end
end
    if sim==1
        Ralldlt=[5,6,11,12,13,17,21,22,23,24,29,32,33,34,36,37,41,42,44,45,46,47,48,49,50,51,52,58,59,60,61,62,66,67,68,69,70,71,73,74,76,78,80,82,83,84,85,86,89,93,94,97,98,100,101,102,106,107,110,112,113,114,115,116,117,119,120,121,122,123,124,126,127,128,129,130,131,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,164,165,166,167,168,169,170,171,172,173,175,177,178,183,184,185,186,188,191,192,196,200,201,202,204,205,206,207,208,225,226,227,228,229,230,231,232,233,234,235,237,238,239,240,241,242,244,245,247,248,250,251,253,254,265,266,267,268,269,270,271,272,273,274,275,276,277,278,279,280,283,284,287,288,292,293,294,295,296,297,298,299,300,302,303,304,305,306,307,308,309,310,311,312,313,314,319,320,325,326,327,328,329,330,331,332,333,334,335,336,337,338,339,340,341,342,343,344,345,346,347,348,349,350,351,352,353,354,355,356,357,358,359,360,361,362,368,376,381,382,383,384,385,386,387,388,389,390,391,392,393,394,395,396,397,398,403,404,405,406,407,408,409,410,411,412,416,417,418,427,428,433,434,435,436,437,438,439,440,441,442,443,444,445,446,447,448,449,450,451,452,453,454,455,456,457,458,459,460,461,462,466,467,468,476,477,478,484,485,486,487,488,489,490,491,492,493,494,495,496,497,498,505,506,507,515,516,518,519,520,521,522,523,524,525,526,527,528,530,531,532,533,534,535,536,537,538,539,540,542,543,544,545,546,547,548,549,550,551,552,553,554,555,556,557,558,559,560,561,562,563,564,565,566,567,568,569,570,572,575,576,577,578,579,580,581,582,583,584,585,586,587,588,589,590,591,592,593,594,595,596,597,598,599,600,601,602,603,604,605,606,607,608,609,610,611,612,613,614,615,616,617,618,619,620,621,622,623,624,625,627,628,630,631,632,633,634,635,636,637,638,639,640,641,642,643,644,645,646,647,648,649,650,651,652,653,654,655,656,657,658,659,660,661,662,663,664,665,666,667,668,681,682,683,684,685,686,687,688,701,702,703,704,705,706,707,708,709,710,711,712,755,756,757,758,759,760,761,762,763,764,765,766,767,768,769,770,775,776,777,778,779,780,781,782,787,788,789,790,791,792,793,794,795,796,797,798,807,808,809,810,811,812,813,814,819,820,821,822,827,828,829,830,835,836,837,838,839,840,841,842,843,844,845,846,847,848,849,850,851,852,853,854,855,856,857,858,859,860,861,862,863,864,865,866,867,868,869,870,871,872,873,874,875,876,877,878,879,880,881,882,883,884,885,886,887,888,889,890,891,892,893,894,895,896,897,898,899,900,901,902,903,904,905,906,907,908,909,910,911,912,913,914,915,916,917,918,919,921,922,923,924,925,926,927,928,929,930,931,932,933,934,935,936,937,938,939,940,941,942,943,944,945,946,947,948,949,950,951,952,953,954,955,956,957,958,959,960,961,962,963,964,965,966,967,968,969,970,971,977,983,984,989,990,991,992,993,994,995,996,997,998,999,1000,1001,1002,1003,1004,1005,1006,1007,1008,1009,1010,1011,1012,1013,1014,1015,1016,1017,1018,1019,1020,1021,1022,1023,1024,1025,1026,1027,1028,1029,1030,1031,1032,1033,1034,1035,1036,1037,1038,1039,1040,1041,1042,1043];
        K(intersect([1:length(K)],Ralldlt))=0;
        Kv(intersect([1:length(Kv)],Ralldlt-length(K)))=0;
    end
    K0=K;
    Kv0=Kv;

    %% dn
    % 基态存在电子得失
    R(1)=K(1)*Ngg(1,10);
    R(3)=K(3)*Ngg(1,10);
    R(6)=K(6)*Ngg(1,3);
    R(7)=K(7)*Ngg(6,2);
    R(8)=K(8)*Ngg(6,5);
    R(9)=K(9)*Ngg(1,8);
    R(12)=K(12)*Ngg(3,6);
    R(15)=K(15)*Ngg(6,10);
    R(16)=K(16)*Ngg(6,4);
    R(17)=K(17)*Ngg(1,2);
    R(19)=K(19)*Ngg(1,8);
    R(22)=K(22)*Ngg(1,3);
    R(23)=K(23)*Ngg(1,12);
    R(29)=K(29)*Ngg(2,13);
    R(32)=K(32)*Ngg(1,10);
    R(33)=K(33)*Ngg(1,5);
    R(34)=K(34)*Ngg(1,5);
    R(37)=K(37)*Ngg(1,5);
    R(39)=K(39)*Ngg(4,7);
    R(44)=K(44)*Ngg(1,10);
    R(45)=K(45)*Ngg(1,2);
    R(48)=K(48)*Ngg(1,5)*n(3);
    R(53)=K(53)*Ngg(1,10);
    R(54)=K(54)*Ngg(1,4);
    R(58)=K(58)*Ngg(1,2);
    R(59)=K(59)*Ngg(1,11);
    R(68)=K(68)*Ngg(6,11);
    R(69)=K(69)*Ngg(1,18);
    R(71)=K(71)*Ngg(1,12)*n(10);
    R(73)=K(73)*Ngg(1,3)*n(10);
    R(74)=K(74)*Ngg(1,19);
    R(75)=K(75)*Ngg(1,20);
    R(76)=K(76)*Ngg(1,21);
    R(95)=K(95)*Ngg(1,22);
    R(97)=K(97)*Ngg(1,23);
    R(98)=K(98)*Ngg(1,23);
    R(99)=K(99)*Ngg(1,24);
    R(100)=K(100)*Ngg(1,24);
    R(101)=K(101)*Ngg(1,25);
    R(102)=K(102)*Ngg(1,25);
    R(106)=K(106)*Ngg(1,25)*n(22);
    R(107)=K(107)*Ngg(1,25)*n(10);
    R(108)=K(108)*Ngg(6,26);
    R(109)=K(109)*Ngg(6,22);
    R(110)=K(110)*Ngg(6,23);
    R(111)=K(111)*Ngg(7,26);
    R(112)=K(112)*Ngg(4,33);
    R(113)=K(113)*Ngg(26,33);
    R(114)=K(114)*Ngg(4,34);
    R(115)=K(115)*Ngg(26,34);
    R(116)=K(116)*Ngg(22,32);
    R(117)=K(117)*Ngg(1,29);
    R(118)=K(118)*Ngg(1,28);
    R(119)=K(119)*Ngg(1,31);
    R(120)=K(120)*Ngg(1,30);
    R(177)=K(177)*Ngg(2,32);
    R(178)=K(178)*Ngg(10,32);
    
    % 基态无电子得失
    R(2)=K(2)*Ngg(1,10);
    R(4)=K(4)*Ngg(1,5);
    R(5)=K(5)*Ngg(1,3);
    R(10)=K(10)*Ngg(7,8);
    R(11)=K(11)*Ngg(4,3)*n(10);
    R(13)=K(13)*Ngg(1,10);
    R(14)=K(14)*Ngg(4,7);
    R(18)=K(18)*Ngg(1,2);
    R(20)=K(20)*Ngg(10,11);
    R(21)=K(21)*Ngg(11,3);
    R(24)=K(24)*Ngg(6,12);
    R(25)=K(25)*Ngg(6,12);
    R(26)=K(26)*Ngg(7,12);
    R(27)=K(27)*Ngg(7,12);
    R(28)=K(28)*Ngg(6,10)*n(10);
    R(30)=K(30)*Ngg(8,13);
    R(31)=K(31)*Ngg(4,13);
    R(35)=K(35)*Ngg(4,8);
    R(36)=K(36)*Ngg(3,8);
    R(38)=K(38)*Ngg(12,13);
    R(40)=K(40)*Ngg(2,6)*n(10);
    R(41)=K(41)*Ngg(3,6)*n(10);
    R(42)=K(42)*Ngg(4,4)*n(10);
    R(43)=K(43)*Ngg(6,8);
    R(46)=K(46)*Ngg(10,14);
    R(47)=K(47)*Ngg(11,12);
    R(49)=K(49)*Ngg(10,15);
    R(50)=K(50)*Ngg(5,6);
    R(51)=K(51)*Ngg(5,7);
    R(52)=K(52)*Ngg(4,15);
    R(55)=K(55)*Ngg(10,16);
    R(56)=K(56)*Ngg(10,16);
    R(57)=K(57)*Ngg(4,8);
    R(60)=K(60)*Ngg(10,17);
    R(61)=K(61)*Ngg(11,12);
    R(62)=K(62)*Ngg(7,3)*n(10);
    R(63)=K(63)*Ngg(8,9);
    R(64)=K(64)*Ngg(12,9);
    R(65)=K(65)*Ngg(4,9);
    R(66)=K(66)*Ngg(4,9);
    R(67)=K(67)*Ngg(4,9);
    R(70)=K(70)*Ngg(10,12)*n(10);
    R(72)=K(72)*Ngg(7,10)*n(10);
    R(77)=K(77)*Ngg(8,10)*n(2);
    R(78)=K(78)*Ngg(2,20);
    R(79)=K(79)*Ngg(2,21);
    R(80)=K(80)*Ngg(2,20)*n(10);
    R(81)=K(81)*Ngg(2,21)*n(10);
    R(82)=K(82)*Ngg(3,19);
    R(83)=K(83)*Ngg(10,19);
    R(84)=K(84)*Ngg(13,19);
    R(85)=K(85)*Ngg(9,19);
    R(86)=K(86)*Ngg(7,19);
    R(87)=K(87)*Ngg(13,20);
    R(88)=K(88)*Ngg(9,20);
    R(89)=K(89)*Ngg(7,20);
    R(90)=K(90)*Ngg(10,21);
    R(91)=K(91)*Ngg(13,21);
    R(92)=K(92)*Ngg(9,21);
    R(93)=K(93)*Ngg(7,21);
    R(94)=K(94)*Ngg(1,22);
    R(96)=K(96)*Ngg(1,22);
    R(103)=K(103)*Ngg(1,23);
    R(104)=K(104)*Ngg(1,25);
    R(105)=K(105)*Ngg(1,24);
    R(121)=K(121)*Ngg(8,23);
    R(122)=K(122)*Ngg(3,28);
    R(123)=K(123)*Ngg(23,28);
    R(124)=K(124)*Ngg(2,28);
    R(125)=K(125)*Ngg(10,28);
    R(126)=K(126)*Ngg(12,26);
    R(127)=K(127)*Ngg(12,23);
    R(128)=K(128)*Ngg(12,25);
    R(129)=K(129)*Ngg(5,29);
    R(130)=K(130)*Ngg(23,31);
    R(131)=K(131)*Ngg(23,30);
    R(132)=K(132)*Ngg(6,25);
    R(133)=K(133)*Ngg(6,25);
    R(134)=K(134)*Ngg(6,24);
    R(135)=K(135)*Ngg(7,25);
    R(136)=K(136)*Ngg(7,24);
    R(137)=K(137)*Ngg(7,27);
    R(138)=K(138)*Ngg(15,23);
    R(139)=K(139)*Ngg(15,23);
    R(140)=K(140)*Ngg(15,25);
    R(141)=K(141)*Ngg(15,25);
    R(142)=K(142)*Ngg(15,27);
    R(143)=K(143)*Ngg(13,23);
    R(144)=K(144)*Ngg(13,25);
    R(145)=K(145)*Ngg(9,23);
    R(146)=K(146)*Ngg(3,32);
    R(147)=K(147)*Ngg(5,33);
    R(148)=K(148)*Ngg(25,33);
    R(149)=K(149)*Ngg(27,33);
    R(150)=K(150)*Ngg(23,34);
    R(151)=K(151)*Ngg(6,23)*n(22);
    R(152)=K(152)*Ngg(8,33);
    R(153)=K(153)*Ngg(8,34);
    R(154)=K(154)*Ngg(13,29);
    R(155)=K(155)*Ngg(9,29);
    R(156)=K(156)*Ngg(29,33);
    R(157)=K(157)*Ngg(29,34);
    R(158)=K(158)*Ngg(7,29);
    R(159)=K(159)*Ngg(6,29);
    R(160)=K(160)*Ngg(12,33);
    R(161)=K(161)*Ngg(12,34);
    R(162)=K(162)*Ngg(13,28);
    R(163)=K(163)*Ngg(9,28);
    R(164)=K(164)*Ngg(28,33);
    R(165)=K(165)*Ngg(28,34);
    R(166)=K(166)*Ngg(7,28);
    R(167)=K(167)*Ngg(6,28);
    R(168)=K(168)*Ngg(21,33);
    R(169)=K(169)*Ngg(21,34);
    R(170)=K(170)*Ngg(20,33);
    R(171)=K(171)*Ngg(20,34);
    R(172)=K(172)*Ngg(19,33);
    R(173)=K(173)*Ngg(19,34);
    R(174)=K(174)*Ngg(10,26);
    R(175)=K(175)*Ngg(2,25);
    R(176)=K(176)*Ngg(4,23)*n(10);
    
    % 振动态
if Vibration==1
    Rv(1)=Kv(1)*Ngg(1,10);
    Rv(2)=Kv(2)*Ngg(1,10);
    Rv(3)=Kv(3)*Ngv(1,1);
    Rv(4)=Kv(4)*Ngv(1,2);
    Rv(5)=Kv(5)*Ngv(1,1);
    Rv(6)=Kv(6)*Ngv(1,2);
    Rv(7)=Kv(7)*Ngv(14,1);
    Rv(8)=Kv(8)*Ngv(14,2);
    Rv(9)=Kv(9)*Ngv(1,1);
    Rv(10)=Kv(10)*Ngv(1,2);
    Rv(11)=Kv(11)*Ngv(16,1);
    Rv(12)=Kv(12)*Ngv(16,2);
    Rv(13)=Kv(13)*Ngv(17,1);
    Rv(14)=Kv(14)*Ngv(17,2);
    Rv(15)=Kv(15)*Ngg(1,10);
    Rv(16)=Kv(16)*Ngg(1,10);
    Rv(17)=Kv(17)*Ngg(1,10);
    Rv(18)=Kv(18)*Ngg(1,10);
    Rv(19)=Kv(19)*Ngv(10,3);
    Rv(20)=Kv(20)*Ngv(2,3);
    Rv(21)=Kv(21)*Ngv(3,3);
    Rv(22)=Kv(22)*Ngv(10,4);
    Rv(23)=Kv(23)*Ngv(2,4);
    Rv(24)=Kv(24)*Ngv(3,4);
    Rv(25)=Kv(25)*Ngv(10,4);
    Rv(26)=Kv(26)*Ngv(2,4);
    Rv(27)=Kv(27)*Ngv(3,4);
    Rv(28)=Kv(28)*Ngv(10,5);
    Rv(29)=Kv(29)*Ngv(2,5);
    Rv(30)=Kv(30)*Ngv(3,5);
    Rv(31)=Kv(31)*Ngg(1,10);
    Rv(32)=Kv(32)*Ngg(1,10);
    Rv(33)=Kv(33)*Ngg(1,10);
    Rv(34)=Kv(34)*Ngg(1,10);
    Rv(35)=Kv(35)*Ngv(1,7);
    Rv(36)=Kv(36)*Ngv(1,8);
    Rv(37)=Kv(37)*Ngv(1,9);
    Rv(38)=Kv(38)*Ngv(1,10);
    Rv(39)=Kv(39)*Ngv(1,7);
    Rv(40)=Kv(40)*Ngv(1,8);
    Rv(41)=Kv(41)*Ngv(1,9);
    Rv(42)=Kv(42)*Ngv(1,10);
    Rv(43)=Kv(43)*Ngv(1,7);
    Rv(44)=Kv(44)*Ngv(1,8);
    Rv(45)=Kv(45)*Ngv(1,9);
    Rv(46)=Kv(46)*Ngv(1,10);
    Rv(47)=Kv(47)*Ngv(1,7);
    Rv(48)=Kv(48)*Ngv(1,8);
    Rv(49)=Kv(49)*Ngv(1,9);
    Rv(50)=Kv(50)*Ngv(1,10);
    Rv(51)=Kv(51)*Ngv(1,7);
    Rv(52)=Kv(52)*Ngv(1,8);
    Rv(53)=Kv(53)*Ngv(1,9);
    Rv(54)=Kv(54)*Ngv(1,10);
    Rv(55)=Kv(55)*Ngv(1,7);
    Rv(56)=Kv(56)*Ngv(1,8);
    Rv(57)=Kv(57)*Ngv(1,9);
    Rv(58)=Kv(58)*Ngv(1,10);
    Rv(59)=Kv(59)*Ngv(10,7);
    Rv(60)=Kv(60)*Ngv(2,7);
    Rv(61)=Kv(61)*Ngv(3,7);
    Rv(62)=Kv(62)*Ngv(10,7);
    Rv(63)=Kv(63)*Ngv(2,7);
    Rv(64)=Kv(64)*Ngv(3,7);
    Rv(65)=Kv(65)*Ngv(10,7);
    Rv(66)=Kv(66)*Ngv(2,7);
    Rv(67)=Kv(67)*Ngv(3,7);
    Rv(68)=Kv(68)*Ngv(10,8);
    Rv(69)=Kv(69)*Ngv(2,8);
    Rv(70)=Kv(70)*Ngv(3,8);
    Rv(71)=Kv(71)*Ngv(10,9);
    Rv(72)=Kv(72)*Ngv(2,9);
    Rv(73)=Kv(73)*Ngv(3,9);
    Rv(74)=Kv(74)*Ngv(10,10);
    Rv(75)=Kv(75)*Ngv(2,10);
    Rv(76)=Kv(76)*Ngv(3,10);
    Rv(77)=Kv(77)*Ngv(10,7);
    Rv(78)=Kv(78)*Ngv(10,8);
    Rv(79)=Kv(79)*Ngv(10,8);
    Rv(80)=Kv(80)*Ngv(10,9);
    Rv(81)=Kv(81)*Ngv(10,9);
    Rv(82)=Kv(82)*Ngv(10,10);
    Rv(83)=Kv(83)*Ngv(10,10);
    Rv(84)=Kv(84)*Ngv(10,4);
    Rv(85)=Kv(85)*Ngv(10,5);
    Rv(86)=Kv(86)*Nvv(3,5);
    Rv(87)=Kv(87)*Ngv(10,3);
    Rv(88)=Kv(88)*Ngv(10,4);
    Rv(89)=Kv(89)*Ngv(10,5);
    Rv(90)=Kv(90)*Ngv(10,6);
    Rv(91)=Kv(91)*Ngv(10,7);
    Rv(92)=Kv(92)*Ngv(10,8);
    Rv(93)=Kv(93)*Ngv(10,9);
    Rv(94)=Kv(94)*Ngv(10,10);
    Rv(95)=Kv(95)*Ngv(4,3);
    Rv(96)=Kv(96)*Ngv(4,4);
    Rv(97)=Kv(97)*Ngv(4,5);
    Rv(98)=Kv(98)*Ngv(4,6);
    Rv(99)=Kv(99)*Ngv(4,7);
    Rv(100)=Kv(100)*Ngv(4,8);
    Rv(101)=Kv(101)*Ngv(4,9);
    Rv(102)=Kv(102)*Ngv(4,10);
    Rv(103)=Kv(103)*Ngv(1,3);
    Rv(104)=Kv(104)*Ngv(1,4);
    Rv(105)=Kv(105)*Ngv(1,5);
    Rv(106)=Kv(106)*Ngv(1,6);
    Rv(107)=Kv(107)*Ngv(1,3);
    Rv(108)=Kv(108)*Ngv(1,4);
    Rv(109)=Kv(109)*Ngv(1,5);
    Rv(110)=Kv(110)*Ngv(1,6);
    Rv(111)=Kv(111)*Ngv(1,3);
    Rv(112)=Kv(112)*Ngv(1,4);
    Rv(113)=Kv(113)*Ngv(1,5);
    Rv(114)=Kv(114)*Ngv(1,6);
    Rv(115)=Kv(115)*Ngv(1,3);
    Rv(116)=Kv(116)*Ngv(1,4);
    Rv(117)=Kv(117)*Ngv(1,5);
    Rv(118)=Kv(118)*Ngv(1,6);
    Rv(119)=Kv(119)*Ngv(1,3);
    Rv(120)=Kv(120)*Ngv(1,4);
    Rv(121)=Kv(121)*Ngv(1,5);
    Rv(122)=Kv(122)*Ngv(1,6);
    Rv(123)=Kv(123)*Ngv(1,3);
    Rv(124)=Kv(124)*Ngv(1,4);
    Rv(125)=Kv(125)*Ngv(1,5);
    Rv(126)=Kv(126)*Ngv(1,6);
    Rv(127)=Kv(127)*Ngv(11,1);
    Rv(128)=Kv(128)*Ngv(11,2);
    Rv(129)=Kv(129)*Ngv(11,3);
    Rv(130)=Kv(130)*Ngv(11,4);
    Rv(131)=Kv(131)*Ngv(11,5);
    Rv(132)=Kv(132)*Ngv(11,6);
    Rv(133)=Kv(133)*Ngv(11,7);
    Rv(134)=Kv(134)*Ngv(11,8);
    Rv(135)=Kv(135)*Ngv(11,9);
    Rv(136)=Kv(136)*Ngv(11,10);
    Rv(137)=Kv(137)*Ngg(6,10)*nv(1);
    Rv(138)=Kv(138)*Ngg(6,10)*nv(2);
    Rv(139)=Kv(139)*Ngg(6,10)*nv(3);
    Rv(140)=Kv(140)*Ngg(6,10)*nv(4);
    Rv(141)=Kv(141)*Ngg(6,10)*nv(5);
    Rv(142)=Kv(142)*Ngg(6,10)*nv(6);
    Rv(143)=Kv(143)*Ngg(6,10)*nv(7);
    Rv(144)=Kv(144)*Ngg(6,10)*nv(8);
    Rv(145)=Kv(145)*Ngg(6,10)*nv(9);
    Rv(146)=Kv(146)*Ngg(6,10)*nv(10);
    Rv(147)=Kv(147)*Ngg(2,6)*nv(1);
    Rv(148)=Kv(148)*Ngg(2,6)*nv(2);
    Rv(149)=Kv(149)*Ngg(2,6)*nv(3);
    Rv(150)=Kv(150)*Ngg(2,6)*nv(4);
    Rv(151)=Kv(151)*Ngg(2,6)*nv(5);
    Rv(152)=Kv(152)*Ngg(2,6)*nv(6);
    Rv(153)=Kv(153)*Ngg(2,6)*nv(7);
    Rv(154)=Kv(154)*Ngg(2,6)*nv(8);
    Rv(155)=Kv(155)*Ngg(2,6)*nv(9);
    Rv(156)=Kv(156)*Ngg(2,6)*nv(10);
    Rv(157)=Kv(157)*Ngg(3,6)*nv(1);
    Rv(158)=Kv(158)*Ngg(3,6)*nv(2);
    Rv(159)=Kv(159)*Ngg(3,6)*nv(3);
    Rv(160)=Kv(160)*Ngg(3,6)*nv(4);
    Rv(161)=Kv(161)*Ngg(3,6)*nv(5);
    Rv(162)=Kv(162)*Ngg(3,6)*nv(6);
    Rv(163)=Kv(163)*Ngg(3,6)*nv(7);
    Rv(164)=Kv(164)*Ngg(3,6)*nv(8);
    Rv(165)=Kv(165)*Ngg(3,6)*nv(9);
    Rv(166)=Kv(166)*Ngg(3,6)*nv(10);
    Rv(167)=Kv(167)*Ngv(14,3);
    Rv(168)=Kv(168)*Ngv(14,4);
    Rv(169)=Kv(169)*Ngv(14,5);
    Rv(170)=Kv(170)*Ngv(14,6);
    Rv(171)=Kv(171)*Ngv(14,7);
    Rv(172)=Kv(172)*Ngv(14,8);
    Rv(173)=Kv(173)*Ngv(14,9);
    Rv(174)=Kv(174)*Ngv(14,10);
    Rv(175)=Kv(175)*Ngv(15,1);
    Rv(176)=Kv(176)*Ngv(15,2);
    Rv(177)=Kv(177)*Ngv(15,3);
    Rv(178)=Kv(178)*Ngv(15,4);
    Rv(179)=Kv(179)*Ngv(15,5);
    Rv(180)=Kv(180)*Ngv(15,6);
    Rv(181)=Kv(181)*Ngv(15,7);
    Rv(182)=Kv(182)*Ngv(15,8);
    Rv(183)=Kv(183)*Ngv(15,9);
    Rv(184)=Kv(184)*Ngv(15,10);
    Rv(185)=Kv(185)*Ngv(16,1);
    Rv(186)=Kv(186)*Ngv(16,2);
    Rv(187)=Kv(187)*Ngv(16,3);
    Rv(188)=Kv(188)*Ngv(16,4);
    Rv(189)=Kv(189)*Ngv(16,5);
    Rv(190)=Kv(190)*Ngv(16,6);
    Rv(191)=Kv(191)*Ngv(16,7);
    Rv(192)=Kv(192)*Ngv(16,8);
    Rv(193)=Kv(193)*Ngv(16,9);
    Rv(194)=Kv(194)*Ngv(16,10);
    Rv(195)=Kv(195)*Ngv(16,3);
    Rv(196)=Kv(196)*Ngv(16,4);
    Rv(197)=Kv(197)*Ngv(16,5);
    Rv(198)=Kv(198)*Ngv(16,6);
    Rv(199)=Kv(199)*Ngv(16,7);
    Rv(200)=Kv(200)*Ngv(16,8);
    Rv(201)=Kv(201)*Ngv(16,9);
    Rv(202)=Kv(202)*Ngv(16,10);
    Rv(203)=Kv(203)*Ngv(17,3);
    Rv(204)=Kv(204)*Ngv(17,4);
    Rv(205)=Kv(205)*Ngv(17,5);
    Rv(206)=Kv(206)*Ngv(17,6);
    Rv(207)=Kv(207)*Ngv(17,7);
    Rv(208)=Kv(208)*Ngv(17,8);
    Rv(209)=Kv(209)*Ngv(17,9);
    Rv(210)=Kv(210)*Ngv(17,10);
    Rv(211)=Kv(211)*Ngg(3,7)*nv(1);
    Rv(212)=Kv(212)*Ngg(3,7)*nv(2);
    Rv(213)=Kv(213)*Ngg(3,7)*nv(3);
    Rv(214)=Kv(214)*Ngg(3,7)*nv(4);
    Rv(215)=Kv(215)*Ngg(3,7)*nv(5);
    Rv(216)=Kv(216)*Ngg(3,7)*nv(6);
    Rv(217)=Kv(217)*Ngg(3,7)*nv(7);
    Rv(218)=Kv(218)*Ngg(3,7)*nv(8);
    Rv(219)=Kv(219)*Ngg(3,7)*nv(9);
    Rv(220)=Kv(220)*Ngg(3,7)*nv(10);
    Rv(221)=Kv(221)*Ngv(1,1);
    Rv(222)=Kv(222)*Ngv(1,2);
    Rv(223)=Kv(223)*Ngv(1,1);
    Rv(224)=Kv(224)*Ngv(1,2);
    Rv(225)=Kv(225)*Ngg(3,4)*nv(1);
    Rv(226)=Kv(226)*Ngg(3,4)*nv(2);
    Rv(227)=Kv(227)*Ngg(3,4)*nv(3);
    Rv(228)=Kv(228)*Ngg(3,4)*nv(4);
    Rv(229)=Kv(229)*Ngg(3,4)*nv(5);
    Rv(230)=Kv(230)*Ngg(3,4)*nv(6);
    Rv(231)=Kv(231)*Ngg(3,4)*nv(7);
    Rv(232)=Kv(232)*Ngg(3,4)*nv(8);
    Rv(233)=Kv(233)*Ngg(3,4)*nv(9);
    Rv(234)=Kv(234)*Ngg(3,4)*nv(10);
    Rv(235)=Kv(235)*Ngv(6,1);
    Rv(236)=Kv(236)*Ngv(6,2);
    Rv(237)=Kv(237)*Ngv(6,3);
    Rv(238)=Kv(238)*Ngv(6,4);
    Rv(239)=Kv(239)*Ngv(6,5);
    Rv(240)=Kv(240)*Ngv(6,6);
    Rv(241)=Kv(241)*Ngv(6,7);
    Rv(242)=Kv(242)*Ngv(6,8);
    Rv(243)=Kv(243)*Ngv(6,9);
    Rv(244)=Kv(244)*Ngv(6,10);
    Rv(245)=Kv(245)*Ngg(6,10)*nv(1);
    Rv(246)=Kv(246)*Ngg(6,10)*nv(2);
    Rv(247)=Kv(247)*Ngg(6,10)*nv(3);
    Rv(248)=Kv(248)*Ngg(6,10)*nv(4);
    Rv(249)=Kv(249)*Ngg(6,10)*nv(5);
    Rv(250)=Kv(250)*Ngg(6,10)*nv(6);
    Rv(251)=Kv(251)*Ngg(6,10)*nv(7);
    Rv(252)=Kv(252)*Ngg(6,10)*nv(8);
    Rv(253)=Kv(253)*Ngg(6,10)*nv(9);
    Rv(254)=Kv(254)*Ngg(6,10)*nv(10);
    Rv(255)=Kv(255)*Ngg(4,4)*nv(1);
    Rv(256)=Kv(256)*Ngg(4,4)*nv(2);
    Rv(257)=Kv(257)*Ngg(4,4)*nv(3);
    Rv(258)=Kv(258)*Ngg(4,4)*nv(4);
    Rv(259)=Kv(259)*Ngg(4,4)*nv(5);
    Rv(260)=Kv(260)*Ngg(4,4)*nv(6);
    Rv(261)=Kv(261)*Ngg(4,4)*nv(7);
    Rv(262)=Kv(262)*Ngg(4,4)*nv(8);
    Rv(263)=Kv(263)*Ngg(4,4)*nv(9);
    Rv(264)=Kv(264)*Ngg(4,4)*nv(10);
    Rv(265)=Kv(265)*Ngg(10,12)*nv(1);
    Rv(266)=Kv(266)*Ngg(10,12)*nv(2);
    Rv(267)=Kv(267)*Ngg(10,12)*nv(3);
    Rv(268)=Kv(268)*Ngg(10,12)*nv(4);
    Rv(269)=Kv(269)*Ngg(10,12)*nv(5);
    Rv(270)=Kv(270)*Ngg(10,12)*nv(6);
    Rv(271)=Kv(271)*Ngg(10,12)*nv(7);
    Rv(272)=Kv(272)*Ngg(10,12)*nv(8);
    Rv(273)=Kv(273)*Ngg(10,12)*nv(9);
    Rv(274)=Kv(274)*Ngg(10,12)*nv(10);
    Rv(275)=Kv(275)*Ngg(1,12)*nv(1);
    Rv(276)=Kv(276)*Ngg(1,12)*nv(2);
    Rv(277)=Kv(277)*Ngg(1,12)*nv(3);
    Rv(278)=Kv(278)*Ngg(1,12)*nv(4);
    Rv(279)=Kv(279)*Ngg(1,12)*nv(5);
    Rv(280)=Kv(280)*Ngg(1,12)*nv(6);
    Rv(281)=Kv(281)*Ngg(1,12)*nv(7);
    Rv(282)=Kv(282)*Ngg(1,12)*nv(8);
    Rv(283)=Kv(283)*Ngg(1,12)*nv(9);
    Rv(284)=Kv(284)*Ngg(1,12)*nv(10);
    Rv(285)=Kv(285)*Ngg(7,10)*nv(1);
    Rv(286)=Kv(286)*Ngg(7,10)*nv(2);
    Rv(287)=Kv(287)*Ngg(7,10)*nv(3);
    Rv(288)=Kv(288)*Ngg(7,10)*nv(4);
    Rv(289)=Kv(289)*Ngg(7,10)*nv(5);
    Rv(290)=Kv(290)*Ngg(7,10)*nv(6);
    Rv(291)=Kv(291)*Ngg(7,10)*nv(7);
    Rv(292)=Kv(292)*Ngg(7,10)*nv(8);
    Rv(293)=Kv(293)*Ngg(7,10)*nv(9);
    Rv(294)=Kv(294)*Ngg(7,10)*nv(10);
    Rv(295)=Kv(295)*Ngg(7,10)*nv(1);
    Rv(296)=Kv(296)*Ngg(7,10)*nv(2);
    Rv(297)=Kv(297)*Ngg(7,10)*nv(3);
    Rv(298)=Kv(298)*Ngg(7,10)*nv(4);
    Rv(299)=Kv(299)*Ngg(7,10)*nv(5);
    Rv(300)=Kv(300)*Ngg(7,10)*nv(6);
    Rv(301)=Kv(301)*Ngg(7,10)*nv(7);
    Rv(302)=Kv(302)*Ngg(7,10)*nv(8);
    Rv(303)=Kv(303)*Ngg(7,10)*nv(9);
    Rv(304)=Kv(304)*Ngg(7,10)*nv(10);
    Rv(305)=Kv(305)*Ngv(10,5);
    Rv(306)=Kv(306)*Ngv(2,5);
    Rv(307)=Kv(307)*Ngv(3,5);
    Rv(308)=Kv(308)*Ngv(10,6);
    Rv(309)=Kv(309)*Ngv(2,6);
    Rv(310)=Kv(310)*Ngv(3,6);
    Rv(311)=Kv(311)*Ngv(10,6);
    Rv(312)=Kv(312)*Ngv(2,6);
    Rv(313)=Kv(313)*Ngv(3,6);
    Rv(314)=Kv(314)*Ngv(10,6);
    Rv(315)=Kv(315)*Ngv(2,6);
    Rv(316)=Kv(316)*Ngv(3,6);
    Rv(317)=Kv(317)*Ngv(10,6);
    Rv(318)=Kv(318)*Ngv(10,6);
    Rv(319)=Kv(319)*Nvv(3,6);
    Rv(320)=Kv(320)*Nvv(4,6);
    Rv(321)=Kv(321)*Nvv(7,7);
    Rv(322)=Kv(322)*Nvv(7,8);
    Rv(323)=Kv(323)*Nvv(7,9);
    Rv(324)=Kv(324)*Nvv(8,8);
    Rv(325)=Kv(325)*Nvv(8,9);
    Rv(326)=Kv(326)*Nvv(9,9);
    Rv(327)=Kv(327)*Ngg(1,12)*nv(11);
    Rv(328)=Kv(328)*Ngg(1,12)*nv(12);
    Rv(329)=Kv(329)*Ngg(1,12)*nv(13);
    Rv(330)=Kv(330)*Nvv(7,10);
    Rv(331)=Kv(331)*Nvv(7,9);
    Rv(332)=Kv(332)*Nvv(8,10);
    Rv(333)=Kv(333)*Ngv(1,7);
    Rv(334)=Kv(334)*Ngv(1,8);
    Rv(335)=Kv(335)*Ngv(1,9);
    Rv(336)=Kv(336)*Ngg(1,10);
    Rv(337)=Kv(337)*Ngg(1,10);
    Rv(338)=Kv(338)*Ngg(1,10);
    Rv(339)=Kv(339)*Ngv(1,11);
    Rv(340)=Kv(340)*Ngv(1,12);
    Rv(341)=Kv(341)*Ngv(1,13);
    Rv(342)=Kv(342)*Ngv(1,11);
    Rv(343)=Kv(343)*Ngv(1,12);
    Rv(344)=Kv(344)*Ngv(1,13);
    Rv(345)=Kv(345)*Ngv(1,11);
    Rv(346)=Kv(346)*Ngv(1,12);
    Rv(347)=Kv(347)*Ngv(1,13);
    Rv(348)=Kv(348)*Ngv(1,11);
    Rv(349)=Kv(349)*Ngv(1,12);
    Rv(350)=Kv(350)*Ngv(1,13);
    Rv(351)=Kv(351)*Ngv(1,11);
    Rv(352)=Kv(352)*Ngv(1,12);
    Rv(353)=Kv(353)*Ngv(1,13);
    Rv(354)=Kv(354)*Ngv(1,11);
    Rv(355)=Kv(355)*Ngv(1,12);
    Rv(356)=Kv(356)*Ngv(1,13);
    Rv(357)=Kv(357)*Ngv(10,11);
    Rv(358)=Kv(358)*Ngv(2,11);
    Rv(359)=Kv(359)*Ngv(3,11);
    Rv(360)=Kv(360)*Ngv(10,11);
    Rv(361)=Kv(361)*Ngv(2,11);
    Rv(362)=Kv(362)*Ngv(3,11);
    Rv(363)=Kv(363)*Ngv(10,11);
    Rv(364)=Kv(364)*Ngv(2,11);
    Rv(365)=Kv(365)*Ngv(3,11);
    Rv(366)=Kv(366)*Ngv(10,12);
    Rv(367)=Kv(367)*Ngv(2,12);
    Rv(368)=Kv(368)*Ngv(3,12);
    Rv(369)=Kv(369)*Ngv(10,12);
    Rv(370)=Kv(370)*Ngv(2,12);
    Rv(371)=Kv(371)*Ngv(3,12);
    Rv(372)=Kv(372)*Ngv(10,12);
    Rv(373)=Kv(373)*Ngv(2,12);
    Rv(374)=Kv(374)*Ngv(3,12);
    Rv(375)=Kv(375)*Ngv(10,12);
    Rv(376)=Kv(376)*Ngv(2,12);
    Rv(377)=Kv(377)*Ngv(3,12);
    Rv(378)=Kv(378)*Ngv(10,13);
    Rv(379)=Kv(379)*Ngv(2,13);
    Rv(380)=Kv(380)*Ngv(3,13);
    Rv(381)=Kv(381)*Ngv(10,13);
    Rv(382)=Kv(382)*Ngv(2,13);
    Rv(383)=Kv(383)*Ngv(3,13);
    Rv(384)=Kv(384)*Ngv(10,8);
    Rv(385)=Kv(385)*Ngv(2,8);
    Rv(386)=Kv(386)*Ngv(3,8);
    Rv(387)=Kv(387)*Ngv(10,8);
    Rv(388)=Kv(388)*Ngv(2,8);
    Rv(389)=Kv(389)*Ngv(3,8);
    Rv(390)=Kv(390)*Ngv(10,8);
    Rv(391)=Kv(391)*Ngv(2,8);
    Rv(392)=Kv(392)*Ngv(3,8);
    Rv(393)=Kv(393)*Ngv(10,8);
    Rv(394)=Kv(394)*Ngv(10,8);
    Rv(395)=Kv(395)*Ngv(10,11);
    Rv(396)=Kv(396)*Nvv(3,11);
    Rv(397)=Kv(397)*Nvv(4,11);
    Rv(398)=Kv(398)*Nvv(5,11);
    Rv(399)=Kv(399)*Ngv(10,12);
    Rv(400)=Kv(400)*Nvv(3,12);
    Rv(401)=Kv(401)*Nvv(4,12);
    Rv(402)=Kv(402)*Nvv(5,12);
    Rv(403)=Kv(403)*Nvv(7,12);
    Rv(404)=Kv(404)*Ngv(10,13);
    Rv(405)=Kv(405)*Nvv(3,13);
    Rv(406)=Kv(406)*Nvv(4,13);
    Rv(407)=Kv(407)*Nvv(5,13);
    Rv(408)=Kv(408)*Nvv(7,13);
    Rv(409)=Kv(409)*Nvv(11,13);
    Rv(410:417)=Rv(87:94);
    Rv(418)=Kv(418)*Ngv(10,11);
    Rv(419)=Kv(419)*Ngv(10,12);
    Rv(420)=Kv(420)*Ngv(10,13);
    Rv(421:423)=Rv(418:420);
    Rv(424)=Kv(424)*Ngv(4,11);
    Rv(425)=Kv(425)*Ngv(4,12);
    Rv(426)=Kv(426)*Ngv(4,13);
    Rv(427)=Kv(427)*Ngv(11,11);
    Rv(428)=Kv(428)*Ngv(11,12);
    Rv(429)=Kv(429)*Ngv(11,13);
    Rv(430)=Kv(430)*Ngg(6,10)*nv(11);
    Rv(431)=Kv(431)*Ngg(6,10)*nv(12);
    Rv(432)=Kv(432)*Ngg(6,10)*nv(13);
    Rv(433)=Kv(433)*Ngg(2,6)*nv(11);
    Rv(434)=Kv(434)*Ngg(2,6)*nv(12);
    Rv(435)=Kv(435)*Ngg(2,6)*nv(13);
    Rv(436)=Kv(436)*Ngg(3,6)*nv(11);
    Rv(437)=Kv(437)*Ngg(3,6)*nv(12);
    Rv(438)=Kv(438)*Ngg(3,6)*nv(13);
    Rv(439:441)=Rv(430:432);
    Rv(442)=Kv(442)*Ngv(14,11);
    Rv(443)=Kv(443)*Ngv(14,12);
    Rv(444)=Kv(444)*Ngv(14,13);
    Rv(445)=Kv(445)*Ngv(15,11);
    Rv(446)=Kv(446)*Ngv(15,12);
    Rv(447)=Kv(447)*Ngv(15,13);
    Rv(448)=Kv(448)*Ngv(16,11);
    Rv(449)=Kv(449)*Ngv(16,12);
    Rv(450)=Kv(450)*Ngv(16,13);
    Rv(451)=Kv(451)*Ngv(16,11);
    Rv(452)=Kv(452)*Ngv(16,12);
    Rv(453)=Kv(453)*Ngv(16,13);
    Rv(454)=Kv(454)*Ngv(17,11);
    Rv(455)=Kv(455)*Ngv(17,12);
    Rv(456)=Kv(456)*Ngv(17,13);
    Rv(457)=Kv(457)*Ngg(10,12)*nv(11);
    Rv(458)=Kv(458)*Ngg(10,12)*nv(12);
    Rv(459)=Kv(459)*Ngg(10,12)*nv(13);
    Rv(460)=Kv(460)*Ngg(7,10)*nv(11);
    Rv(461)=Kv(461)*Ngg(7,10)*nv(12);
    Rv(462)=Kv(462)*Ngg(7,10)*nv(13);
    Rv(463)=Kv(463)*Ngg(7,3)*nv(11);
    Rv(464)=Kv(464)*Ngg(7,3)*nv(12);
    Rv(465)=Kv(465)*Ngg(7,3)*nv(13);
    Rv(466:468)=Rv(460:462);
    Rv(469)=Kv(469)*Ngg(3,4)*nv(11);
    Rv(470)=Kv(470)*Ngg(3,4)*nv(12);
    Rv(471)=Kv(471)*Ngg(3,4)*nv(13);
    Rv(472)=Kv(472)*Ngv(6,11);
    Rv(473)=Kv(473)*Ngv(6,12);
    Rv(474)=Kv(474)*Ngv(6,13);
    Rv(475)=Kv(475)*Ngg(4,4)*nv(11);
    Rv(476)=Kv(476)*Ngg(4,4)*nv(12);
    Rv(477)=Kv(477)*Ngg(4,4)*nv(13);
    Rv(478)=Kv(478)*Ngg(1,3)*nv(1);
    Rv(479)=Kv(479)*Ngg(1,3)*nv(2);
    Rv(480)=Kv(480)*Ngg(1,3)*nv(3);
    Rv(481)=Kv(481)*Ngg(1,3)*nv(4);
    Rv(482)=Kv(482)*Ngg(1,3)*nv(5);
    Rv(483)=Kv(483)*Ngg(1,3)*nv(6);
    Rv(484)=Kv(484)*Ngg(1,3)*nv(7);
    Rv(485)=Kv(485)*Ngg(1,3)*nv(8);
    Rv(486)=Kv(486)*Ngg(1,3)*nv(9);
    Rv(487)=Kv(487)*Ngg(1,3)*nv(10);
    Rv(488)=Kv(488)*Ngg(1,3)*nv(11);
    Rv(489)=Kv(489)*Ngg(1,3)*nv(12);
    Rv(490)=Kv(490)*Ngg(1,3)*nv(13);
    Rv(491)=Kv(491)*Ngg(1,10);
    Rv(492)=Kv(492)*Ngg(1,10);
    Rv(493)=Kv(493)*Ngg(1,10);
    Rv(494)=Kv(494)*Ngg(1,10);
    Rv(495)=Kv(495)*Ngv(1,10);
    Rv(496)=Kv(496)*Ngv(1,14);
    Rv(497)=Kv(497)*Ngv(1,15);
    Rv(498)=Kv(498)*Ngv(1,16);
    Rv(499)=Kv(499)*Ngv(1,14);
    Rv(500)=Kv(500)*Ngv(1,15);
    Rv(501)=Kv(501)*Ngv(1,16);
    Rv(502)=Kv(502)*Ngv(1,17);
    Rv(503)=Kv(503)*Ngv(1,14);
    Rv(504)=Kv(504)*Ngv(1,15);
    Rv(505)=Kv(505)*Ngv(1,16);
    Rv(506)=Kv(506)*Ngv(1,17);
    Rv(507)=Kv(507)*Ngv(1,14);
    Rv(508)=Kv(508)*Ngv(1,15);
    Rv(509)=Kv(509)*Ngv(1,16);
    Rv(510)=Kv(510)*Ngv(1,17);
    Rv(511)=Kv(511)*Ngv(1,14);
    Rv(512)=Kv(512)*Ngv(1,15);
    Rv(513)=Kv(513)*Ngv(1,16);
    Rv(514)=Kv(514)*Ngv(1,17);
    Rv(515)=Kv(515)*Ngv(1,14);
    Rv(516)=Kv(516)*Ngv(1,15);
    Rv(517)=Kv(517)*Ngv(1,16);
    Rv(518)=Kv(518)*Ngv(1,17);
    Rv(519)=Kv(519)*Ngv(1,14);
    Rv(520)=Kv(520)*Ngv(1,15);
    Rv(521)=Kv(521)*Ngv(1,16);
    Rv(522)=Kv(522)*Ngv(1,17);
    Rv(523)=Kv(523)*Ngv(10,14);
    Rv(524)=Kv(524)*Ngv(2,14);
    Rv(525)=Kv(525)*Ngv(3,14);
    Rv(526)=Kv(526)*Ngv(10,15);
    Rv(527)=Kv(527)*Ngv(2,15);
    Rv(528)=Kv(528)*Ngv(3,15);
    Rv(529)=Kv(529)*Ngv(10,16);
    Rv(530)=Kv(530)*Ngv(2,16);
    Rv(531)=Kv(531)*Ngv(3,16);
    Rv(532)=Kv(532)*Ngv(10,17);
    Rv(533)=Kv(533)*Ngv(2,17);
    Rv(534)=Kv(534)*Ngv(3,17);
    Rv(535)=Kv(535)*Ngv(10,14);
    Rv(536)=Kv(536)*Ngv(10,14);
    Rv(537)=Kv(537)*Ngv(10,15);
    Rv(538)=Kv(538)*Ngv(10,15);
    Rv(539)=Kv(539)*Ngv(10,16);
    Rv(540)=Kv(540)*Ngv(10,16);
    Rv(541)=Kv(541)*Ngv(10,17);
    Rv(542)=Kv(542)*Ngv(10,17);
    Rv(543)=Kv(543)*Nvv(7,10);
    Rv(544)=Kv(544)*Nvv(8,10);
    Rv(545)=Kv(545)*Nvv(8,14);
    Rv(546)=Kv(546)*Nvv(9,10);
    Rv(547)=Kv(547)*Nvv(9,14);
    Rv(548)=Kv(548)*Nvv(9,15);
    Rv(549)=Kv(549)*Nvv(10,10);
    Rv(550)=Kv(550)*Nvv(10,14);
    Rv(551)=Kv(551)*Nvv(10,15);
    Rv(552)=Kv(552)*Nvv(10,16);
    Rv(553)=Kv(553)*Nvv(14,14);
    Rv(554)=Kv(554)*Nvv(14,15);
    Rv(555)=Kv(555)*Nvv(14,16);
    Rv(556)=Kv(556)*Nvv(15,15);
    Rv(557)=Kv(557)*Nvv(15,16);
    Rv(558)=Kv(558)*Nvv(16,16);
    Rv(559)=Kv(559)*Nvv(7,14);
    Rv(560)=Kv(560)*Nvv(7,15);
    Rv(561)=Kv(561)*Nvv(7,16);
    Rv(562)=Kv(562)*Nvv(7,17);
    Rv(563)=Kv(563)*Nvv(8,14);
    Rv(564)=Kv(564)*Nvv(8,15);
    Rv(565)=Kv(565)*Nvv(8,16);
    Rv(566)=Kv(566)*Nvv(8,17);
    Rv(567)=Kv(567)*Nvv(9,14);
    Rv(568)=Kv(568)*Nvv(9,15);
    Rv(569)=Kv(569)*Nvv(9,16);
    Rv(570)=Kv(570)*Nvv(9,17);
    Rv(571)=Kv(571)*Nvv(10,15);
    Rv(572)=Kv(572)*Nvv(10,16);
    Rv(573)=Kv(573)*Nvv(10,17);
    Rv(574)=Kv(574)*Nvv(14,16);
    Rv(575)=Kv(575)*Nvv(14,17);
    Rv(576)=Kv(576)*Nvv(15,17);
    Rv(577)=Kv(577)*Ngv(10,14);
    Rv(578)=Kv(578)*Ngv(10,15);
    Rv(579)=Kv(579)*Ngv(10,16);
    Rv(580)=Kv(580)*Ngv(10,17);
    Rv(581)=Kv(581)*Ngv(10,14);
    Rv(582)=Kv(582)*Ngv(10,15);
    Rv(583)=Kv(583)*Ngv(10,16);
    Rv(584)=Kv(584)*Ngv(10,17);
    Rv(585)=Kv(585)*Ngv(4,14);
    Rv(586)=Kv(586)*Ngv(4,15);
    Rv(587)=Kv(587)*Ngv(4,16);
    Rv(588)=Kv(588)*Ngv(4,17);
    Rv(589)=Kv(589)*Ngv(11,14);
    Rv(590)=Kv(590)*Ngv(11,15);
    Rv(591)=Kv(591)*Ngv(11,16);
    Rv(592)=Kv(592)*Ngv(11,17);
    Rv(593)=Kv(593)*Ngg(6,10)*nv(14);
    Rv(594)=Kv(594)*Ngg(6,10)*nv(15);
    Rv(595)=Kv(595)*Ngg(6,10)*nv(16);
    Rv(596)=Kv(596)*Ngg(6,10)*nv(17);
    Rv(597)=Kv(597)*Ngg(2,6)*nv(14);
    Rv(598)=Kv(598)*Ngg(2,6)*nv(15);
    Rv(599)=Kv(599)*Ngg(2,6)*nv(16);
    Rv(600)=Kv(600)*Ngg(2,6)*nv(17);
    Rv(601)=Kv(601)*Ngg(3,6)*nv(14);
    Rv(602)=Kv(602)*Ngg(3,6)*nv(15);
    Rv(603)=Kv(603)*Ngg(3,6)*nv(16);
    Rv(604)=Kv(604)*Ngg(3,6)*nv(17);
    Rv(605)=Kv(605)*Ngg(6,10)*nv(14);
    Rv(606)=Kv(606)*Ngg(6,10)*nv(15);
    Rv(607)=Kv(607)*Ngg(6,10)*nv(16);
    Rv(608)=Kv(608)*Ngg(6,10)*nv(17);
    Rv(609)=Kv(609)*Ngg(1,12)*nv(14);
    Rv(610)=Kv(610)*Ngg(1,12)*nv(15);
    Rv(611)=Kv(611)*Ngg(1,12)*nv(16);
    Rv(612)=Kv(612)*Ngg(1,12)*nv(17);
    Rv(613)=Kv(613)*Ngv(14,14);
    Rv(614)=Kv(614)*Ngv(14,15);
    Rv(615)=Kv(615)*Ngv(14,16);
    Rv(616)=Kv(616)*Ngv(14,17);
    Rv(617)=Kv(617)*Ngv(15,14);
    Rv(618)=Kv(618)*Ngv(15,15);
    Rv(619)=Kv(619)*Ngv(15,16);
    Rv(620)=Kv(620)*Ngv(15,17);
    Rv(621)=Kv(621)*Ngv(16,14);
    Rv(622)=Kv(622)*Ngv(16,15);
    Rv(623)=Kv(623)*Ngv(16,16);
    Rv(624)=Kv(624)*Ngv(16,17);
    Rv(625)=Kv(625)*Ngv(16,14);
    Rv(626)=Kv(626)*Ngv(16,15);
    Rv(627)=Kv(627)*Ngv(16,16);
    Rv(628)=Kv(628)*Ngv(16,17);
    Rv(629)=Kv(629)*Ngv(17,14);
    Rv(630)=Kv(630)*Ngv(17,15);
    Rv(631)=Kv(631)*Ngv(17,16);
    Rv(632)=Kv(632)*Ngv(17,17);
    Rv(633)=Kv(633)*Ngg(10,12)*nv(14);
    Rv(634)=Kv(634)*Ngg(10,12)*nv(15);
    Rv(635)=Kv(635)*Ngg(10,12)*nv(16);
    Rv(636)=Kv(636)*Ngg(10,12)*nv(17);
    Rv(637)=Kv(637)*Ngg(7,10)*nv(14);
    Rv(638)=Kv(638)*Ngg(7,10)*nv(15);
    Rv(639)=Kv(639)*Ngg(7,10)*nv(16);
    Rv(640)=Kv(640)*Ngg(7,10)*nv(17);
    Rv(641)=Kv(641)*Ngg(3,7)*nv(14);
    Rv(642)=Kv(642)*Ngg(3,7)*nv(15);
    Rv(643)=Kv(643)*Ngg(3,7)*nv(16);
    Rv(644)=Kv(644)*Ngg(3,7)*nv(17);
    Rv(645)=Kv(645)*Ngg(7,10)*nv(14);
    Rv(646)=Kv(646)*Ngg(7,10)*nv(15);
    Rv(647)=Kv(647)*Ngg(7,10)*nv(16);
    Rv(648)=Kv(648)*Ngg(7,10)*nv(17);
    Rv(649)=Kv(649)*Ngg(3,4)*nv(14);
    Rv(650)=Kv(650)*Ngg(3,4)*nv(15);
    Rv(651)=Kv(651)*Ngg(3,4)*nv(16);
    Rv(652)=Kv(652)*Ngg(3,4)*nv(17);
    Rv(653)=Kv(653)*Ngv(6,14);
    Rv(654)=Kv(654)*Ngv(6,15);
    Rv(655)=Kv(655)*Ngv(6,16);
    Rv(656)=Kv(656)*Ngv(6,17);
    Rv(657)=Kv(657)*Ngg(4,4)*nv(14);
    Rv(658)=Kv(658)*Ngg(4,4)*nv(15);
    Rv(659)=Kv(659)*Ngg(4,4)*nv(16);
    Rv(660)=Kv(660)*Ngg(4,4)*nv(17);
    Rv(661)=Kv(661)*Ngg(1,3)*nv(14);
    Rv(662)=Kv(662)*Ngg(1,3)*nv(15);
    Rv(663)=Kv(663)*Ngg(1,3)*nv(16);
    Rv(664)=Kv(664)*Ngg(1,3)*nv(17);
    Rv(665)=Kv(665)*Ngg(1,10);
    Rv(666)=Kv(666)*Ngv(1,18);
    Rv(667)=Kv(667)*Ngv(1,18);
    Rv(668)=Kv(668)*Ngv(1,18);
    Rv(669)=Kv(669)*Ngv(1,18);
    Rv(670)=Kv(670)*Ngv(1,18);
    Rv(671)=Kv(671)*Ngv(1,18);
    Rv(672)=Kv(672)*Ngv(10,7);
    Rv(673)=Kv(673)*Ngv(2,7);
    Rv(674)=Kv(674)*Ngv(3,7);
    Rv(675)=Kv(675)*Ngv(10,18);
    Rv(676)=Kv(676)*Ngv(2,18);
    Rv(677)=Kv(677)*Ngv(3,18);
    Rv(678)=Kv(678)*Ngv(10,18);
    Rv(679)=Kv(679)*Ngv(2,18);
    Rv(680)=Kv(680)*Ngv(3,18);
    Rv(681)=Kv(681)*Ngv(10,18);
    Rv(682)=Kv(682)*Ngv(10,18);
    Rv(683)=Kv(683)*Ngv(4,18);
    Rv(684)=Kv(684)*Ngv(11,18);
    Rv(685)=Kv(685)*Ngg(6,10)*nv(18);
    Rv(686)=Kv(686)*Ngg(2,6)*nv(18);
    Rv(687)=Kv(687)*Ngg(3,6)*nv(18);
    Rv(688)=Kv(688)*Ngg(6,10)*nv(18);
    Rv(689)=Kv(689)*Ngg(1,12)*nv(18);
    Rv(690)=Kv(690)*Ngv(14,18);
    Rv(691)=Kv(691)*Ngv(15,18);
    Rv(692)=Kv(692)*Ngv(16,18);
    Rv(693)=Kv(693)*Ngv(16,18);
    Rv(694)=Kv(694)*Ngv(17,18);
    Rv(695)=Kv(695)*Ngg(10,12)*nv(18);
    Rv(696)=Kv(696)*Ngg(7,10)*nv(18);
    Rv(697)=Kv(697)*Ngg(3,7)*nv(18);
    Rv(698)=Kv(698)*Ngg(7,10)*nv(18);
    Rv(699)=Kv(699)*Ngg(3,4)*nv(18);
    Rv(700)=Kv(700)*Ngv(6,18);
    Rv(701)=Kv(701)*Ngg(4,4)*nv(18);
    Rv(702)=Kv(702)*Ngg(1,3)*nv(18);
    Rv(703)=Kv(703)*Ngv(10,18);
    Rv(704)=Kv(704)*Ngv(8,1)*n(2);
    Rv(705)=Kv(705)*Ngv(8,2)*n(2);
    Rv(706)=Kv(706)*Ngv(8,3)*n(2);
    Rv(707)=Kv(707)*Ngv(8,4)*n(2);
    Rv(708)=Kv(708)*Ngv(8,5)*n(2);
    Rv(709)=Kv(709)*Ngv(8,6)*n(2);
    Rv(710)=Kv(710)*Ngv(8,7)*n(2);
    Rv(711)=Kv(711)*Ngv(8,8)*n(2);
    Rv(712)=Kv(712)*Ngv(8,9)*n(2);
    Rv(713)=Kv(713)*Ngv(8,10)*n(2);
    Rv(714)=Kv(714)*Ngv(8,11)*n(2);
    Rv(715)=Kv(715)*Ngv(8,12)*n(2);
    Rv(716)=Kv(716)*Ngv(8,13)*n(2);
    Rv(717)=Kv(717)*Ngv(8,14)*n(2);
    Rv(718)=Kv(718)*Ngv(8,15)*n(2);
    Rv(719)=Kv(719)*Ngv(8,16)*n(2);
    Rv(720)=Kv(720)*Ngv(8,17)*n(2);
    Rv(721)=Kv(721)*Ngv(8,18)*n(2);
    Rv(722)=Kv(722)*Ngv(20,1)*n(2);
    Rv(723)=Kv(723)*Ngv(20,2)*n(2);
    Rv(724)=Kv(724)*Ngv(20,3)*n(2);
    Rv(725)=Kv(725)*Ngv(20,4)*n(2);
    Rv(726)=Kv(726)*Ngv(20,5)*n(2);
    Rv(727)=Kv(727)*Ngv(20,6)*n(2);
    Rv(728)=Kv(728)*Ngv(20,7)*n(2);
    Rv(729)=Kv(729)*Ngv(20,8)*n(2);
    Rv(730)=Kv(730)*Ngv(20,9)*n(2);
    Rv(731)=Kv(731)*Ngv(20,10)*n(2);
    Rv(732)=Kv(732)*Ngv(20,11)*n(2);
    Rv(733)=Kv(733)*Ngv(20,12)*n(2);
    Rv(734)=Kv(734)*Ngv(20,13)*n(2);
    Rv(735)=Kv(735)*Ngv(20,14)*n(2);
    Rv(736)=Kv(736)*Ngv(20,15)*n(2);
    Rv(737)=Kv(737)*Ngv(20,16)*n(2);
    Rv(738)=Kv(738)*Ngv(20,17)*n(2);
    Rv(739)=Kv(739)*Ngv(20,18)*n(2);
    Rv(740)=Kv(740)*Ngv(21,1)*n(2);
    Rv(741)=Kv(741)*Ngv(21,2)*n(2);
    Rv(742)=Kv(742)*Ngv(21,3)*n(2);
    Rv(743)=Kv(743)*Ngv(21,4)*n(2);
    Rv(744)=Kv(744)*Ngv(21,5)*n(2);
    Rv(745)=Kv(745)*Ngv(21,6)*n(2);
    Rv(746)=Kv(746)*Ngv(21,7)*n(2);
    Rv(747)=Kv(747)*Ngv(21,8)*n(2);
    Rv(748)=Kv(748)*Ngv(21,9)*n(2);
    Rv(749)=Kv(749)*Ngv(21,10)*n(2);
    Rv(750)=Kv(750)*Ngv(21,11)*n(2);
    Rv(751)=Kv(751)*Ngv(21,12)*n(2);
    Rv(752)=Kv(752)*Ngv(21,13)*n(2);
    Rv(753)=Kv(753)*Ngv(21,14)*n(2);
    Rv(754)=Kv(754)*Ngv(21,15)*n(2);
    Rv(755)=Kv(755)*Ngv(21,16)*n(2);
    Rv(756)=Kv(756)*Ngv(21,17)*n(2);
    Rv(757)=Kv(757)*Ngv(21,18)*n(2);
    Rv(758)=Kv(758)*Ngv(19,1);
    Rv(759)=Kv(759)*Ngv(19,2);
    Rv(760)=Kv(760)*Ngv(19,3);
    Rv(761)=Kv(761)*Ngv(19,4);
    Rv(762)=Kv(762)*Ngv(19,5);
    Rv(763)=Kv(763)*Ngv(19,6);
    Rv(764)=Kv(764)*Ngv(19,7);
    Rv(765)=Kv(765)*Ngv(19,8);
    Rv(766)=Kv(766)*Ngv(19,9);
    Rv(767)=Kv(767)*Ngv(19,10);
    Rv(768)=Kv(768)*Ngv(19,11);
    Rv(769)=Kv(769)*Ngv(19,12);
    Rv(770)=Kv(770)*Ngv(19,13);
    Rv(771)=Kv(771)*Ngv(19,14);
    Rv(772)=Kv(772)*Ngv(19,15);
    Rv(773)=Kv(773)*Ngv(19,16);
    Rv(774)=Kv(774)*Ngv(19,17);
    Rv(775)=Kv(775)*Ngv(19,18);
    Rv(776)=Kv(776)*Ngv(21,1);
    Rv(777)=Kv(777)*Ngv(21,2);
    Rv(778)=Kv(778)*Ngv(21,3);
    Rv(779)=Kv(779)*Ngv(21,4);
    Rv(780)=Kv(780)*Ngv(21,5);
    Rv(781)=Kv(781)*Ngv(21,6);
    Rv(782)=Kv(782)*Ngv(21,7);
    Rv(783)=Kv(783)*Ngv(21,8);
    Rv(784)=Kv(784)*Ngv(21,9);
    Rv(785)=Kv(785)*Ngv(21,10);
    Rv(786)=Kv(786)*Ngv(21,11);
    Rv(787)=Kv(787)*Ngv(21,12);
    Rv(788)=Kv(788)*Ngv(21,13);
    Rv(789)=Kv(789)*Ngv(21,14);
    Rv(790)=Kv(790)*Ngv(21,15);
    Rv(791)=Kv(791)*Ngv(21,16);
    Rv(792)=Kv(792)*Ngv(21,17);
    Rv(793)=Kv(793)*Ngv(21,18);
    Rv(794)=Kv(794)*Ngv(28,1);
    Rv(795)=Kv(795)*Ngv(28,2);
    Rv(796)=Kv(796)*Ngv(28,3);
    Rv(797)=Kv(797)*Ngv(28,4);
    Rv(798)=Kv(798)*Ngv(28,5);
    Rv(799)=Kv(799)*Ngv(28,6);
    Rv(800)=Kv(800)*Ngv(28,7);
    Rv(801)=Kv(801)*Ngv(28,8);
    Rv(802)=Kv(802)*Ngv(28,9);
    Rv(803)=Kv(803)*Ngv(28,10);
    Rv(804)=Kv(804)*Ngv(28,11);
    Rv(805)=Kv(805)*Ngv(28,12);
    Rv(806)=Kv(806)*Ngv(28,13);
    Rv(807)=Kv(807)*Ngv(28,14);
    Rv(808)=Kv(808)*Ngv(28,15);
    Rv(809)=Kv(809)*Ngv(28,16);
    Rv(810)=Kv(810)*Ngv(28,17);
    Rv(811)=Kv(811)*Ngv(28,18);
    Rv(812)=Kv(812)*Ngv(26,1);
    Rv(813)=Kv(813)*Ngv(26,2);
    Rv(814)=Kv(814)*Ngv(26,3);
    Rv(815)=Kv(815)*Ngv(26,4);
    Rv(816)=Kv(816)*Ngv(26,5);
    Rv(817)=Kv(817)*Ngv(26,6);
    Rv(818)=Kv(818)*Ngv(26,7);
    Rv(819)=Kv(819)*Ngv(26,8);
    Rv(820)=Kv(820)*Ngv(26,9);
    Rv(821)=Kv(821)*Ngv(26,10);
    Rv(822)=Kv(822)*Ngv(26,11);
    Rv(823)=Kv(823)*Ngv(26,12);
    Rv(824)=Kv(824)*Ngv(26,13);
    Rv(825)=Kv(825)*Ngv(26,14);
    Rv(826)=Kv(826)*Ngv(26,15);
    Rv(827)=Kv(827)*Ngv(26,16);
    Rv(828)=Kv(828)*Ngv(26,17);
    Rv(829)=Kv(829)*Ngv(26,18);
    Rv(830)=Kv(830)*Ngv(23,1)*n(4);
    Rv(831)=Kv(831)*Ngv(23,2)*n(4);
    Rv(832)=Kv(832)*Ngv(23,3)*n(4);
    Rv(833)=Kv(833)*Ngv(23,4)*n(4);
    Rv(834)=Kv(834)*Ngv(23,5)*n(4);
    Rv(835)=Kv(835)*Ngv(23,6)*n(4);
    Rv(836)=Kv(836)*Ngv(23,7)*n(4);
    Rv(837)=Kv(837)*Ngv(23,8)*n(4);
    Rv(838)=Kv(838)*Ngv(23,9)*n(4);
    Rv(839)=Kv(839)*Ngv(23,10)*n(4);
    Rv(840)=Kv(840)*Ngv(23,11)*n(4);
    Rv(841)=Kv(841)*Ngv(23,12)*n(4);
    Rv(842)=Kv(842)*Ngv(23,13)*n(4);
    Rv(843)=Kv(843)*Ngv(23,14)*n(4);
    Rv(844)=Kv(844)*Ngv(23,15)*n(4);
    Rv(845)=Kv(845)*Ngv(23,16)*n(4);
    Rv(846)=Kv(846)*Ngv(23,17)*n(4);
    Rv(847)=Kv(847)*Ngv(23,18)*n(4);
    Rv(848)=Kv(848)*Ngv(32,1);
    Rv(849)=Kv(849)*Ngv(32,2);
    Rv(850)=Kv(850)*Ngv(32,3);
    Rv(851)=Kv(851)*Ngv(32,4);
    Rv(852)=Kv(852)*Ngv(32,5);
    Rv(853)=Kv(853)*Ngv(32,6);
    Rv(854)=Kv(854)*Ngv(32,7);
    Rv(855)=Kv(855)*Ngv(32,8);
    Rv(856)=Kv(856)*Ngv(32,9);
    Rv(857)=Kv(857)*Ngv(32,10);
    Rv(858)=Kv(858)*Ngv(32,11);
    Rv(859)=Kv(859)*Ngv(32,12);
    Rv(860)=Kv(860)*Ngv(32,13);
    Rv(861)=Kv(861)*Ngv(32,14);
    Rv(862)=Kv(862)*Ngv(32,15);
    Rv(863)=Kv(863)*Ngv(32,16);
    Rv(864)=Kv(864)*Ngv(32,17);
    Rv(865)=Kv(865)*Ngv(32,18);
    
%   R(1,number_of_reaction+1:number_of_reactions)=Rv(1,1:number_of_vreaction);
    
end
%% Wall losses
W(1:nsp)=0;
W(8)=-1/10*area_volume*n(8)*sqrt(1.6e-19*Te/mass(8));
W(12)=-1/10*area_volume*n(12)*sqrt(1.6e-19*Te/mass(12));
W(14)=-1/10*area_volume*n(14)*sqrt(1.6e-19*Te/mass(14));
W(16)=-1/10*area_volume*n(16)*sqrt(1.6e-19*Te/mass(16));
W(17)=-1/10*area_volume*n(17)*sqrt(1.6e-19*Te/mass(17));
%W(18)=-1/2*area_volume*n(18)*sqrt(1.6e-19*Te/mass(18));
W(19)=-1/10*area_volume*n(19)*sqrt(1.6e-19*Te/mass(19));
W(20)=-1/10*area_volume*n(20)*sqrt(1.6e-19*Te/mass(20));
W(21)=-1/10*area_volume*n(21)*sqrt(1.6e-19*Te/mass(21));
W(28)=-1/10*area_volume*n(28)*sqrt(1.6e-19*Te/mass(28));
W(29)=-1/10*area_volume*n(29)*sqrt(1.6e-19*Te/mass(29));
W(30)=-1/10*area_volume*n(30)*sqrt(1.6e-19*Te/mass(30));
W(31)=-1/10*area_volume*n(31)*sqrt(1.6e-19*Te/mass(31));
if sim==1
    W([14 17 19 29 30 31])=0;
end
W(2)=-W(14)-2*W(19)-W(20);
W(3)=-W(12)-W(18)-0.01*W(16);
W(4)=-0.98*W(16);
W(10)=-W(8)-W(18)-W(20)-2*W(21);
W(22)=-W(28);
W(23)=-W(29);
W(24)=-W(30);
W(25)=-W(31);
W(11)=-W(17);
W(1)=W(8)+W(12)+W(14)+W(16)+W(17)+W(18)+W(19)+W(20)+W(21)+W(28)+W(29)+W(30)+W(31);
Wv=0;
if Vibration==1
    Wv(1:18)=-0.023.*nv(1:18);
%     Wv(7:10)=-0.023.*nv(7:10);
%     Wv(14:17)=-0.023.*nv(14:17);
    W(1,length(n)+1:nsp)=Wv;
    %     W(19)=-1/20*area_volume*n(19)*sqrt(8*1.6e-19*Tg/(pi*mass(19)));%借鉴comsol中关于Cl和Ar的解释，黏附系数取1，下同
    %     W(20)=-1/20*area_volume*n(20)*sqrt(8*1.6e-19*Tg/(pi*mass(20)));
    %     W(21)=-1/20*area_volume*n(21)*sqrt(8*1.6e-19*Tg/(pi*mass(21)));
    %     W(22)=-1/20*area_volume*n(22)*sqrt(8*1.6e-19*Tg/(pi*mass(22)));
    %     W(23)=-1/20*area_volume*n(23)*sqrt(8*1.6e-19*Tg/(pi*mass(23)));
    %     W(24)=-1/20*area_volume*n(24)*sqrt(8*1.6e-19*Tg/(pi*mass(24)));
    %     W(25)=-1/20*area_volume*n(25)*sqrt(8*1.6e-19*Tg/(pi*mass(25)));
    %     W(26)=-1/20*area_volume*n(26)*sqrt(8*1.6e-19*Tg/(pi*mass(26)));
    %     W(27)=-1/20*area_volume*n(27)*sqrt(8*1.6e-19*Tg/(pi*mass(27)));
    %     W(28)=-1/20*area_volume*n(28)*sqrt(8*1.6e-19*Tg/(pi*mass(28)));
    %     W(2)=W(2)+K(111)*W(1)*dt*W(25)*dt+K(112)*W(1)*dt*W(26)*dt+K(113)*W(1)*dt*W(27)*dt+K(114)*W(1)*dt*W(28)*dt;
    %     W(4)=W(4)+K(111)*W(1)*dt*W(25)*dt+K(112)*W(1)*dt*W(26)*dt+K(113)*W(1)*dt*W(27)*dt+K(114)*W(1)*dt*W(28)*dt;
    %     W(10)=W(10)-W(19)-W(20)-W(21)-W(22)-W(23)-W(24)-W(25)-W(26)-W(27)-W(28)-(K(111)*W(1)*dt*W(25)*dt+K(112)*W(1)*dt*W(26)*dt+K(113)*W(1)*dt*W(27)*dt+K(114)*W(1)*dt*W(28)*dt);
end
if sim==1
    Wv([6,12,13,18])=0;
    W(1,length(n)+1:nsp)=Wv;
end
W(10)=W(10)-sum(Wv);
if W(1)<=-0.1*n(1)/dt
    W(2:length(n))=-W(2:length(n))*0.1*n(1)/dt/W(1);
    W(1)=-0.1*n(1)/dt;
end

% e: species 1
isp=1;
X(isp,1)=R(1);
X(isp,3)=-R(3);
X(isp,6)=-R(6);
X(isp,7)=R(7);
X(isp,8)=R(8);
X(isp,9)=-R(9);
X(isp,12)=R(12);
X(isp,15)=R(15);
X(isp,16)=R(16);
X(isp,17)=-R(17);
X(isp,19)=-R(19);
X(isp,22)=R(22);
X(isp,23)=-R(23);
X(isp,29)=R(29);
X(isp,32)=R(32);
X(isp,33)=R(33);
X(isp,34)=-R(34);
X(isp,37)=-R(37);
X(isp,39)=R(39);
X(isp,44)=R(44);
X(isp,45)=R(45);
X(isp,48)=-R(48);
X(isp,53)=R(53);
X(isp,54)=R(54);
X(isp,58)=R(58);
X(isp,59)=R(59);
X(isp,68)=R(68);
X(isp,69)=-R(69);
X(isp,71)=-R(71);
X(isp,73)=-R(73);
X(isp,74)=-R(74);
X(isp,75)=-R(75);
X(isp,76)=-R(76);
X(isp,95)=R(95);
X(isp,97)=-R(97);
X(isp,98)=R(98);
X(isp,99)=-R(99);
X(isp,100)=R(100);
X(isp,101)=R(101);
X(isp,102)=R(102);
X(isp,106)=-R(106);
X(isp,107)=-R(107);
X(isp,108)=R(108);
X(isp,109)=R(109);
X(isp,110)=R(110);
X(isp,111)=R(111);
X(isp,112)=R(112);
X(isp,113)=R(113);
X(isp,114)=R(114);
X(isp,115)=R(115);
X(isp,116)=R(116);
X(isp,117)=-R(117);
X(isp,118)=-R(118);
X(isp,119)=-R(119);
X(isp,120)=-R(120);
X(isp,177)=R(177);
X(isp,178)=R(178);
Xgv(isp,3:6)=Vibration.*Rv(3:6);
Xgv(isp,9:10)=Vibration.*Rv(9:10);
Xgv(isp,35:38)=Vibration.*Rv(35:38);
Xgv(isp,43:46)=-Vibration.*Rv(43:46);
Xgv(isp,47:58)=Vibration.*Rv(47:58);
Xgv(isp,107:110)=-Vibration.*Rv(107:110);
Xgv(isp,111:126)=Vibration.*Rv(111:126);
Xgv(isp,223:224)=-Vibration.*Rv(223:224);
Xgv(isp,235:244)=Vibration.*Rv(235:244);
Xgv(isp,275:284)=-Vibration.*Rv(275:284);
Xgv(isp,327:329)=-Vibration.*Rv(327:329);
Xgv(isp,339:350)=Vibration.*Rv(339:350);
Xgv(isp,354:356)=-Vibration.*Rv(354:356);
Xgv(isp,472:474)=Vibration.*Rv(472:474);
Xgv(isp,478:490)=-Vibration.*Rv(478:490);
Xgv(isp,499:514)=Vibration.*Rv(499:514);
Xgv(isp,519:522)=-Vibration.*Rv(519:522);
Xgv(isp,609:612)=-Vibration.*Rv(609:612);
Xgv(isp,653:656)=Vibration.*Rv(653:656);
Xgv(isp,661:664)=-Vibration.*Rv(661:664);
Xgv(isp,666:669)=Vibration.*Rv(666:669);
Xgv(isp,671)=-Vibration.*Rv(671);
Xgv(isp,689)=-Vibration.*Rv(689);
Xgv(isp,700)=Vibration.*Rv(700);
Xgv(isp,702)=-Vibration.*Rv(702);
Xgv(isp,848:865)=Vibration.*Rv(848:865);

% if(~only_e_balance)
    % CO: species 2
    isp=2;
    X(isp,2)=R(2);
    X(isp,3)=R(3);
    X(isp,7)=-R(7);
    X(isp,9)=R(9);
    X(isp,10)=R(10);
    X(isp,17)=-R(17);
    X(isp,18)=-R(18);
    X(isp,20)=2*R(20);
    X(isp,21)=R(21);
    X(isp,29)=-R(29);
    X(isp,35)=R(35);
    X(isp,45)=-R(45);
    X(isp,46)=R(46);
    X(isp,53)=R(53);
    X(isp,55)=R(55);
    X(isp,58)=-R(58);
    X(isp,60)=R(60);
    X(isp,68)=R(68);
    X(isp,74)=2*R(74);
    X(isp,75)=R(75);
    X(isp,78)=-R(78);
    X(isp,79)=-R(79);
    X(isp,80)=-R(80);
    X(isp,81)=-R(81);
    X(isp,82)=2*R(82);
    X(isp,83)=R(83);
    X(isp,84)=2*R(84);
    X(isp,85)=2*R(85);
    X(isp,86)=2*R(86);
    X(isp,87)=R(87);
    X(isp,88)=R(88);
    X(isp,89)=R(89);
    X(isp,124)=-R(124);
    X(isp,152)=R(152);
    X(isp,153)=R(153);
    X(isp,154)=R(154);
    X(isp,170)=R(170);
    X(isp,171)=R(171);
    X(isp,172)=2*R(172);
    X(isp,173)=2*R(173);
    X(isp,174)=R(174);
    X(isp,175)=-R(175);
    Xgv(isp,7:10)=Vibration.*Rv(7:10);
    Xgv(isp,13:14)=Vibration.*Rv(13:14);
    Xgv(isp,39:46)=Vibration.*Rv(39:46);
    Xgv(isp,55:58)=Vibration.*Rv(55:58);
    Xgv(isp,87:110)=Vibration.*Rv(87:110);
    Xgv(isp,123:126)=Vibration.*Rv(123:126);
    Xgv(isp,127:136)=2*Vibration.*Rv(127:136);
    Xgv(isp,167:174)=Vibration.*Rv(167:174);
    Xgv(isp,185:194)=Vibration.*Rv(185:194);
    Xgv(isp,203:210)=Vibration.*Rv(203:210);
    Xgv(isp,221:224)=Vibration.*Rv(221:224);
    Xgv(isp,348:356)=Vibration.*Rv(348:356);
    Xgv(isp,410:426)=Vibration.*Rv(410:426);
    Xgv(isp,427:429)=2*Vibration.*Rv(427:429);
    Xgv(isp,442:444)=Vibration.*Rv(442:444);
    Xgv(isp,448:450)=Vibration.*Rv(448:450);
    Xgv(isp,454:456)=Vibration.*Rv(454:456);
    Xgv(isp,511:522)=Vibration.*Rv(511:522);
    Xgv(isp,577:588)=Vibration.*Rv(577:588);
    Xgv(isp,589:592)=2*Vibration.*Rv(589:592);
    Xgv(isp,613:616)=Vibration.*Rv(613:616);
    Xgv(isp,621:624)=Vibration.*Rv(621:624);
    Xgv(isp,629:632)=Vibration.*Rv(629:632);
    Xgv(isp,669:671)=Vibration.*Rv(669:671);
    Xgv(isp,681:683)=Vibration.*Rv(681:683);
    Xgv(isp,684)=2*Vibration.*Rv(684);
    Xgv(isp,690)=Vibration.*Rv(690);
    Xgv(isp,692)=Vibration.*Rv(692);
    Xgv(isp,694)=Vibration.*Rv(694);
    Xgv(isp,722:757)=-Vibration.*Rv(722:757);
    Xgv(isp,758:775)=Vibration.*Rv(758:775);
    Xgv(isp,812:829)=Vibration.*Rv(812:829);
    
    % O2: species 3
    isp=3;
    X(isp,4)=R(4);
    X(isp,5)=-R(5);
    X(isp,6)=-R(6);
    X(isp,8)=R(8)+R(8);
    X(isp,10)=R(10);
    X(isp,11)=-R(11);
    X(isp,14)=R(14);
    X(isp,16)=R(16);
    X(isp,19)=R(19);
    X(isp,21)=-R(21);
    X(isp,22)=-R(22);
    X(isp,24)=R(24);
    X(isp,26)=2*R(26);
    X(isp,27)=R(27);
    X(isp,34)=R(34);
    X(isp,36)=-R(36);
    X(isp,38)=R(38);
    X(isp,42)=R(42);
    X(isp,49)=R(49);
    X(isp,51)=R(51);
    X(isp,52)=R(52);
    X(isp,61)=R(61);
    X(isp,63)=R(63);
    X(isp,64)=2*R(64);
    X(isp,65)=R(65);
    X(isp,66)=R(66);
    X(isp,69)=R(69);
    X(isp,71)=R(71);
    X(isp,73)=-R(73);
    X(isp,82)=-R(82);
    X(isp,85)=R(85);
    X(isp,86)=R(86);
    X(isp,88)=R(88);
    X(isp,89)=R(89);
    X(isp,92)=R(92);
    X(isp,93)=R(93);
    X(isp,113)=R(113);
    X(isp,114)=R(114);
    X(isp,119)=R(119);
    X(isp,122)=-R(122);
    X(isp,127)=R(127);
    X(isp,128)=R(128);
    X(isp,129)=R(129);
    X(isp,135)=R(135);
    X(isp,137)=R(137);
    X(isp,139)=R(139);
    X(isp,141)=R(141);
    X(isp,146)=-R(146);
    X(isp,147)=R(147);
    X(isp,155)=R(155);
    X(isp,158)=R(158);
    X(isp,163)=R(163);
    Xgv(isp,95:102)=Vibration.*Rv(95:102);
    Xgv(isp,175:184)=Vibration.*Rv(175:184);
    Xgv(isp,225:234)=-Vibration.*Rv(225:234);
    Xgv(isp,255:264)=Vibration.*Rv(255:264);
    Xgv(isp,275:284)=Vibration.*Rv(275:284);
    Xgv(isp,327:329)=Vibration.*Rv(327:329);
    Xgv(isp,424:426)=Vibration.*Rv(424:426);
    Xgv(isp,445:447)=Vibration.*Rv(445:447);
    Xgv(isp,469:471)=-Vibration.*Rv(469:471);
    Xgv(isp,475:477)=Vibration.*Rv(475:477);
    Xgv(isp,478:490)=-Vibration.*Rv(478:490);
    Xgv(isp,585:588)=Vibration.*Rv(585:588);
    Xgv(isp,609:612)=Vibration.*Rv(609:612);
    Xgv(isp,617:620)=Vibration.*Rv(617:620);
    Xgv(isp,649:652)=-Vibration.*Rv(649:652);
    Xgv(isp,657:660)=Vibration.*Rv(657:660);
    Xgv(isp,661:664)=-Vibration.*Rv(661:664);
    Xgv(isp,683)=Vibration.*Rv(683);
    Xgv(isp,689)=Vibration.*Rv(689);
    Xgv(isp,691)=Vibration.*Rv(691);
    Xgv(isp,699)=-Vibration.*Rv(699);
    Xgv(isp,701)=Vibration.*Rv(701);
    Xgv(isp,702)=-Vibration.*Rv(702);
    
    % O: species 4
    isp=4;
    X(isp,2)=R(2);
    X(isp,4)=R(4);
    X(isp,5)=R(5)+R(5);
    X(isp,6)=R(6);
    X(isp,9)=R(9);
    X(isp,10)=R(10);
    X(isp,11)=-R(11);
    X(isp,12)=R(12);
    X(isp,14)=-R(14);
    X(isp,15)=R(15);
    X(isp,16)=-R(16);
    X(isp,18)=R(18);
    X(isp,21)=R(21);
    X(isp,23)=2*R(23);
    X(isp,24)=R(24);
    X(isp,25)=3*R(25);
    X(isp,27)=2*R(27);
    X(isp,30)=R(30);
    X(isp,31)=-R(31);
    X(isp,33)=R(33);
    X(isp,35)=-R(35);
    X(isp,37)=R(37);
    X(isp,38)=R(38);
    X(isp,39)=-R(39);
    X(isp,42)=-2*R(42);
    X(isp,43)=R(43);
    X(isp,44)=R(44);
    X(isp,47)=R(47);
    X(isp,50)=R(50);
    X(isp,52)=-R(52);
    X(isp,54)=-R(54);
    X(isp,56)=R(56);
    X(isp,57)=-R(57);
    X(isp,58)=R(58);
    X(isp,65)=-R(65);
    X(isp,66)=-R(66);
    X(isp,67)=-R(67);
    X(isp,84)=R(84);
    X(isp,87)=R(87);
    X(isp,91)=R(91);
    X(isp,102)=R(102);
    X(isp,103)=R(103);
    X(isp,104)=R(104);
    X(isp,105)=R(105);
    X(isp,112)=-R(112);
    X(isp,114)=-R(114);
    X(isp,117)=R(117);
    X(isp,120)=R(120);
    X(isp,126)=R(126);
    X(isp,133)=R(133);
    X(isp,138)=R(138);
    X(isp,152)=R(152);
    X(isp,153)=R(153);
    X(isp,154)=R(154);
    X(isp,156)=R(156);
    X(isp,157)=R(157);
    X(isp,158)=R(158);
    X(isp,159)=R(159);
    X(isp,160)=2*R(160);
    X(isp,161)=2*R(161);
    X(isp,162)=R(162);
    X(isp,166)=2*R(166);
    X(isp,167)=R(167);
    X(isp,169)=R(169);
    X(isp,173)=R(173);
    X(isp,176)=-R(176);
    Xgv(isp,5:6)=Vibration.*Rv(5:6);
    Xgv(isp,11:12)=Vibration.*Rv(11:12);
    Xgv(isp,39:42)=Vibration.*Rv(39:42);
    Xgv(isp,51:54)=Vibration.*Rv(51:54);
    Xgv(isp,87:94)=Vibration.*Rv(87:94);
    Xgv(isp,95:102)=-Vibration.*Rv(95:102);
    Xgv(isp,103:106)=Vibration.*Rv(103:106);
    Xgv(isp,119:122)=Vibration.*Rv(119:122);
    Xgv(isp,195:202)=Vibration.*Rv(195:202);
    Xgv(isp,221:222)=Vibration.*Rv(221:222);
    Xgv(isp,225:234)=-Vibration.*Rv(225:234);
    Xgv(isp,235:244)=Vibration.*Rv(235:244);
    Xgv(isp,255:264)=-2*Vibration.*Rv(255:264);
    Xgv(isp,345:347)=Vibration.*Rv(345:347);
    Xgv(isp,351:353)=Vibration.*Rv(351:353);
    Xgv(isp,410:423)=Vibration.*Rv(410:423);
    Xgv(isp,424:426)=-Vibration.*Rv(424:426);
    Xgv(isp,451:453)=Vibration.*Rv(451:453);
    Xgv(isp,469:471)=-Vibration.*Rv(469:471);
    Xgv(isp,472:474)=Vibration.*Rv(472:474);
    Xgv(isp,475:477)=-2*Vibration.*Rv(475:477);
    Xgv(isp,507:510)=Vibration.*Rv(507:510);
    Xgv(isp,515:518)=Vibration.*Rv(515:518);
    Xgv(isp,577:584)=Vibration.*Rv(577:584);
    Xgv(isp,585:588)=-Vibration.*Rv(585:588);
    Xgv(isp,625:628)=Vibration.*Rv(625:628);
    Xgv(isp,649:652)=-Vibration.*Rv(649:652);
    Xgv(isp,653:656)=Vibration.*Rv(653:656);
    Xgv(isp,657:660)=-2*Vibration.*Rv(657:660);
    Xgv(isp,668)=Vibration.*Rv(668);
    Xgv(isp,670)=Vibration.*Rv(670);
    Xgv(isp,681:682)=Vibration.*Rv(681:682);
    Xgv(isp,683)=-Vibration.*Rv(683);
    Xgv(isp,693)=Vibration.*Rv(693);
    Xgv(isp,699)=-Vibration.*Rv(699);
    Xgv(isp,700)=Vibration.*Rv(700);
    Xgv(isp,701)=-2*Vibration.*Rv(701);
    Xgv(isp,830:847)=-Vibration.*Rv(830:847);
    
    % O3: species 5
    isp=5;
    X(isp,4)=-R(4);
    X(isp,8)=-R(8);
    X(isp,11)=R(11);
    X(isp,33)=-R(33);
    X(isp,34)=-R(34);
    X(isp,37)=-R(37);
    X(isp,39)=R(39);
    X(isp,48)=-R(48);
    X(isp,50)=-R(50);
    X(isp,51)=-R(51);
    X(isp,115)=R(115);
    X(isp,129)=-R(129);
    X(isp,140)=R(140);
    X(isp,142)=R(142);
    X(isp,147)=-R(147);
    Xgv(isp,225:234)=Vibration.*Rv(225:234);
    Xgv(isp,469:471)=Vibration.*Rv(469:471);
    Xgv(isp,649:652)=Vibration.*Rv(649:652);
    Xgv(isp,699)=Vibration.*Rv(699);
    
    % O—: species 6
    isp=6;
    X(isp,3)=R(3);
    X(isp,6)=R(6);
    X(isp,7)=-R(7);
    X(isp,8)=-R(8);
    X(isp,12)=-R(12);
    X(isp,14)=R(14);
    X(isp,15)=-R(15);
    X(isp,16)=-R(16);
    X(isp,17)=R(17);
    X(isp,24)=-R(24);
    X(isp,25)=-R(25);
    X(isp,28)=-R(28);
    X(isp,34)=R(34);
    X(isp,40)=-R(40);
    X(isp,41)=-R(41);
    X(isp,43)=-R(43);
    X(isp,50)=-R(50);
    X(isp,66)=R(66);
    X(isp,68)=-R(68);
    X(isp,97)=R(97);
    X(isp,99)=R(99);
    X(isp,108)=-R(108);
    X(isp,109)=-R(109);
    X(isp,110)=-R(110);
    X(isp,132)=-R(132);
    X(isp,133)=-R(133);
    X(isp,134)=-R(134);
    X(isp,151)=-R(151);
    X(isp,159)=-R(159);
    X(isp,167)=-R(167);
    Xgv(isp,43:46)=Vibration.*Rv(43:46);
    Xgv(isp,107:110)=Vibration.*Rv(107:110);
    Xgv(isp,137:166)=-Vibration.*Rv(137:166);
    Xgv(isp,223:224)=Vibration.*Rv(223:224);
    Xgv(isp,235:254)=-Vibration.*Rv(235:254);
    Xgv(isp,354:356)=Vibration.*Rv(354:356);
    Xgv(isp,430:441)=-Vibration.*Rv(430:441);
    Xgv(isp,472:474)=-Vibration.*Rv(472:474);
    Xgv(isp,519:522)=Vibration.*Rv(519:522);
    Xgv(isp,593:608)=-Vibration.*Rv(593:608);
    Xgv(isp,653:656)=-Vibration.*Rv(653:656);
    Xgv(isp,671)=Vibration.*Rv(671);
    Xgv(isp,685:688)=-Vibration.*Rv(685:688);
    Xgv(isp,700)=-Vibration.*Rv(700);
     
    % O2-: species 7
    isp=7;
    X(isp,10)=-R(10);
    X(isp,14)=-R(14);
    X(isp,26)=-R(26);
    X(isp,27)=-R(27);
    X(isp,31)=R(31);
    X(isp,37)=R(37);
    X(isp,39)=-R(39);
    X(isp,51)=-R(51);
    X(isp,52)=R(52);
    X(isp,62)=-R(62);
    X(isp,72)=-R(72);
    X(isp,73)=R(73);
    X(isp,86)=-R(86);
    X(isp,89)=-R(89);
    X(isp,93)=-R(93);
    X(isp,111)=-R(111);
    X(isp,132)=R(132);
    X(isp,135)=-R(135);
    X(isp,136)=-R(136);
    X(isp,137)=-R(137);
    X(isp,146)=R(146);
    X(isp,158)=-R(158);
    X(isp,166)=-R(166);
    Xgv(isp,211:220)=-Vibration.*Rv(211:220);
    Xgv(isp,285:304)=-Vibration.*Rv(285:304);
    Xgv(isp,460:468)=-Vibration.*Rv(460:468);
    Xgv(isp,478:490)=Vibration.*Rv(478:490);
    Xgv(isp,637:648)=-Vibration.*Rv(637:648);
    Xgv(isp,661:664)=Vibration.*Rv(661:664);
    Xgv(isp,696:698)=-Vibration.*Rv(696:698);
    Xgv(isp,702)=Vibration.*Rv(702);
    
    % CO2+; species 8
    isp=8;
    X(isp,1)=R(1);
    X(isp,9)=-R(9);
    X(isp,10)=-R(10);
    X(isp,19)=-R(19);
    X(isp,30)=-R(30);
    X(isp,35)=-R(35);
    X(isp,36)=-R(36);
    X(isp,43)=-R(43);
    X(isp,46)=R(46);
    X(isp,56)=R(56);
    X(isp,57)=-R(57);
    X(isp,63)=-R(63);
    X(isp,77)=-R(77);
    X(isp,90)=R(90);
    X(isp,121)=-R(121);
    X(isp,125)=R(125);
    X(isp,152)=-R(152);
    X(isp,153)=-R(153);
    Xgv(isp,3:4)=Vibration.*Rv(3:4);
    Xgv(isp,7:8)=Vibration.*Rv(7:8);
    Xgv(isp,11:12)=Vibration.*Rv(11:12);
    Xgv(isp,35:38)=Vibration.*Rv(35:38);
    Xgv(isp,111:114)=Vibration.*Rv(111:114);
    Xgv(isp,167:174)=Vibration.*Rv(167:174);
    Xgv(isp,195:202)=Vibration.*Rv(195:202);
    Xgv(isp,339:341)=Vibration.*Rv(339:341);
    Xgv(isp,442:444)=Vibration.*Rv(442:444);
    Xgv(isp,451:453)=Vibration.*Rv(451:453);
    Xgv(isp,499:502)=Vibration.*Rv(499:502);
    Xgv(isp,613:616)=Vibration.*Rv(613:616);
    Xgv(isp,625:628)=Vibration.*Rv(625:628);
    Xgv(isp,666)=Vibration.*Rv(666);
    Xgv(isp,690)=Vibration.*Rv(690);
    Xgv(isp,693)=Vibration.*Rv(693);
    Xgv(isp,704:721)=-Vibration.*Rv(704:721);
    Xgv(isp,776:811)=Vibration.*Rv(776:811);
    
    % CO4-; species 9
    isp=9;
    X(isp,62)=R(62);
    X(isp,63)=-R(63);
    X(isp,64)=-R(64);
    X(isp,65)=-R(65);
    X(isp,66)=-R(66);
    X(isp,67)=-R(67);
    X(isp,72)=R(72);
    X(isp,85)=-R(85);
    X(isp,88)=-R(88);
    X(isp,92)=-R(92);
    X(isp,145)=-R(145);
    X(isp,155)=-R(155);
    X(isp,163)=-R(163);
    Xgv(isp,211:220)=Vibration.*Rv(211:220);
    Xgv(isp,285:304)=Vibration.*Rv(285:304);
    Xgv(isp,460:468)=Vibration.*Rv(460:468);
    Xgv(isp,637:648)=Vibration.*Rv(637:648);
    Xgv(isp,696:698)=Vibration.*Rv(696:698);
    
    % CO2; species 10
    isp=10;
    X(isp,1)=-R(1);
    X(isp,2)=-R(2);
    X(isp,3)=-R(3);
    X(isp,7)=R(7);
    X(isp,20)=-R(20);
    X(isp,28)=-R(28);
    X(isp,29)=2*R(29);
    X(isp,30)=2*R(30);
    X(isp,31)=R(31);
    X(isp,32)=-R(32);
    X(isp,36)=R(36);
    X(isp,38)=R(38);
    X(isp,40)=-R(40);
    X(isp,41)=-R(41);
    X(isp,43)=R(43);
    X(isp,44)=-R(44);
    X(isp,46)=-R(46);
    X(isp,49)=-R(49);
    X(isp,53)=-R(53);
    X(isp,55)=-R(55);
    X(isp,56)=-R(56);
    X(isp,57)=R(57);
    X(isp,60)=-R(60);
    X(isp,62)=-R(62);
    X(isp,63)=2*R(63);
    X(isp,64)=R(64);
    X(isp,66)=R(66);
    X(isp,67)=R(67);
    X(isp,69)=R(69);
    X(isp,70)=-R(70);
    X(isp,72)=-R(72);
    X(isp,75)=R(75);
    X(isp,76)=2*R(76);
    X(isp,77)=-R(77);
    X(isp,78)=R(78);
    X(isp,79)=R(79);
    X(isp,80)=R(80);
    X(isp,81)=R(81);
    X(isp,84)=R(84);
    X(isp,85)=R(85);
    X(isp,87)=2*R(87);
    X(isp,88)=2*R(88);
    X(isp,89)=R(89);
    X(isp,90)=R(90);
    X(isp,92)=3*R(92);
    X(isp,93)=3*R(93);
    X(isp,94)=2*R(94);
    X(isp,121)=R(121);
    X(isp,125)=-R(125);
    X(isp,143)=R(143);
    X(isp,144)=R(144);
    X(isp,145)=R(145);
    X(isp,155)=R(155);
    X(isp,162)=R(162);
    X(isp,163)=R(163);
    X(isp,168)=2*R(168);
    X(isp,169)=2*R(169);
    X(isp,170)=R(170);
    X(isp,171)=R(171);
    X(isp,174)=-R(174);
    X(isp,175)=R(175);
    Xgv(isp,1:2)=-Vibration.*Rv(1:2);
    Xgv(isp,15:18)=-Vibration.*Rv(15:18);
    Xgv(isp,19:24)=Vibration.*Rv(19:24);
    Xgv(isp,31:34)=-Vibration.*Rv(31:34);
    Xgv(isp,77:85)=-Vibration.*Rv(77:85);
    Xgv(isp,245:254)=-Vibration.*Rv(245:254);
    Xgv(isp,295:304)=-Vibration.*Rv(295:304);
    Xgv(isp,317:318)=-Vibration.*Rv(317:318);
    Xgv(isp,321:323)=Vibration.*Rv(321:323);
    Xgv(isp,336:338)=-Vibration.*Rv(336:338);
    Xgv(isp,393:394)=-Vibration.*Rv(393:394);
    Xgv(isp,395)=-Vibration.*Rv(395);
    Xgv(isp,399)=-Vibration.*Rv(399);
    Xgv(isp,404)=-Vibration.*Rv(404);
    Xgv(isp,410:417)=-Vibration.*Rv(410:417);
    Xgv(isp,421:423)=-Vibration.*Rv(421:423);
    Xgv(isp,439:441)=-Vibration.*Rv(439:441);
    Xgv(isp,466:468)=-Vibration.*Rv(466:468);
    Xgv(isp,491:494)=-Vibration.*Rv(491:494);
    Xgv(isp,535:542)=-Vibration.*Rv(535:542);
    Xgv(isp,543)=Vibration.*Rv(543);
    Xgv(isp,577:580)=-Vibration.*Rv(577:580);
    Xgv(isp,605:608)=-Vibration.*Rv(605:608);
    Xgv(isp,645:648)=-Vibration.*Rv(645:648);
    Xgv(isp,665)=-Vibration.*Rv(665);
    Xgv(isp,681)=-Vibration.*Rv(681);
    Xgv(isp,688)=-Vibration.*Rv(688);
    Xgv(isp,698)=-Vibration.*Rv(698);
    Xgv(isp,703)=-Vibration.*Rv(703);
    Xgv(isp,722:757)=Vibration.*Rv(722:757);
    Xgv(isp,776:793)=Vibration.*Rv(776:793);
    
    % C: species 11
    isp=11;
    X(isp,17)=R(17);
    X(isp,18)=R(18);
    X(isp,19)=R(19);
    X(isp,20)=-R(20);
    X(isp,21)=-R(21);
    X(isp,32)=R(32);
    X(isp,47)=-R(47);
    X(isp,59)=-R(59);
    X(isp,61)=-R(61);
    X(isp,68)=-R(68);
    Xgv(isp,47:50)=Vibration.*Rv(47:50);
    Xgv(isp,115:118)=Vibration.*Rv(115:118);
    Xgv(isp,127:136)=-Vibration.*Rv(127:136);
    Xgv(isp,342:344)=Vibration.*Rv(342:344);
    Xgv(isp,427:429)=-Vibration.*Rv(427:429);
    Xgv(isp,503:506)=Vibration.*Rv(503:506);
    Xgv(isp,589:592)=-Vibration.*Rv(589:592);
    Xgv(isp,667)=Vibration.*Rv(667);
    Xgv(isp,684)=-Vibration.*Rv(684);
    
    % O2+: species 12
    isp=12;
    X(isp,22)=R(22);
    X(isp,23)=-R(23);
    X(isp,24)=-R(24);
    X(isp,25)=-R(25);
    X(isp,26)=-R(26);
    X(isp,27)=-R(27);
    X(isp,32)=R(32);
    X(isp,33)=R(33);
    X(isp,35)=R(35);
    X(isp,36)=R(36);
    X(isp,38)=-R(38);
    X(isp,47)=-R(47);
    X(isp,55)=R(55);
    X(isp,61)=-R(61);
    X(isp,64)=-R(64);
    X(isp,70)=-R(70);
    X(isp,71)=-R(71);
    X(isp,82)=R(82);
    X(isp,122)=R(122);
    X(isp,126)=-R(126);
    X(isp,127)=-R(127);
    X(isp,128)=-R(128);
    X(isp,160)=-R(160);
    X(isp,161)=-R(161);
    Xgv(isp,47:50)=Vibration.*Rv(47:50);
    Xgv(isp,115:118)=Vibration.*Rv(115:118);
    Xgv(isp,185:194)=Vibration.*Rv(185:194);
    Xgv(isp,265:284)=-Vibration.*Rv(265:284);
    Xgv(isp,327:329)=-Vibration.*Rv(327:329);
    Xgv(isp,342:344)=Vibration.*Rv(342:344);
    Xgv(isp,448:450)=Vibration.*Rv(448:450);
    Xgv(isp,457:459)=-Vibration.*Rv(457:459);
    Xgv(isp,503:506)=Vibration.*Rv(503:506);
    Xgv(isp,609:612)=-Vibration.*Rv(609:612);
    Xgv(isp,621:624)=Vibration.*Rv(621:624);
    Xgv(isp,633:636)=-Vibration.*Rv(633:636);
    Xgv(isp,667)=Vibration.*Rv(667);
    Xgv(isp,689)=-Vibration.*Rv(689);
    Xgv(isp,692)=Vibration.*Rv(692);
    Xgv(isp,695)=Vibration.*Rv(695);
    
    % CO3-: species 13
    isp=13;
    X(isp,28)=R(28);
    X(isp,29)=-R(29);
    X(isp,30)=-R(30);
    X(isp,31)=-R(31);
    X(isp,38)=-R(38);
    X(isp,40)=R(40);
    X(isp,41)=R(41);
    X(isp,49)=R(49);
    X(isp,65)=R(65);
    X(isp,84)=-R(84);
    X(isp,87)=-R(87);
    X(isp,91)=-R(91);
    X(isp,143)=-R(143);
    X(isp,144)=-R(144);
    X(isp,154)=-R(154);
    X(isp,162)=-R(162);
    Xgv(isp,137:166)=Vibration.*Rv(137:166);
    Xgv(isp,175:184)=Vibration.*Rv(175:184);
    Xgv(isp,245:254)=Vibration.*Rv(245:254);
    Xgv(isp,430:441)=Vibration.*Rv(430:441);
    Xgv(isp,445:447)=Vibration.*Rv(445:447);
    Xgv(isp,593:608)=Vibration.*Rv(593:608);
    Xgv(isp,617:620)=Vibration.*Rv(617:620);
    Xgv(isp,685:688)=Vibration.*Rv(685:688);
    Xgv(isp,691)=Vibration.*Rv(691);
        
    % CO+: species 14
    isp=14;
    X(isp,44)=R(44);
    X(isp,45)=R(45);
    X(isp,46)=-R(46);
    X(isp,47)=R(47);
    X(isp,60)=R(60);
    X(isp,83)=R(83);
    X(isp,124)=R(124);
    Xgv(isp,5:6)=Vibration.*Rv(5:6);
    Xgv(isp,7:8)=-Vibration.*Rv(7:8);
    Xgv(isp,13:14)=Vibration.*Rv(13:14);
    Xgv(isp,51:54)=Vibration.*Rv(51:54);
    Xgv(isp,119:122)=Vibration.*Rv(119:122);
    Xgv(isp,167:174)=-Vibration.*Rv(167:174);
    Xgv(isp,203:210)=Vibration.*Rv(203:210);
    Xgv(isp,345:347)=Vibration.*Rv(345:347);
    Xgv(isp,442:444)=-Vibration.*Rv(442:444);
    Xgv(isp,454:456)=Vibration.*Rv(454:456);
    Xgv(isp,507:510)=Vibration.*Rv(507:510);
    Xgv(isp,613:616)=-Vibration.*Rv(613:616);
    Xgv(isp,629:632)=Vibration.*Rv(629:632);
    Xgv(isp,668)=Vibration.*Rv(668);
    Xgv(isp,690)=-Vibration.*Rv(690);
    Xgv(isp,694)=Vibration.*Rv(694);
    Xgv(isp,758:775)=Vibration.*Rv(758:775);
    
    % O3-: species 15
    isp=15;
    X(isp,48)=R(48);
    X(isp,49)=-R(49);
    X(isp,50)=R(50);
    X(isp,51)=R(51);
    X(isp,52)=-R(52);
    X(isp,67)=R(67);
    X(isp,136)=R(136);
    X(isp,138)=-R(138);
    X(isp,139)=-R(139);
    X(isp,140)=-R(140);
    X(isp,141)=-R(141);
    X(isp,142)=-R(142);
    Xgv(isp,175:184)=-Vibration.*Rv(175:184);
    Xgv(isp,445:447)=-Vibration.*Rv(445:447);
    Xgv(isp,617:620)=-Vibration.*Rv(617:620);
    Xgv(isp,691)=-Vibration.*Rv(691);

    % O+: species 16
    isp=16;
    X(isp,53)=R(53);
    X(isp,54)=R(54);
    X(isp,55)=-R(55);
    X(isp,56)=-R(56);
    X(isp,57)=R(57);
    Xgv(isp,9:10)=Vibration.*Rv(9:10);
    Xgv(isp,11:12)=-Vibration.*Rv(11:12);
    Xgv(isp,55:58)=Vibration.*Rv(55:58);
    Xgv(isp,123:126)=Vibration.*Rv(123:126);
    Xgv(isp,185:202)=-Vibration.*Rv(185:202);
    Xgv(isp,348:350)=Vibration.*Rv(348:350);
    Xgv(isp,448:450)=-Vibration.*Rv(448:450);
    Xgv(isp,451:453)=-Vibration.*Rv(451:453);
    Xgv(isp,511:514)=Vibration.*Rv(511:514);
    Xgv(isp,621:628)=-Vibration.*Rv(621:628);
    Xgv(isp,669)=Vibration.*Rv(669);
    Xgv(isp,692)=-Vibration.*Rv(692);
    Xgv(isp,693)=-Vibration.*Rv(693);
    
    % C+: species 17
    isp=17;
    X(isp,58)=R(58);
    X(isp,59)=R(59);
    X(isp,60)=-R(60);
    X(isp,61)=R(61);
    Xgv(isp,13:14)=-Vibration.*Rv(13:14);
    Xgv(isp,203:210)=-Vibration.*Rv(203:210);
    Xgv(isp,454:456)=-Vibration.*Rv(454:456);
    Xgv(isp,629:632)=-Vibration.*Rv(629:632);
    Xgv(isp,694)=-Vibration.*Rv(694);

    % CO4+: species 18
    isp=18;
    X(isp,69)=-R(69);
    X(isp,70)=R(70);
    Xgv(isp,265:274)=Vibration.*Rv(265:274);
    Xgv(isp,457:459)=Vibration.*Rv(457:459);
    Xgv(isp,633:636)=Vibration.*Rv(633:636);
    Xgv(isp,695)=Vibration.*Rv(695);
    
    % C2O2+: species 19
    isp=19;
    X(isp,74)=-R(74);
    X(isp,78)=R(78);
    X(isp,80)=R(80);
    X(isp,82)=-R(82);
    X(isp,83)=-R(83);
    X(isp,84)=-R(84);
    X(isp,85)=-R(85);
    X(isp,86)=-R(86);
    X(isp,172)=-R(172);
    X(isp,173)=-R(173);
    Xgv(isp,722:739)=Vibration.*Rv(722:739);
    Xgv(isp,758:775)=-Vibration.*Rv(758:775);
    
    % C2O3+: species 20
    isp=20;
    X(isp,75)=-R(75);
    X(isp,78)=-R(78);
    X(isp,79)=R(79);
    X(isp,80)=-R(80);
    X(isp,81)=R(81);
    X(isp,87)=-R(87);
    X(isp,88)=-R(88);
    X(isp,89)=-R(89);
    X(isp,170)=-R(170);
    X(isp,171)=-R(171);
    Xgv(isp,722:739)=-Vibration.*Rv(722:739);
    Xgv(isp,740:757)=Vibration.*Rv(740:757);
    
    % C2O4+: species 21
    isp=21;
    X(isp,76)=-R(76);
    X(isp,77)=R(77);
    X(isp,79)=-R(79);
    X(isp,81)=-R(81);
    X(isp,90)=-R(90);
    X(isp,91)=-R(91);
    X(isp,92)=-R(92);
    X(isp,93)=-R(93);
    X(isp,168)=-R(168);
    X(isp,169)=-R(169);
    Xgv(isp,704:721)=Vibration.*Rv(704:721);
    Xgv(isp,740:757)=-Vibration.*Rv(740:757);
    Xgv(isp,776:793)=-Vibration.*Rv(776:793);
    
    % N2: species 22
    isp=22;
    X(isp,95)=-R(95);
    X(isp,96)=-R(96);
    X(isp,99)=R(99);
    X(isp,105)=R(105);
    X(isp,109)=-R(109);
    X(isp,113)=R(113);
    X(isp,115)=R(115);
    X(isp,120)=R(120);
    X(isp,122)=R(122);
    X(isp,123)=R(123);
    X(isp,124)=R(124);
    X(isp,125)=R(125);
    X(isp,136)=R(136);
    X(isp,162)=R(162);
    X(isp,163)=R(163);
    X(isp,166)=R(166);
    X(isp,167)=R(167);
    Xgv(isp,794:811)=Vibration.*Rv(794:811);
    
    % NO: species 23
    isp=23;
    X(isp,97)=-R(97);
    X(isp,98)=-R(98);
    X(isp,103)=-R(103);
    X(isp,104)=R(104);
    X(isp,108)=R(108);
    X(isp,110)=-R(110);
    X(isp,116)=R(116);
    X(isp,121)=-R(121);
    X(isp,123)=-R(123);
    X(isp,127)=-R(127);
    X(isp,130)=-R(130);
    X(isp,131)=-R(131);
    X(isp,132)=R(132);
    X(isp,134)=R(134);
    X(isp,138)=-R(138);
    X(isp,139)=-R(139);
    X(isp,143)=-R(143);
    X(isp,145)=-R(145);
    X(isp,146)=R(146);
    X(isp,148)=R(148);
    X(isp,150)=-R(150);
    X(isp,151)=-R(151);
    X(isp,154)=R(154);
    X(isp,155)=R(155);
    X(isp,159)=R(159);
    X(isp,174)=R(174);
    X(isp,175)=R(175);
    X(isp,176)=-R(176);
    X(isp,177)=R(177);
    X(isp,178)=R(178);
    Xgv(isp,812:829)=Vibration.*Rv(812:829);
    Xgv(isp,830:847)=-Vibration.*Rv(830:847);
    Xgv(isp,848:865)=Vibration.*Rv(848:865);
    
    % N2O: species 24
    isp=24;
    X(isp,99)=-R(99);
    X(isp,100)=-R(100);
    X(isp,105)=-R(105);
    X(isp,109)=R(109);
    X(isp,131)=R(131);
    X(isp,134)=-R(134);
    X(isp,136)=-R(136);
    
    % NO2: species 25
    isp=25;
    X(isp,101)=-R(101);
    X(isp,102)=-R(102);
    X(isp,104)=-R(104);
    X(isp,106)=-R(106);
    X(isp,107)=-R(107);
    X(isp,110)=R(110);
    X(isp,111)=R(111);
    X(isp,114)=R(114);
    X(isp,128)=-R(128);
    X(isp,130)=R(130);
    X(isp,132)=-R(132);
    X(isp,133)=-R(133);
    X(isp,135)=-R(135);
    X(isp,140)=-R(140);
    X(isp,141)=-R(141);
    X(isp,144)=-R(144);
    X(isp,148)=-R(148);
    X(isp,149)=R(149);
    X(isp,150)=R(150);
    X(isp,152)=R(152);
    X(isp,156)=R(156);
    X(isp,160)=R(160);
    X(isp,164)=R(164);
    X(isp,168)=R(168);
    X(isp,169)=R(169);
    X(isp,170)=R(170);
    X(isp,171)=R(171);
    X(isp,172)=R(172);
    X(isp,173)=R(173);
    X(isp,175)=-R(175);
    X(isp,176)=R(176);
    Xgv(isp,830:847)=Vibration.*Rv(830:847);
    
    % N: species 26
    isp=26;
    X(isp,96)=2*R(96);
    X(isp,97)=R(97);
    X(isp,103)=R(103);
    X(isp,108)=-R(108);
    X(isp,111)=-R(111);
    X(isp,113)=-R(113);
    X(isp,115)=-R(115);
    X(isp,117)=R(117);
    X(isp,118)=2*R(118);
    X(isp,119)=R(119);
    X(isp,126)=-R(126);
    X(isp,156)=R(156);
    X(isp,157)=R(157);
    X(isp,158)=R(158);
    X(isp,164)=2*R(164);
    X(isp,165)=2*R(165);
    X(isp,174)=-R(174);
    Xgv(isp,812:829)=-Vibration.*Rv(812:829);
    
    % NO3: species 27
    isp=27;
    X(isp,112)=R(112);
    X(isp,137)=-R(137);
    X(isp,142)=-R(142);
    X(isp,149)=-R(149);
    X(isp,153)=R(153);
    X(isp,157)=R(157);
    X(isp,161)=R(161);
    X(isp,165)=R(165);
    
    % N2+: species 28
    isp=28;
    X(isp,95)=R(95);
    X(isp,118)=-R(118);
    X(isp,122)=-R(122);
    X(isp,123)=-R(123);
    X(isp,124)=-R(124);
    X(isp,125)=-R(125);
    X(isp,162)=-R(162);
    X(isp,163)=-R(163);
    X(isp,164)=-R(164);
    X(isp,165)=-R(165);
    X(isp,166)=-R(166);
    X(isp,167)=-R(167);
    Xgv(isp,794:811)=-Vibration.*Rv(794:811);
    
    % NO+: species 29
    isp=29;
    X(isp,98)=R(98);
    X(isp,102)=R(102);
    X(isp,117)=-R(117);
    X(isp,121)=R(121);
    X(isp,123)=R(123);
    X(isp,126)=R(126);
    X(isp,127)=R(127);
    X(isp,129)=-R(129);
    X(isp,130)=R(130);
    X(isp,131)=R(131);
    X(isp,154)=-R(154);
    X(isp,155)=-R(155);
    X(isp,156)=-R(156);
    X(isp,157)=-R(157);
    X(isp,158)=-R(158);
    X(isp,159)=-R(159);
    
    % N2O+: species 30
    isp=30;
    X(isp,100)=R(100);
    X(isp,120)=-R(120);
    X(isp,131)=-R(131);
    
    % NO2+: species 31
    isp=31;
    X(isp,101)=R(101);
    X(isp,119)=-R(119);
    X(isp,128)=R(128);
    X(isp,129)=R(129);
    X(isp,130)=-R(130);
    X(isp,164)=-R(164);
    
    % NO-: species 32
    isp=32;
    X(isp,116)=-R(116);
    X(isp,134)=R(134);
    X(isp,146)=-R(146);
    X(isp,177)=-R(177);
    X(isp,178)=-R(178);
    Xgv(isp,848:865)=-Vibration.*Rv(848:865);
    
    % NO2-: species 33
    isp=33;
    X(isp,106)=R(106);
    X(isp,107)=R(107);
    X(isp,112)=-R(112);
    X(isp,113)=-R(113);
    X(isp,133)=R(133);
    X(isp,135)=R(135);
    X(isp,139)=R(139);
    X(isp,140)=R(140);
    X(isp,143)=R(143);
    X(isp,147)=-R(147);
    X(isp,148)=-R(148);
    X(isp,149)=-R(149);
    X(isp,150)=R(150);
    X(isp,151)=R(151);
    X(isp,152)=-R(152);
    X(isp,156)=-R(156);
    X(isp,160)=-R(160);
    X(isp,168)=-R(168);
    X(isp,170)=-R(170);
    X(isp,172)=-R(172);
    
    % NO3-: species 34
    isp=34;
    X(isp,114)=-R(114);
    X(isp,115)=-R(115);
    X(isp,137)=R(137);
    X(isp,138)=R(138);
    X(isp,141)=R(141);
    X(isp,142)=R(142);
    X(isp,144)=R(144);
    X(isp,145)=R(145);
    X(isp,147)=R(147);
    X(isp,148)=R(148);
    X(isp,149)=R(149);
    X(isp,150)=-R(150);
    X(isp,153)=-R(153);
    X(isp,157)=-R(157);
    X(isp,161)=-R(161);
    X(isp,165)=-R(165);
    X(isp,169)=-R(169);
    X(isp,171)=-R(171);
    X(isp,173)=-R(173);
   
if Vibration==1
    % CO2e1: species 35, vspecies 1
    isp=1;
    Xvv(isp,1)=Rv(1);
    Xvv(isp,3)=-Rv(3);
    Xvv(isp,5)=-Rv(5);
    Xvv(isp,7)=-Rv(7);
    Xvv(isp,9)=-Rv(9);
    Xvv(isp,11)=-Rv(11);
    Xvv(isp,13)=-Rv(13);
    Xvv(isp,127)=-Rv(127);
    Xvv(isp,137)=-Rv(137);
    Xvv(isp,147)=-Rv(147);
    Xvv(isp,157)=-Rv(157);
    Xvv(isp,175)=-Rv(175);
    Xvv(isp,185)=-Rv(185);
    Xvv(isp,211)=-Rv(211);
    Xvv(isp,221)=-Rv(221);
    Xvv(isp,223)=-Rv(223);
    Xvv(isp,265)=-Rv(265);
    Xvv(isp,285)=-Rv(285);
    Xvv(isp,704)=-Rv(704);
    Xvv(isp,794)=-Rv(794);
    Xvv(isp,812)=-Rv(812);
    
    % CO2e2: species 36, vspecies 2
    isp=2;
    Xvv(isp,2)=Rv(2);
    Xvv(isp,4)=-Rv(4);
    Xvv(isp,6)=-Rv(6);
    Xvv(isp,8)=-Rv(8);
    Xvv(isp,10)=-Rv(10);
    Xvv(isp,12)=-Rv(12);
    Xvv(isp,14)=-Rv(14);
    Xvv(isp,128)=-Rv(128);
    Xvv(isp,138)=-Rv(138);
    Xvv(isp,148)=-Rv(148);
    Xvv(isp,158)=-Rv(158);
    Xvv(isp,176)=-Rv(176);
    Xvv(isp,186)=-Rv(186);
    Xvv(isp,212)=-Rv(212);
    Xvv(isp,222)=-Rv(222);
    Xvv(isp,224)=-Rv(224);
    Xvv(isp,266)=-Rv(266);
    Xvv(isp,286)=-Rv(286);
    Xvv(isp,705)=-Rv(705);
    Xvv(isp,795)=-Rv(795);
    Xvv(isp,813)=-Rv(813);
    
    % CO2va: species 37, vspecies 3
    isp=3;
    Xvv(isp,15)=Rv(15);
    Xvv(isp,19)=-Rv(19);
    Xvv(isp,20)=-Rv(20);
    Xvv(isp,21)=-Rv(21);
    Xvv(isp,25:30)=Rv(25:30);
    Xvv(isp,59)=Rv(59);
    Xvv(isp,60)=Rv(60);
    Xvv(isp,61)=Rv(61);
    Xvv(isp,77)=Rv(77);
    Xvv(isp,78)=Rv(78);
    Xvv(isp,80)=Rv(80);
    Xvv(isp,82)=Rv(82);
    Xvv(isp,84)=2*Rv(84);
    Xvv(isp,85)=Rv(85);
    Xvv(isp,86)=-Rv(86);
    Xvv(isp,87)=-Rv(87);
    Xvv(isp,95)=-Rv(95);
    Xvv(isp,103)=-Rv(103);
    Xvv(isp,107)=-Rv(107);
    Xvv(isp,111)=-Rv(111);
    Xvv(isp,115)=-Rv(115);
    Xvv(isp,119)=-Rv(119);
    Xvv(isp,123)=-Rv(123);
    Xvv(isp,129)=-Rv(129);
    Xvv(isp,139)=-Rv(139);
    Xvv(isp,149)=-Rv(149);
    Xvv(isp,159)=-Rv(159);
    Xvv(isp,167)=-Rv(167);
    Xvv(isp,177)=-Rv(177);
    Xvv(isp,187)=-Rv(187);
    Xvv(isp,195)=-Rv(195);
    Xvv(isp,203)=-Rv(203);
    Xvv(isp,213)=-Rv(213);
    Xvv(isp,267)=-Rv(267);
    Xvv(isp,287)=-Rv(287);
    Xvv(isp,318)=Rv(318);
    Xvv(isp,319)=-Rv(319);
    Xvv(isp,394:395)=Rv(394:395);
    Xvv(isp,396)=-Rv(396);
    Xvv(isp,399)=Rv(399);
    Xvv(isp,400)=-Rv(400);
    Xvv(isp,404)=Rv(404);
    Xvv(isp,405)=-Rv(405);
    Xvv(isp,535)=Rv(535);
    Xvv(isp,537)=Rv(537);
    Xvv(isp,539)=Rv(539);
    Xvv(isp,541)=Rv(541);
    Xvv(isp,675:677)=Rv(675:677);
    Xvv(isp,703)=Rv(703);
    Xvv(isp,706)=-Rv(706);
    Xvv(isp,796)=-Rv(796);
    Xvv(isp,814)=-Rv(814);
   
    % CO2vb: species 38, vspecies 4
    isp=4;
    Xvv(isp,16)=Rv(16);
    Xvv(isp,22:27)=-Rv(22:27);
    Xvv(isp,62:64)=Rv(62:64);
    Xvv(isp,77)=Rv(77);
    Xvv(isp,79)=Rv(79);
    Xvv(isp,81)=Rv(81);
    Xvv(isp,83)=Rv(83);
    Xvv(isp,84)=-Rv(84);
    Xvv(isp,85)=Rv(85);
    Xvv(isp,86)=2*Rv(86);
    Xvv(isp,88)=-Rv(88);
    Xvv(isp,96)=-Rv(96);
    Xvv(isp,104)=-Rv(104);
    Xvv(isp,108)=-Rv(108);
    Xvv(isp,112)=-Rv(112);
    Xvv(isp,116)=-Rv(116);
    Xvv(isp,120)=-Rv(120);
    Xvv(isp,124)=-Rv(124);
    Xvv(isp,130)=-Rv(130);
    Xvv(isp,140)=-Rv(140);
    Xvv(isp,150)=-Rv(150);
    Xvv(isp,160)=-Rv(160);
    Xvv(isp,168)=-Rv(168);
    Xvv(isp,178)=-Rv(178);
    Xvv(isp,188)=-Rv(188);
    Xvv(isp,196)=-Rv(196);
    Xvv(isp,204)=-Rv(204);
    Xvv(isp,214)=-Rv(214);
    Xvv(isp,268)=-Rv(268);
    Xvv(isp,288)=-Rv(288);
    Xvv(isp,305:310)=Rv(305:310);
    Xvv(isp,317)=2*Rv(317);
    Xvv(isp,319)=Rv(319);
    Xvv(isp,320)=-Rv(320);
    Xvv(isp,393)=Rv(393);
    Xvv(isp,396)=Rv(396);
    Xvv(isp,397)=-Rv(397);
    Xvv(isp,400)=Rv(400);
    Xvv(isp,401)=-Rv(401);
    Xvv(isp,405)=Rv(405);
    Xvv(isp,406)=-Rv(406);
    Xvv(isp,536)=Rv(536);
    Xvv(isp,538)=Rv(538);
    Xvv(isp,540)=Rv(540);
    Xvv(isp,542)=Rv(542);
    Xvv(isp,678:680)=Rv(678:680);
    Xvv(isp,703)=Rv(703);
    Xvv(isp,707)=-Rv(707);
    Xvv(isp,797)=-Rv(797);
    Xvv(isp,815)=-Rv(815);
    
    % CO2vc: species 39, vspecies 5
    isp=5;
    Xvv(isp,17)=Rv(17);
    Xvv(isp,28)=-Rv(28);
    Xvv(isp,29)=-Rv(29);
    Xvv(isp,30)=-Rv(30);
    Xvv(isp,65)=Rv(65);
    Xvv(isp,66)=Rv(66);
    Xvv(isp,67)=Rv(67);
    Xvv(isp,85:86)=-Rv(85:86);
    Xvv(isp,89)=-Rv(89);
    Xvv(isp,97)=-Rv(97);
    Xvv(isp,105)=-Rv(105);
    Xvv(isp,109)=-Rv(109);
    Xvv(isp,113)=-Rv(113);
    Xvv(isp,117)=-Rv(117);
    Xvv(isp,121)=-Rv(121);
    Xvv(isp,125)=-Rv(125);
    Xvv(isp,131)=-Rv(131);
    Xvv(isp,141)=-Rv(141);
    Xvv(isp,151)=-Rv(151);
    Xvv(isp,161)=-Rv(161);
    Xvv(isp,169)=-Rv(169);
    Xvv(isp,179)=-Rv(179);
    Xvv(isp,189)=-Rv(189);
    Xvv(isp,197)=-Rv(197);
    Xvv(isp,205)=-Rv(205);
    Xvv(isp,215)=-Rv(215);
    Xvv(isp,269)=-Rv(269);
    Xvv(isp,289)=-Rv(289);
    Xvv(isp,305:307)=-Rv(305:307);
    Xvv(isp,311:313)=Rv(311:313);
    Xvv(isp,318:319)=Rv(318:319);
    Xvv(isp,320)=2*Rv(320);
    Xvv(isp,357:359)=Rv(357:359);
    Xvv(isp,366:368)=Rv(366:368);
    Xvv(isp,397)=Rv(397);
    Xvv(isp,398)=-Rv(398);
    Xvv(isp,401)=Rv(401);
    Xvv(isp,402)=-Rv(402);
    Xvv(isp,406)=Rv(406);
    Xvv(isp,407)=-Rv(407);
    Xvv(isp,708)=-Rv(708);
    Xvv(isp,798)=-Rv(798);
    Xvv(isp,816)=-Rv(816);
    
    % CO2vd: species 40, vspecies 6
    isp=6;
    Xvv(isp,18)=Rv(18);
    Xvv(isp,90)=-Rv(90);
    Xvv(isp,98)=-Rv(98);
    Xvv(isp,106)=-Rv(106);
    Xvv(isp,110)=-Rv(110);
    Xvv(isp,114)=-Rv(114);
    Xvv(isp,118)=-Rv(118);
    Xvv(isp,122)=-Rv(122);
    Xvv(isp,126)=-Rv(126);
    Xvv(isp,132)=-Rv(132);
    Xvv(isp,142)=-Rv(142);
    Xvv(isp,152)=-Rv(152);
    Xvv(isp,162)=-Rv(162);
    Xvv(isp,170)=-Rv(170);
    Xvv(isp,180)=-Rv(180);
    Xvv(isp,190)=-Rv(190);
    Xvv(isp,198)=-Rv(198);
    Xvv(isp,206)=-Rv(206);
    Xvv(isp,216)=-Rv(216);
    Xvv(isp,270)=-Rv(270);
    Xvv(isp,290)=-Rv(290);
    Xvv(isp,308:320)=-Rv(308:320);
    Xvv(isp,360:362)=Rv(360:362);
    Xvv(isp,369:371)=Rv(369:371);
    Xvv(isp,398:401)=Rv(398:401);
    Xvv(isp,402)=2*Rv(402);
    Xvv(isp,404:406)=Rv(404:406);
    Xvv(isp,407)=2*Rv(407);
    Xvv(isp,709)=-Rv(709);
    Xvv(isp,799)=-Rv(799);
    Xvv(isp,817)=-Rv(817);
    
    % CO2v1: species 41, vspecies 7
    isp=7;
    Xvv(isp,31)=Rv(31);
    Xvv(isp,35)=-Rv(35);
    Xvv(isp,39)=-Rv(39);
    Xvv(isp,43)=-Rv(43);
    Xvv(isp,47)=-Rv(47);
    Xvv(isp,51)=-Rv(51);
    Xvv(isp,55)=-Rv(55);
    Xvv(isp,59)=-Rv(59);
    Xvv(isp,60)=-Rv(60);
    Xvv(isp,61)=-Rv(61);
    Xvv(isp,62)=-Rv(62);
    Xvv(isp,63)=-Rv(63);
    Xvv(isp,64)=-Rv(64);
    Xvv(isp,65)=-Rv(65);
    Xvv(isp,66)=-Rv(66);
    Xvv(isp,67)=-Rv(67);
    Xvv(isp,68)=Rv(68);
    Xvv(isp,69)=Rv(69);
    Xvv(isp,70)=Rv(70);
    Xvv(isp,77)=-Rv(77);
    Xvv(isp,78)=Rv(78);
    Xvv(isp,79)=Rv(79);
    Xvv(isp,91)=-Rv(91);
    Xvv(isp,99)=-Rv(99);
    Xvv(isp,133)=-Rv(133);
    Xvv(isp,143)=-Rv(143);
    Xvv(isp,153)=-Rv(153);
    Xvv(isp,163)=-Rv(163);
    Xvv(isp,171)=-Rv(171);
    Xvv(isp,181)=-Rv(181);
    Xvv(isp,191)=-Rv(191);
    Xvv(isp,199)=-Rv(199);
    Xvv(isp,207)=-Rv(207);
    Xvv(isp,217)=-Rv(217);
    Xvv(isp,271)=-Rv(271);
    Xvv(isp,291)=-Rv(291);
    Xvv(isp,314:316)=Rv(314:316);
    Xvv(isp,321)=-2*Rv(321);
    Xvv(isp,322:323)=-Rv(322:323);
    Xvv(isp,324:325)=Rv(324:325);
    Xvv(isp,330:331)=-Rv(330:331);
    Xvv(isp,333)=-Rv(333);
    Xvv(isp,363:365)=Rv(363:365);
    Xvv(isp,372:374)=Rv(372:374);
    Xvv(isp,395:398)=Rv(395:398);
    Xvv(isp,403)=-Rv(403);
    Xvv(isp,408)=-Rv(408);
    Xvv(isp,543)=-Rv(543);
    Xvv(isp,544:545)=Rv(544:545);
    Xvv(isp,559:562)=-Rv(559:562);
    Xvv(isp,672:674)=-Rv(672:674);
    Xvv(isp,710)=-Rv(710);
    Xvv(isp,800)=-Rv(800);
    Xvv(isp,818)=-Rv(818);
    
    % CO2v2: species 42, vspecies 8
    isp=8;
    Xvv(isp,32)=Rv(32);
    Xvv(isp,36)=-Rv(36);
    Xvv(isp,40)=-Rv(40);
    Xvv(isp,44)=-Rv(44);
    Xvv(isp,48)=-Rv(48);
    Xvv(isp,52)=-Rv(52);
    Xvv(isp,56)=-Rv(56);
    Xvv(isp,68)=-Rv(68);
    Xvv(isp,69)=-Rv(69);
    Xvv(isp,70)=-Rv(70);
    Xvv(isp,71)=Rv(71);
    Xvv(isp,72)=Rv(72);
    Xvv(isp,73)=Rv(73);
    Xvv(isp,78)=-Rv(78);
    Xvv(isp,79)=-Rv(79);
    Xvv(isp,80)=Rv(80);
    Xvv(isp,81)=Rv(81);
    Xvv(isp,92)=-Rv(92);
    Xvv(isp,100)=-Rv(100);
    Xvv(isp,134)=-Rv(134);
    Xvv(isp,144)=-Rv(144);
    Xvv(isp,154)=-Rv(154);
    Xvv(isp,164)=-Rv(164);
    Xvv(isp,172)=-Rv(172);
    Xvv(isp,182)=-Rv(182);
    Xvv(isp,192)=-Rv(192);
    Xvv(isp,200)=-Rv(200);
    Xvv(isp,208)=-Rv(208);
    Xvv(isp,218)=-Rv(218);
    Xvv(isp,272)=-Rv(272);
    Xvv(isp,292)=-Rv(292);
    Xvv(isp,321)=Rv(321);
    Xvv(isp,322)=-Rv(322);
    Xvv(isp,324)=-2*Rv(324);
    Xvv(isp,325)=-Rv(325);
    Xvv(isp,326)=Rv(326);
    Xvv(isp,330)=Rv(330);
    Xvv(isp,331)=2*Rv(331);
    Xvv(isp,332)=-Rv(332);
    Xvv(isp,333)=Rv(333);
    Xvv(isp,334)=-Rv(334);
    Xvv(isp,384:392)=-Rv(384:392);
    Xvv(isp,393:394)=-Rv(393:394);
    Xvv(isp,544:545)=-Rv(544:545);
    Xvv(isp,546:548)=Rv(546:548);
    Xvv(isp,559:562)=Rv(559:562);
    Xvv(isp,563:566)=-Rv(563:566);
    Xvv(isp,711)=-Rv(711);
    Xvv(isp,801)=-Rv(801);
    Xvv(isp,819)=-Rv(819);
    
    % CO2v3: species 43, vspecies 9
    isp=9;
    Xvv(isp,33)=Rv(33);
    Xvv(isp,37)=-Rv(37);
    Xvv(isp,41)=-Rv(41);
    Xvv(isp,45)=-Rv(45);
    Xvv(isp,49)=-Rv(49);
    Xvv(isp,53)=-Rv(53);
    Xvv(isp,57)=-Rv(57);
    Xvv(isp,71)=-Rv(71);
    Xvv(isp,72)=-Rv(72);
    Xvv(isp,73)=-Rv(73);
    Xvv(isp,74)=Rv(74);
    Xvv(isp,75)=Rv(75);
    Xvv(isp,76)=Rv(76);
    Xvv(isp,80)=-Rv(80);
    Xvv(isp,81)=-Rv(81);
    Xvv(isp,82)=Rv(82);
    Xvv(isp,83)=Rv(83);
    Xvv(isp,93)=-Rv(93);
    Xvv(isp,101)=-Rv(101);
    Xvv(isp,135)=-Rv(135);
    Xvv(isp,145)=-Rv(145);
    Xvv(isp,155)=-Rv(155);
    Xvv(isp,165)=-Rv(165);
    Xvv(isp,173)=-Rv(173);
    Xvv(isp,183)=-Rv(183);
    Xvv(isp,193)=-Rv(193);
    Xvv(isp,201)=-Rv(201);
    Xvv(isp,209)=-Rv(209);
    Xvv(isp,219)=-Rv(219);
    Xvv(isp,273)=-Rv(273);
    Xvv(isp,293)=-Rv(293);
    Xvv(isp,322)=Rv(322);
    Xvv(isp,323)=-Rv(323);
    Xvv(isp,324)=Rv(324);
    Xvv(isp,325)=-Rv(325);
    Xvv(isp,326)=-2*Rv(326);
    Xvv(isp,330)=Rv(330);
    Xvv(isp,331)=-Rv(331);
    Xvv(isp,332)=2*Rv(332);
    Xvv(isp,334)=Rv(334);
    Xvv(isp,335)=-Rv(335);
    Xvv(isp,546:548)=-Rv(546:548);
    Xvv(isp,549:552)=Rv(549:552);
    Xvv(isp,563:566)=Rv(563:566);
    Xvv(isp,567:570)=-Rv(567:570);
    Xvv(isp,712)=-Rv(712);
    Xvv(isp,802)=-Rv(802);
    Xvv(isp,820)=-Rv(820);
    
    % CO2v4: species 44, vspecies 10
    isp=10;
    Xvv(isp,34)=Rv(34);
    Xvv(isp,38)=-Rv(38);
    Xvv(isp,42)=-Rv(42);
    Xvv(isp,46)=-Rv(46);
    Xvv(isp,50)=-Rv(50);
    Xvv(isp,54)=-Rv(54);
    Xvv(isp,58)=-Rv(58);
    Xvv(isp,74)=-Rv(74);
    Xvv(isp,75)=-Rv(75);
    Xvv(isp,76)=-Rv(76);
    Xvv(isp,82)=-Rv(82);
    Xvv(isp,83)=-Rv(83);
    Xvv(isp,94)=-Rv(94);
    Xvv(isp,102)=-Rv(102);
    Xvv(isp,136)=-Rv(136);
    Xvv(isp,146)=-Rv(146);
    Xvv(isp,156)=-Rv(156);
    Xvv(isp,166)=-Rv(166);
    Xvv(isp,174)=-Rv(174);
    Xvv(isp,184)=-Rv(184);
    Xvv(isp,194)=-Rv(194);
    Xvv(isp,202)=-Rv(202);
    Xvv(isp,210)=-Rv(210);
    Xvv(isp,220)=-Rv(220);
    Xvv(isp,274)=-Rv(274);
    Xvv(isp,294)=-Rv(294);
    Xvv(isp,323)=Rv(323);
    Xvv(isp,325:326)=Rv(325:326);
    Xvv(isp,330)=-Rv(330);
    Xvv(isp,332)=-Rv(332);
    Xvv(isp,335)=Rv(335);
    Xvv(isp,495)=-Rv(495);
    Xvv(isp,535:536)=Rv(535:536);
    Xvv(isp,543:544)=-Rv(543:544);
    Xvv(isp,546)=-Rv(546);
    Xvv(isp,549)=-2*Rv(549);
    Xvv(isp,550:552)=-Rv(550:552);
    Xvv(isp,553:555)=Rv(553:555);
    Xvv(isp,559)=Rv(559);
    Xvv(isp,563)=Rv(563);
    Xvv(isp,567)=2*Rv(567);
    Xvv(isp,568:570)=Rv(568:570);
    Xvv(isp,571:573)=-Rv(571:573);
    Xvv(isp,713)=-Rv(713);
    Xvv(isp,803)=-Rv(803);
    Xvv(isp,821)=-Rv(821);
    
    % CO2v1a: species 45, vspecies 11
    isp=11;
    Xvv(isp,336)=Rv(336);
    Xvv(isp,339)=-Rv(339);
    Xvv(isp,342)=-Rv(342);
    Xvv(isp,345)=-Rv(345);
    Xvv(isp,348)=-Rv(348);
    Xvv(isp,351)=-Rv(351);
    Xvv(isp,354)=-Rv(354);
    Xvv(isp,357:365)=-Rv(357:365);
    Xvv(isp,375:380)=Rv(375:380);
    Xvv(isp,384:386)=Rv(384:386);
    Xvv(isp,393)=Rv(393);
    Xvv(isp,395:398)=-Rv(395:398);
    Xvv(isp,403)=2*Rv(403);
    Xvv(isp,408)=Rv(408);
    Xvv(isp,409)=-Rv(409);
    Xvv(isp,418)=-Rv(418);
    Xvv(isp,424)=-Rv(424);
    Xvv(isp,427)=-Rv(427);
    Xvv(isp,430)=-Rv(430);
    Xvv(isp,433)=-Rv(433);
    Xvv(isp,436)=-Rv(436);
    Xvv(isp,442)=-Rv(442);
    Xvv(isp,445)=-Rv(445);
    Xvv(isp,448)=-Rv(448);
    Xvv(isp,451)=-Rv(451);
    Xvv(isp,454)=-Rv(454);
    Xvv(isp,457)=-Rv(457);
    Xvv(isp,460)=-Rv(460);
    Xvv(isp,463)=-Rv(463);
    Xvv(isp,714)=-Rv(714);
    Xvv(isp,804)=-Rv(804);
    Xvv(isp,822)=-Rv(822);
    
    % CO2v1b: species 46, vspecies 12
    isp=12;
    Xvv(isp,337)=Rv(337);
    Xvv(isp,340)=-Rv(340);
    Xvv(isp,343)=-Rv(343);
    Xvv(isp,346)=-Rv(346);
    Xvv(isp,349)=-Rv(349);
    Xvv(isp,352)=-Rv(352);
    Xvv(isp,355)=-Rv(355);
    Xvv(isp,366:377)=-Rv(366:377);
    Xvv(isp,381:383)=Rv(381:383);
    Xvv(isp,387:389)=Rv(387:389);
    Xvv(isp,394)=Rv(394);
    Xvv(isp,399:403)=-Rv(399:403);
    Xvv(isp,408)=Rv(408);
    Xvv(isp,409)=2*Rv(409);
    Xvv(isp,419)=-Rv(419);
    Xvv(isp,425)=-Rv(425);
    Xvv(isp,428)=-Rv(428);
    Xvv(isp,431)=-Rv(431);
    Xvv(isp,434)=-Rv(434);
    Xvv(isp,437)=-Rv(437);
    Xvv(isp,443)=-Rv(443);
    Xvv(isp,446)=-Rv(446);
    Xvv(isp,449)=-Rv(449);
    Xvv(isp,452)=-Rv(452);
    Xvv(isp,455)=-Rv(455);
    Xvv(isp,458)=-Rv(458);
    Xvv(isp,461)=-Rv(461);
    Xvv(isp,464)=-Rv(464);
    Xvv(isp,715)=-Rv(715);
    Xvv(isp,805)=-Rv(805);
    Xvv(isp,823)=-Rv(823);
    
    % CO2v1c: species 47, vspecies 13
    isp=13;
    Xvv(isp,338)=Rv(338);
    Xvv(isp,341)=-Rv(341);
    Xvv(isp,344)=-Rv(344);
    Xvv(isp,347)=-Rv(347);
    Xvv(isp,350)=-Rv(350);
    Xvv(isp,353)=-Rv(353);
    Xvv(isp,356)=-Rv(356);
    Xvv(isp,378:383)=-Rv(378:383);
    Xvv(isp,390:392)=Rv(390:392);
    Xvv(isp,404:409)=-Rv(404:409);
    Xvv(isp,420)=-Rv(420);
    Xvv(isp,426)=-Rv(426);
    Xvv(isp,429)=-Rv(429);
    Xvv(isp,432)=-Rv(432);
    Xvv(isp,435)=-Rv(435);
    Xvv(isp,438)=-Rv(438);
    Xvv(isp,444)=-Rv(444);
    Xvv(isp,447)=-Rv(447);
    Xvv(isp,450)=-Rv(450);
    Xvv(isp,453)=-Rv(453);
    Xvv(isp,456)=-Rv(456);
    Xvv(isp,459)=-Rv(459);
    Xvv(isp,462)=-Rv(462);
    Xvv(isp,465)=-Rv(465);
    Xvv(isp,716)=-Rv(716);
    Xvv(isp,806)=-Rv(806);
    Xvv(isp,824)=-Rv(824);
    
    % CO2v5: species 48, vspecies 14
    isp=14;
    Xvv(isp,491)=Rv(491);
    Xvv(isp,495)=Rv(495);
    Xvv(isp,496)=-Rv(496);
    Xvv(isp,499)=-Rv(499);
    Xvv(isp,503)=-Rv(503);
    Xvv(isp,507)=-Rv(507);
    Xvv(isp,511)=-Rv(511);
    Xvv(isp,515)=-Rv(515);
    Xvv(isp,519)=-Rv(519);
    Xvv(isp,535:536)=-Rv(535:536);
    Xvv(isp,537:538)=Rv(537:538);
    Xvv(isp,543:544)=Rv(543:544);
    Xvv(isp,545)=-Rv(545);
    Xvv(isp,546)=Rv(546);
    Xvv(isp,547)=-Rv(547);
    Xvv(isp,549)=Rv(549);
    Xvv(isp,550)=-Rv(550);
    Xvv(isp,553)=-2*Rv(553);
    Xvv(isp,554:555)=-Rv(554:555);
    Xvv(isp,556:557)=Rv(556:557);
    Xvv(isp,559)=-Rv(559);
    Xvv(isp,560)=Rv(560);
    Xvv(isp,563)=-Rv(563);
    Xvv(isp,564)=Rv(564);
    Xvv(isp,567)=-Rv(567);
    Xvv(isp,568)=Rv(568);
    Xvv(isp,571)=2*Rv(571);
    Xvv(isp,572:573)=Rv(572:573);
    Xvv(isp,574:575)=-Rv(574:575);
    Xvv(isp,581)=-Rv(581);
    Xvv(isp,585)=-Rv(585);
    Xvv(isp,589)=-Rv(589);
    Xvv(isp,593)=-Rv(593);
    Xvv(isp,597)=-Rv(597);
    Xvv(isp,601)=-Rv(601);
    Xvv(isp,613)=-Rv(613);
    Xvv(isp,617)=-Rv(617);
    Xvv(isp,621)=-Rv(621);
    Xvv(isp,625)=-Rv(625);
    Xvv(isp,629)=-Rv(629);
    Xvv(isp,633)=-Rv(633);
    Xvv(isp,637)=-Rv(637);
    Xvv(isp,641)=-Rv(641);
    Xvv(isp,717)=-Rv(717);
    Xvv(isp,807)=-Rv(807);
    Xvv(isp,825)=-Rv(825);
    
    % CO2v6: species 49, vspecies 15
    isp=15;
    Xvv(isp,492)=Rv(492);
    Xvv(isp,496)=Rv(496);
    Xvv(isp,497)=-Rv(497);
    Xvv(isp,500)=-Rv(500);
    Xvv(isp,504)=-Rv(504);
    Xvv(isp,508)=-Rv(508);
    Xvv(isp,512)=-Rv(512);
    Xvv(isp,516)=-Rv(516);
    Xvv(isp,520)=-Rv(520);
    Xvv(isp,537:538)=-Rv(537:538);
    Xvv(isp,539:540)=Rv(539:540);
    Xvv(isp,545)=Rv(545);
    Xvv(isp,547)=Rv(547);
    Xvv(isp,548)=-Rv(548);
    Xvv(isp,550)=Rv(550);
    Xvv(isp,551)=-Rv(551);
    Xvv(isp,553)=Rv(553);
    Xvv(isp,554)=-Rv(554);
    Xvv(isp,556)=-2*Rv(556);
    Xvv(isp,557)=-Rv(557);
    Xvv(isp,558)=Rv(558);
    Xvv(isp,560)=-Rv(560);
    Xvv(isp,561)=Rv(561);
    Xvv(isp,564)=-Rv(564);
    Xvv(isp,565)=Rv(565);
    Xvv(isp,568)=-Rv(568);
    Xvv(isp,569)=Rv(569);
    Xvv(isp,571)=-Rv(571);
    Xvv(isp,572)=Rv(572);
    Xvv(isp,574)=2*Rv(574);
    Xvv(isp,575)=Rv(575);
    Xvv(isp,576)=-Rv(576);
    Xvv(isp,582)=-Rv(582);
    Xvv(isp,586)=-Rv(586);
    Xvv(isp,590)=-Rv(590);
    Xvv(isp,594)=-Rv(594);
    Xvv(isp,598)=-Rv(598);
    Xvv(isp,602)=-Rv(602);
    Xvv(isp,614)=-Rv(614);
    Xvv(isp,618)=-Rv(618);
    Xvv(isp,622)=-Rv(622);
    Xvv(isp,626)=-Rv(626);
    Xvv(isp,630)=-Rv(630);
    Xvv(isp,634)=-Rv(634);
    Xvv(isp,638)=-Rv(638);
    Xvv(isp,642)=-Rv(642);
    Xvv(isp,718)=-Rv(718);
    Xvv(isp,808)=-Rv(808);
    Xvv(isp,826)=-Rv(826);
    
    % CO2v7: species 50, vspecies 16
    isp=16;
    Xvv(isp,493)=Rv(493);
    Xvv(isp,497)=Rv(497);
    Xvv(isp,498)=-Rv(498);
    Xvv(isp,501)=-Rv(501);
    Xvv(isp,505)=-Rv(505);
    Xvv(isp,509)=-Rv(509);
    Xvv(isp,513)=-Rv(513);
    Xvv(isp,517)=-Rv(517);
    Xvv(isp,521)=-Rv(521);
    Xvv(isp,539:540)=-Rv(539:540);
    Xvv(isp,541:542)=Rv(541:542);
    Xvv(isp,548)=Rv(548);
    Xvv(isp,551)=Rv(551);
    Xvv(isp,552)=-Rv(552);
    Xvv(isp,554)=Rv(554);
    Xvv(isp,555)=-Rv(555);
    Xvv(isp,556)=Rv(556);
    Xvv(isp,557)=-Rv(557);
    Xvv(isp,558)=-2*Rv(558);
    Xvv(isp,561)=-Rv(561);
    Xvv(isp,562)=Rv(562);
    Xvv(isp,565)=-Rv(565);
    Xvv(isp,566)=Rv(566);
    Xvv(isp,569)=-Rv(569);
    Xvv(isp,570)=Rv(570);
    Xvv(isp,572)=-Rv(572);
    Xvv(isp,573)=Rv(573);
    Xvv(isp,574)=-Rv(574);
    Xvv(isp,575)=Rv(575);
    Xvv(isp,576)=2*Rv(576);
    Xvv(isp,583)=-Rv(583);
    Xvv(isp,587)=-Rv(587);
    Xvv(isp,591)=-Rv(591);
    Xvv(isp,595)=-Rv(595);
    Xvv(isp,599)=-Rv(599);
    Xvv(isp,603)=-Rv(603);
    Xvv(isp,615)=-Rv(615);
    Xvv(isp,619)=-Rv(619);
    Xvv(isp,623)=-Rv(623);
    Xvv(isp,627)=-Rv(627);
    Xvv(isp,631)=-Rv(631);
    Xvv(isp,635)=-Rv(635);
    Xvv(isp,639)=-Rv(639);
    Xvv(isp,643)=-Rv(643);
    Xvv(isp,719)=-Rv(719);
    Xvv(isp,809)=-Rv(809);
    Xvv(isp,827)=-Rv(827);
    
    % CO2v8: species 51, vspecies 17
    isp=17;
    Xvv(isp,494)=Rv(494);
    Xvv(isp,498)=Rv(498);
    Xvv(isp,502)=-Rv(502);
    Xvv(isp,506)=-Rv(506);
    Xvv(isp,510)=-Rv(510);
    Xvv(isp,514)=-Rv(514);
    Xvv(isp,518)=-Rv(518);
    Xvv(isp,522)=-Rv(522);
    Xvv(isp,541:542)=-Rv(541:542);
    Xvv(isp,552)=Rv(552);
    Xvv(isp,555)=Rv(555);
    Xvv(isp,557)=Rv(557);
    Xvv(isp,558)=Rv(558);
    Xvv(isp,562)=-Rv(562);
    Xvv(isp,566)=-Rv(566);
    Xvv(isp,570)=-Rv(570);
    Xvv(isp,573)=-Rv(573);
    Xvv(isp,575)=-Rv(575);
    Xvv(isp,576)=-Rv(576);
    Xvv(isp,584)=-Rv(584);
    Xvv(isp,588)=-Rv(588);
    Xvv(isp,592)=-Rv(592);
    Xvv(isp,596)=-Rv(596);
    Xvv(isp,600)=-Rv(600);
    Xvv(isp,604)=-Rv(604);
    Xvv(isp,616)=-Rv(616);
    Xvv(isp,620)=-Rv(620);
    Xvv(isp,624)=-Rv(624);
    Xvv(isp,628)=-Rv(628);
    Xvv(isp,632)=-Rv(632);
    Xvv(isp,636)=-Rv(636);
    Xvv(isp,640)=-Rv(640);
    Xvv(isp,644)=-Rv(644);
    Xvv(isp,720)=-Rv(720);
    Xvv(isp,810)=-Rv(810);
    Xvv(isp,828)=-Rv(828);
    
    % CO2vk: species 52, vspecies 18
    isp=18;
    Xvv(isp,665)=Rv(665);
    Xvv(isp,666:671)=-Rv(666:671);
    Xvv(isp,672:674)=Rv(672:674);
    Xvv(isp,675:680)=-Rv(675:680);
    Xvv(isp,682:687)=-Rv(682:687);
    Xvv(isp,690:697)=-Rv(690:697);
    Xvv(isp,703)=-Rv(703);
    Xvv(isp,721)=-Rv(721);
    Xvv(isp,811)=-Rv(811);
    Xvv(isp,829)=-Rv(829);
end

X(1:length(n),number_of_reaction+1:number_of_reactions)=Xgv;
X(length(n)+1:nsp,number_of_reaction+1:number_of_reactions)=Xvv;

%Wall losses

X(:,number_of_reactions+1)=W';

% Rate(1)=-sum(Xgv(10,[15:18 31:34 336:338 491:494 665]));
% Rate(2)=-sum(Xgv(10,[1 2]));
% Rate(3)=-X(10,2);
% Rate(4)=-sum(X(10,[1 44 53]));
% Rate(5)=-X(10,3);
% Rate(6)=sum(Xgv(2,[39:42 103:106 351:353 515:518 670]));
% Rate(7)=sum(Xgv(1,[35:38 51:58 111:114 119:126 339:341 345:350 499:502 507:514 666 668 669]));
% Rate(8)=sum(Xgv(2,[43:46 107:110 354:356 519:522 671]));
% 
% RateK(1)=sum(Kv([15:18 31:34 336:338 491:494 665]));
% RateK(2)=sum(Kv([1 2]));
% RateK(3)=K(2);
% RateK(4)=sum(K([1 44 53]));
% RateK(5)=K(3);
% RateK(6)=sum(Kv([39:42 103:106 351:353 515:518 670]));
% RateK(7)=sum(Kv([35:38 51:58 111:114 119:126 339:341 345:350 499:502 507:514 666 668 669]));
% RateK(8)=sum(Kv([43:46 107:110 354:356 519:522 671]));

% for i=1:nsp
%     if n(i)<=1e-3 && sum(X(i,:))<=0
%         X(i,:)=0;
%     end
% end

if n(1)<=1e6 && sum(X(1,:))<0
%     inde=X(1,1:number_of_reactions)<0;
%     X(:,inde)=0;
    X(1,:)=0;
end
p(2:length(n)) = n(2:end)<=1e-3 & sum(X(2:length(n),:),2)<0;
p(length(n)+1:nsp) = nv(1:end)<=1e-3 & sum(X(length(n)+1:nsp,:),2)<0;
[~,ind]=find(p'.*X<0);
X(:,ind)=0;

%Map reactions to wavelengths
%This can be use to generate synthetic spectra
% lambda=zeros(1,length(R));
% %Example: Reaction 4 radiates at 850nm
% lambda(4)= 850; %此处非准确数值

% if sim==1
%     n([14 19 34 27 33 15 29 32 18 17 31 30])=0;
% end
% 能量变化

dS_ape=2/(3*e*n(1))*S_abs*dt;
dS_el=-2*me/mass(10)*K(13)*n(10)*(Te-Tg)*dt-2*me/mass(22)*K(94)*n(22)*(Te-Tg)*dt-Vibration*2*me/mass(10)*K(13)*sum(nv)*(Te-Tg)*dt;
dS_epro=-Te/n(1).*(sum(X(1,:)))*dt;
dS_rea=-2/3*(13.8*K(1)*n(10)+11.46*K(2)*n(10)+6*K(5)*n(3)+9.2*K(17)*n(2)+11.1*K(18)*n(2)+12.06*K(22)*n(3)+19.5*K(44)*n(10)+14.01*K(45)*n(2)+19.1*K(53)*n(10)+13.6*K(54)*n(4)+22*K(58)*n(2)+11.2*K(59)*n(11)+...
    +15.58*K(95)*n(22)+9.753*K(96)*n(22)+9.26*K(98)*n(23)+12.89*K(100)*n(24)+10*K(101)*n(25)+12.9*K(102)*n(25)+6.5*K(103)*n(23)+3.11*K(104)*n(25)+1.67*K(105)*n(24))*dt...
    -Vibration*2/3*(7*Kv(1)*n(10)+10.5*Kv(2)*n(10)+6.8*Kv(3)*nv(1)+3.2*Kv(4)*nv(2)+12.5*Kv(5)*nv(1)+9*Kv(6)*nv(2)+12.1*Kv(9)*nv(1)+8.6*Kv(10)*nv(2)+...
        +0.083*Kv(15)*n(10)+0.167*Kv(16)*n(10)+0.252*K(17)*n(10)+0.341*K(18)*n(10)+0.291*Kv(31)*n(10)+0.58*Kv(32)*n(10)+0.86*Kv(33)*n(10)+1.14*Kv(34)*n(10)+...
        +13.8*Kv(35)*nv(7)+13.8*Kv(36)*nv(8)+13.8*Kv(37)*nv(9)+13.8*Kv(38)*nv(10)+11.171*Kv(39)*nv(7)+10.88*Kv(40)*nv(8)+10.6*Kv(41)*nv(9)+10.32*Kv(42)*nv(10)+...
        +19.21*Kv(51)*nv(7)+18.92*Kv(52)*nv(8)+18.64*Kv(53)*nv(9)+18.36*Kv(54)*nv(10)+18.81*Kv(55)*nv(7)+18.52*Kv(56)*nv(8)+18.24*Kv(57)*nv(9)+17.96*Kv(58)*nv(10)+...
        +11.37*Kv(103)*nv(3)+11.295*Kv(104)*nv(4)+11.208*Kv(105)*nv(5)+11.121*Kv(106)*nv(6)+13.8*Kv(111)*nv(3)+13.8*Kv(112)*nv(4)+13.8*Kv(113)*nv(5)+13.8*Kv(114)*nv(6)+...
        +19.417*Kv(119)*nv(3)+19.335*Kv(120)*nv(4)+19.248*Kv(121)*nv(5)+19.161*Kv(122)*nv(6)+19.017*Kv(123)*nv(3)+18.933*Kv(124)*nv(4)+18.848*Kv(125)*nv(5)+18.761*Kv(126)*nv(6)+...
        +0.29*Kv(333)*nv(7)+0.28*Kv(334)*nv(8)+0.28*Kv(335)*nv(9)+0.339*Kv(336)*n(10)+0.422*Kv(337)*n(10)+0.505*Kv(338)*n(10)+13.8*Kv(339)*nv(11)+13.8*Kv(340)*nv(12)+13.8*Kv(341)*nv(13)+...
        +19.161*Kv(345)*nv(11)+19.078*Kv(346)*nv(12)+18.995*Kv(347)*nv(13)+18.761*Kv(348)*nv(11)+18.678*Kv(349)*nv(12)+18.595*Kv(350)*nv(13)+11.121*Kv(354)*nv(11)+11.038*Kv(355)*nv(12)+10.955*Kv(356)*nv(13)+...
        +1.43*Kv(491)*n(10)+1.7*Kv(492)*n(10)+1.97*Kv(493)*n(10)+2.24*Kv(494)*n(10)+0.29*Kv(495)*nv(10)+0.27*Kv(496)*nv(14)+0.27*Kv(497)*nv(15)+0.27*Kv(498)*nv(16)+...
        +13.8*K(1)*(nv(14)+nv(15)+nv(16)+nv(17))+18.07*Kv(507)*nv(14)+17.8*Kv(508)*nv(15)+17.53*Kv(509)*nv(16)+17.26*Kv(510)*nv(17)+17.67*Kv(511)*nv(14)+17.4*Kv(512)*nv(15)+17.13*Kv(513)*nv(16)+16.86*Kv(514)*nv(17)+...
        10.03*Kv(515)*nv(14)+9.76*Kv(516)*nv(15)+9.49*Kv(517)*nv(16)+9.22*Kv(518)*nv(17)+2.5*Kv(665)*n(10)+13.8*Kv(666)*nv(18)+17*Kv(668)*nv(18)+16.6*Kv(669)*nv(18)+8.96*Kv(670)*nv(18))*dt;
dS_wall=-2/(3*n(1))*(-W(1))*(Te/2*log(mass(10)/(2*pi*me))+5/2*Te)*dt;%(Te/2*log(mg/(2*pi*me))+5/2*Te)为1个电子-离子对的能量
% if Vibration==1
%     dS_el=-2*me/mass(10)*K(13)*(n(10)+n(18)+n(19)+n(20)+n(21)+n(22)+n(23)+n(24)+n(25)+n(26)+n(27))*(Te-Tg)*dt;
%     dS_rea = dS_rea - 2/3*(7*K(73)*n(10)+10.5*K(74)*n(10)+6.8*K(75)*n(19)+3.2*K(76)*n(20)+12.5*K(77)*n(19)+9*K(78)*n(20)+12.1*K(81)*n(19)+8.6*K(82)*n(20)+...
%         +0.083*K(87)*n(10)+0.167*K(88)*n(10)+0.252*K(89)*n(10)+0.341*K(90)*n(10)+0.291*K(103)*n(10)+0.58*K(104)*n(10)+0.86*K(105)*n(10)+1.14*K(106)*n(10)+...
%         +13.8*K(107)*n(25)+13.8*K(108)*n(26)+13.8*K(109)*n(27)+13.8*K(110)*n(28)+11.171*K(111)*n(25)+10.88*K(112)*n(26)+10.6*K(113)*n(27)+10.32*K(114)*n(28)+...
%         +19.21*K(123)*n(25)+18.92*K(124)*n(26)+18.64*K(125)*n(27)+18.36*K(126)*n(28)+18.81*K(127)*n(25)+18.52*K(128)*n(26)+18.24*K(129)*n(27)+17.96*K(130)*n(28)+...
%         +11.37*K(175)*n(21)+11.295*K(176)*n(22)+11.208*K(177)*n(23)+11.121*K(178)*n(24)+13.8*K(183)*n(21)+13.8*K(184)*n(22)+13.8*K(185)*n(23)+13.8*K(186)*n(24)+...
%         +19.417*K(191)*n(21)+19.335*K(192)*n(22)+19.248*K(193)*n(23)+19.161*K(194)*n(24)+19.017*K(195)*n(21)+18.933*K(196)*n(22)+18.848*K(197)*n(23)+18.761*K(198)*n(24))*dt;
% end
dTe=dS_ape+dS_el+dS_epro+dS_rea+dS_wall;

end

