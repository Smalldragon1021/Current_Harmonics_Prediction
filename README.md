# Current_Harmonics_Prediction

This code runs in MATLAB.

creator: Xiaolong Zhang, xiaolongzh88@gmail.com

Test V6.1 voltage
Test V6.2 current with LC and Resistive load
Test V6.3 include filter inductor resistance into model
Test V6.4 motor current prediction
Test V6.5 motor current prediction
          read in ia/ib/ic/vab in seperate Excel files
Test V6.6 add motor impedance nonlinear model (2022.4)
Test V6.8 use new filter parameters and new analytical impedance model (2022.11.15)
          add R_c to reference voltage (modulation index) setting function
          change capacitor bank to Y connection
          add R_c to the calculation of Veq
Test V6.9 add capacitor's ESL (including cable inductance)
Test V6.10 Sawtooth carrier - use phase angle representation to calculate harmonics
           add deadtime PWM spectrum #9,#10,#11,#12,#13 and verified
           dead time loss is included in M calculation
Test V6.11 Use only one phase in each sequence circuit calculation and then inverse
           transform back to all phases
           Phase shift in Vs1 due to deadtime is included in M calculation by vetorial sum
                            distortion, changed deadtime expression, verified with test
Test V6.13 Improved linear circuit impedance and harmonic calculation
Test V6.14 Compare frequency-invariant method (ita0 = vsa./Zeq, line 1670-1689) vs.
           proposed method (it1 = vs1./Zeq1; it2 = vs2./Zeq2; ita = it1 + it2), constrasted plot in Figure 45, prediction error in Figure 46.
Test V6.15 Simplified harmonic synthesis using complex-number representation: Cmn_a = Amn_ph.*exp(j*Theta_a)... (line 1296-1233)
