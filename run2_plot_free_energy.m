global  lamda_e N1 Np ep epM epA za2 zp_e
% If A is charged and B is neutral, then for AAABBB 5 chemical
% bondings exist, so N1=3-1=2. If the configuration is ABABAB, then N1=0.
% However, if B < A, then the minimum (n_A-n_B-1) is not 0.
% For example,  AAABB -> ABABA=0 and AAAB -> ABAA=1.

dbstop if error

for zp_e =[-1] %  valence for polymer charged part
    for za2 =[2] %  valence for counterion
        for Np=1e4 % total length of the chain
            for lamda_e=[1]% charged fraction[0.1:0.05:0.5] [0.6:0.1:1]
                for N1=ceil(Np*lamda_e*1)-1 % number of the charged chemical bonding
                    for ep=[0] % e adspsilon_AC - strength of polymer charged-counterion interaction
                        for epM=0 % epsilon_B - strength of polymer neutral-neutral interaction
                            for epA = 0 % epsilon_A - strength of polymer charged-charged interaction
                                plot_free_energy(7,0,1); % balance2_s(lb_max,lb_min,deltalb)
                            end
                        end
                    end
                end
            end
        end
    end
end
% copy of the polt_salt_free.m function to draw as the order of parameters



