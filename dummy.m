%% Simulation
clear;

%% setup

load('n_M100_fixed_small_arrival.mat','n','G','M','psi','C','l','mis','sim_time');

% restricted or unrestricted
%PSI = realmax*ones(n,1);
PSI = psi;

%w_idx_c = Calculate_whittleIndices(M,PSI,C,l,mis,'clique');
w_idx_g = Calculate_whittleIndices(M,PSI,C,l,mis,'graph');

%initialize packets in queue's to 0 OR some random value
X_wc = zeros(n,1); %jobs in queue
X_wg = zeros(n,1);
X_nsw = zeros(n,1);

%Compute the Cost, Transmissions, input => arrivals and packets dropped per
%time slot
out_sa = slottedAloha(n,G,X_wc,M,psi,l,C,sim_time);
Cost_sa = out_sa.cost;
Transmissions_sa = out_sa.transmissions;
input_sa = out_sa.inputs;
Packets_dropped_sa = out_sa.packets_dropped;

disp('Slotted ALOHA completed');

cot_sa = zeros(sim_time,1);
rapd_sa = zeros(sim_time,1);

%run the MWIS algorithm
out_mwis = MWIS_algo(n,G,M,psi,l,C,sim_time,input_sa);
Cost_mwis = out_mwis.cost;
Transmissions_mwis = out_mwis.transmissions;
Packets_dropped_mwis = out_mwis.packets_dropped;

cpt_mwis = Cost_mwis./Transmissions_mwis;

disp('MWIS completed');

cot_mwis = zeros(sim_time,1);
rapd_mwis = zeros(sim_time,1);

whittle_idx_wc = zeros(n,1);
Cost_wc = zeros(sim_time,1);
Transmissions_wc = zeros(sim_time,1);
Packets_dropped_wc = zeros(sim_time,1);
cot_wc = zeros(sim_time,1);
rapd_wc = zeros(sim_time,1);

whittle_idx_wg = zeros(n,1);
Cost_wg = zeros(sim_time,1);
Transmissions_wg = zeros(sim_time,1);
Packets_dropped_wg = zeros(sim_time,1);
cot_wg = zeros(sim_time,1);
rapd_wg = zeros(sim_time,1);

whittle_idx_nsw = zeros(n,1);
Cost_nsw = zeros(sim_time,1);
Transmissions_nsw = zeros(sim_time,1);
Packets_dropped_nsw = zeros(sim_time,1);
cot_nsw = zeros(sim_time,1);
rapd_nsw = zeros(sim_time,1);

 %% Actual sim

%initialize packets in queue's to 0 OR some random value
X_wc = zeros(n,1); %jobs in queue
X_wg = zeros(n,1);
X_nsw = zeros(n,1);

tic;
for q=1:sim_time
    for m=1:n
        whittle_idx_wc(m) = w_idx_c(m,X_wc(m)+1);
    end
    active_id_wc = Active_ID(G,whittle_idx_wc);
    
    for m=1:n
        whittle_idx_wg(m) = w_idx_g(m,X_wg(m)+1);
    end
    active_id_wg = Active_ID(G,whittle_idx_wg);
    
    whittle_idx_nsw = EquationSolve_NSP(X_nsw,whittle_idx_nsw,G,M,PSI,C,l);
    active_id_nsw = Active_ID(G,whittle_idx_nsw);
    
    Cost_wc(q) = 0;
    Transmissions_wc(q) = 0;
    Cost_wg(q) = 0;
    Transmissions_wg(q) = 0;
    Cost_nsw(q) = 0;
    Transmissions_nsw(q) = 0;
    for m=1:n
        Cost_wc(q) = Cost_wc(q) + C*X_wc(m) + active_id_wc(m)*f(min(X_wc(m),PSI(m)));
        Transmissions_wc(q) = Transmissions_wc(q) + active_id_wc(m)*(min(X_wc(m),PSI(m)));
        Packets_dropped_wc(q) = Packets_dropped_wc(q) + max(0,X_wc(m) - active_id_wc(m)*(min(X_wc(m),PSI(m))) + input_sa(q,m) - M(m));
        X_wc(m) = min(X_wc(m) - active_id_wc(m)*(min(X_wc(m),PSI(m))) + input_sa(q,m),M(m));

        Cost_wg(q) = Cost_wg(q) + C*X_wg(m) + active_id_wg(m)*f(min(X_wg(m),PSI(m)));
        Transmissions_wg(q) = Transmissions_wg(q) + active_id_wg(m)*(min(X_wg(m),PSI(m)));
        Packets_dropped_wg(q) = Packets_dropped_wg(q) + max(0,X_wg(m) - active_id_wg(m)*(min(X_wg(m),PSI(m))) + input_sa(q,m) - M(m));
        X_wg(m) = min(X_wg(m) - active_id_wg(m)*(min(X_wg(m),PSI(m))) + input_sa(q,m),M(m));

        
        Cost_nsw(q) = Cost_nsw(q) + C*X_nsw(m) + active_id_nsw(m)*f(min(X_nsw(m),PSI(m)));
        Transmissions_nsw(q) = Transmissions_nsw(q) + active_id_nsw(m)*(min(X_nsw(m),PSI(m)));
        Packets_dropped_nsw(q) = Packets_dropped_nsw(q) + max(0,X_nsw(m) - active_id_nsw(m)*(min(X_nsw(m),PSI(m))) + input_sa(q,m) - M(m));
        X_nsw(m) = min(X_nsw(m) - active_id_nsw(m)*(min(X_nsw(m),PSI(m))) + input_sa(q,m),M(m));
    end
   
    %Compute the average cost
    if q == 1
        cot_wc(q) = Cost_wc(q);
        cot_wg(q) = Cost_wg(q);
        cot_nsw(q) = Cost_nsw(q);
        cot_sa(q) = Cost_sa(q);
        cot_mwis(q) = Cost_mwis(q);
        rapd_wc(q) = Packets_dropped_wc(q);
        rapd_wg(q) = Packets_dropped_wg(q);
        rapd_sa(q) = Packets_dropped_sa(q);
        rapd_mwis(q) = Packets_dropped_mwis(q);
        rapd_nsw(q) = Packets_dropped_nsw(q);
    else
        cot_wc(q) = (cot_wc(q-1)*(q-1) + Cost_wc(q))/q;
        cot_wg(q) = (cot_wg(q-1)*(q-1) + Cost_wg(q))/q;
        cot_nsw(q) = (cot_nsw(q-1)*(q-1) + Cost_nsw(q))/q;
        cot_sa(q) = (cot_sa(q-1)*(q-1) + Cost_sa(q))/q;
        cot_mwis(q) = (cot_mwis(q-1)*(q-1) + Cost_mwis(q))/q;
        rapd_wc(q) = (rapd_wc(q-1)*(q-1) + Packets_dropped_wc(q))/q;
        rapd_wg(q) = (rapd_wg(q-1)*(q-1) + Packets_dropped_wg(q))/q;
        rapd_sa(q) = (rapd_sa(q-1)*(q-1) + Packets_dropped_sa(q))/q;
        rapd_mwis(q) = (rapd_mwis(q-1)*(q-1) + Packets_dropped_mwis(q))/q;
        rapd_nsw(q) = (rapd_nsw(q-1)*(q-1) + Packets_dropped_nsw(q))/q;
    end
end
toc;

%% Display/store the results

figure();
X_axis = 1:sim_time;
plot(X_axis,cot_wc,'b.',X_axis,cot_wg,'go',X_axis,cot_sa,'rx',X_axis,cot_mwis,'c+',X_axis,cot_nsw,'m*');
legend('Clique Whittle Policy','Graphical Whittle Policy','slotted ALOHA','MWS algorithm','Non-stationary Whittle Policy','Location','southeast');
%title('Average Cost: n\_M100\_small\_arrival.mat, restricted');
title("Average cost");

figure();
plot(X_axis,rapd_wc,'b.',X_axis,rapd_wg,'go',X_axis,rapd_sa,'rx',X_axis,rapd_mwis,'c+',X_axis,rapd_nsw,'m*');
legend('Clique Whittle Policy','Graphical Whittle Policy','slotted ALOHA','Maximum weight policy','Non-stationary Whittle Policy','Location','southeast');
%title('Average Cost: n\_M100\_small\_arrival.mat, restricted');
title("Packets Dropped");

%plot(X_axis,cot_w);
%save('small_arrival_restricted.mat');
%save("NSP_"+num2str(file)+num2str(res_un)+".mat");
%save('C20_NSP.mat','cot_w','-ASCII','-append')
