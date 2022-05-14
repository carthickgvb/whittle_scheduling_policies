% Running all schemes:
clear;
tic;
%% setup
size_vec = ["large","large","small"];
trans_vec = ["unrestricted","restricted","restricted"];
for i = 1:3

%arrival_size = "large";
%transmission = "unrestricted";

arrival_size = size_vec(i);
transmission = trans_vec(i);

load('../mat files/n_M100_fixed_'+arrival_size+'_arrival.mat','n','G','M','psi','C','l','mis','sim_time');
sim_time=10000;

% restricted or unrestricted
if transmission == "unrestricted"
    PSI = realmax*ones(n,1);
else
    PSI = psi;
end

% jobs in queue
X_wc = zeros(n,1);
X_wg = zeros(n,1);
X_nsw = zeros(n,1);
X_cella = zeros(n,1);

% Average cost: cot, Transmissions, packets dropped, rapd

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

whittle_idx_cella = zeros(n,1);
Cost_cella = zeros(sim_time,1);
Transmissions_cella = zeros(sim_time,1);
Packets_dropped_cella = zeros(sim_time,1);
cot_cella = zeros(sim_time,1);
rapd_cella = zeros(sim_time,1);

%% whittle indices for clique, graph

w_idx_c = Calculate_whittleIndices(M,PSI,C,l,mis,'clique');
w_idx_g = Calculate_whittleIndices(M,PSI,C,l,mis,'graph');

%% slotted aloha

out_sa = slottedAloha(n,G,X_wc,M,psi,l,C,sim_time);
Cost_sa = out_sa.cost;
Transmissions_sa = out_sa.transmissions;
input_sa = out_sa.inputs;
Packets_dropped_sa = out_sa.packets_dropped;
cot_sa = zeros(sim_time,1);
rapd_sa = zeros(sim_time,1);

disp('Slotted ALOHA completed');

%% MWS algorithm

out_mwis = MWIS_algo(n,G,M,psi,l,C,sim_time,input_sa);
Cost_mwis = out_mwis.cost;
Transmissions_mwis = out_mwis.transmissions;
Packets_dropped_mwis = out_mwis.packets_dropped;
cpt_mwis = Cost_mwis./Transmissions_mwis;
cot_mwis = zeros(sim_time,1);
rapd_mwis = zeros(sim_time,1);

disp('MWIS completed');

%% Lyapunov algorithm

out_neelay = MWIS_algo3(n,G,M,psi,l,C,sim_time,input_sa);
Cost_neelay = out_neelay.cost;
Transmissions_neelay = out_neelay.transmissions;
Packets_dropped_neelay = out_neelay.packets_dropped;
cpt_neelay = Cost_neelay./Transmissions_neelay;
cot_neelay = zeros(sim_time,1);
rapd_neelay = zeros(sim_time,1);

disp('Lyapunov Drift algorithm completed');

%% whittle based: clique, graph, non-stationary and Cella

% for cella
[a,b]=findMIS_new(G);
rr=size(b,1);

whittle_idx_nsw_all = zeros(n,sim_time);
for q=1:sim_time
    % whittle-clique
    for m=1:n
        whittle_idx_wc(m) = w_idx_c(m,X_wc(m)+1);
    end
    active_id_wc = Active_ID(G,whittle_idx_wc);
    
    % whittle-graph
    for m=1:n
        whittle_idx_wg(m) = w_idx_g(m,X_wg(m)+1);
    end
    active_id_wg = Active_ID(G,whittle_idx_wg);
    
    % whittle-ns
    whittle_idx_nsw = EquationSolve_SP(X_nsw,whittle_idx_nsw,G,M,PSI,C,l,1);
    active_id_nsw = Active_ID(G,whittle_idx_nsw);
    whittle_idx_nsw_all(:,q) = whittle_idx_nsw;
    
    % Cella
    active_id_cella=(b(mod(q-1,rr)+1,:))';
    
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
        
        Cost_cella(q) = Cost_cella(q) + C*X_cella(m) + active_id_cella(m)*f(min(X_cella(m),PSI(m)));
        Transmissions_cella(q) = Transmissions_cella(q) + active_id_cella(m)*(min(X_cella(m),PSI(m)));
        Packets_dropped_cella(q) = Packets_dropped_cella(q) + max(0,X_cella(m) - active_id_cella(m)*(min(X_cella(m),PSI(m))) + input_sa(q,m) - M(m));
        X_cella(m) = min(X_cella(m) - active_id_cella(m)*(min(X_cella(m),PSI(m))) + input_sa(q,m),M(m));
    end
   
    %Compute the average cost
    if q == 1
        cot_wc(q) = Cost_wc(q);
        cot_wg(q) = Cost_wg(q);
        cot_nsw(q) = Cost_nsw(q);
        cot_cella(q) = Cost_cella(q);
        cot_sa(q) = Cost_sa(q);
        cot_mwis(q) = Cost_mwis(q);
        cot_neelay(q) = Cost_neelay(q);
        
        rapd_wc(q) = Packets_dropped_wc(q);
        rapd_wg(q) = Packets_dropped_wg(q);
        rapd_nsw(q) = Packets_dropped_nsw(q);
        rapd_cella(q) = Packets_dropped_cella(q);
        rapd_sa(q) = Packets_dropped_sa(q);
        rapd_mwis(q) = Packets_dropped_mwis(q);
        rapd_neelay(q) = Packets_dropped_mwis(q);
    else
        cot_wc(q) = (cot_wc(q-1)*(q-1) + Cost_wc(q))/q;
        cot_wg(q) = (cot_wg(q-1)*(q-1) + Cost_wg(q))/q;
        cot_nsw(q) = (cot_nsw(q-1)*(q-1) + Cost_nsw(q))/q;
        cot_cella(q) = (cot_cella(q-1)*(q-1) + Cost_cella(q))/q;
        cot_sa(q) = (cot_sa(q-1)*(q-1) + Cost_sa(q))/q;
        cot_mwis(q) = (cot_mwis(q-1)*(q-1) + Cost_mwis(q))/q;
        cot_neelay(q) = (cot_neelay(q-1)*(q-1) + Cost_neelay(q))/q;
        
        rapd_wc(q) = (rapd_wc(q-1)*(q-1) + Packets_dropped_wc(q))/q;
        rapd_wg(q) = (rapd_wg(q-1)*(q-1) + Packets_dropped_wg(q))/q;
        rapd_nsw(q) = (rapd_nsw(q-1)*(q-1) + Packets_dropped_nsw(q))/q;
        rapd_cella(q) = (rapd_cella(q-1)*(q-1) + Packets_dropped_cella(q))/q;
        rapd_sa(q) = (rapd_sa(q-1)*(q-1) + Packets_dropped_sa(q))/q;
        rapd_mwis(q) = (rapd_mwis(q-1)*(q-1) + Packets_dropped_mwis(q))/q;
        rapd_neelay(q) = (rapd_neelay(q-1)*(q-1) + Packets_dropped_neelay(q))/q;
    end
end
disp("All schemes done");
toc;

%% Display/store results

figure('Position', get(0, 'Screensize'));
X_axis = 1:sim_time;
plot(X_axis,cot_wc,'bo',X_axis,cot_wg,'gx',X_axis,cot_nsw,'r+',X_axis,cot_cella,'c*',X_axis,cot_sa,'ms',X_axis,cot_mwis,'yd',X_axis,cot_neelay,'kv');
legend('Clique Whittle Policy','Graphical Whittle Policy','New stationary Whittle Policy','Cella algorithm','slotted ALOHA','MWS algorithm','Lyapunov drift algorithm','Location','southeast','fontsize',16);
ylabel("Average Cost",'fontweight','bold','fontsize',16);
xlabel("Time",'fontweight','bold','fontsize',16);
title(arrival_size+" arrivals, "+transmission+" case",'fontsize',16);
saveas(gcf,"../../results/new_stationary_1/AC"+"_"+num2str(arrival_size)+"_"+num2str(transmission)+".png");

figure('Position', get(0, 'Screensize'));
X_axis = 1:sim_time;
plot(X_axis,rapd_wc,'bo',X_axis,rapd_wg,'gx',X_axis,rapd_nsw,'r+',X_axis,rapd_cella,'c*',X_axis,rapd_sa,'ms',X_axis,rapd_mwis,'yd',X_axis,rapd_neelay,'kv');
legend('Clique Whittle Policy','Graphical Whittle Policy','New stationary Whittle Policy','Cella algorithm','slotted ALOHA','MWS algorithm','Lyapunov drift algorithm','Location','southeast','fontsize',16);
ylabel("Packets Dropped",'fontweight','bold','fontsize',16);
xlabel("Time",'fontweight','bold','fontsize',16);
title(arrival_size+" arrivals, "+transmission+" case",'fontsize',16);
saveas(gcf,"../../results/new_stationary_1/PD"+"_"+num2str(arrival_size)+"_"+num2str(transmission)+".png");
end
%save("matfiles/"+arrival_size+"_"+transmission+".mat");