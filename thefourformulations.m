%% running 
% Non stationary type-1: nos1
% Non stationary type-2: nos2
% New stationary type-1: nes1
% New stationary type-2: nes2

clear; 
close all;
tic;
%% setup
size_vec = ["large","large","small","small"];
trans_vec = ["unrestricted","restricted","unrestricted","restricted"];
cost_vec = [10,20,30,50,80,100];
M_range_vec = [0,20,50,75,100,1000];
for i=2:length(cost_vec)

arrival_size = "small ";
transmission = "restricted";

%arrival_size = size_vec(i);
%transmission = trans_vec(i);

load('../mat files/n_M100_fixed_'+arrival_size+'_arrival.mat','n','G','M','psi','C','l','mis','sim_time');
%sim_time=10000;
%C = cost_vec(i);
a = M_range_vec(i-1);
b = M_range_vec(i);
M = floor(((randi(b,n,1)-1)/(b-1))*(b-a)+a);

% restricted or unrestricted
if transmission == "unrestricted"
    PSI = realmax*ones(n,1);
else
    PSI = psi;
end

% jobs in queue
X_nos1 = zeros(n,1);
X_nos2 = zeros(n,1);
X_nes1 = zeros(n,1);
X_nes2 = zeros(n,1);

% Average cost: cot, Transmissions, packets dropped, rapd

whittle_idx_nos1 = zeros(n,1);
Cost_nos1 = zeros(sim_time,1);
Transmissions_nos1 = zeros(sim_time,1);
Packets_dropped_nos1 = zeros(sim_time,1);
cot_nos1 = zeros(sim_time,1);
rapd_nos1 = zeros(sim_time,1);

whittle_idx_nos2 = zeros(n,1);
Cost_nos2 = zeros(sim_time,1);
Transmissions_nos2 = zeros(sim_time,1);
Packets_dropped_nos2 = zeros(sim_time,1);
cot_nos2 = zeros(sim_time,1);
rapd_nos2 = zeros(sim_time,1);

whittle_idx_nes1 = zeros(n,1);
Cost_nes1 = zeros(sim_time,1);
Transmissions_nes1 = zeros(sim_time,1);
Packets_dropped_nes1 = zeros(sim_time,1);
cot_nes1 = zeros(sim_time,1);
rapd_nes1 = zeros(sim_time,1);

whittle_idx_nes2 = zeros(n,1);
Cost_nes2 = zeros(sim_time,1);
Transmissions_nes2 = zeros(sim_time,1);
Packets_dropped_nes2 = zeros(sim_time,1);
cot_nes2 = zeros(sim_time,1);
rapd_nes2 = zeros(sim_time,1);

% slotted aloha

out_sa = slottedAloha(n,G,X_nos1,M,psi,l,C,sim_time);
Cost_sa = out_sa.cost;
Transmissions_sa = out_sa.transmissions;
input_sa = out_sa.inputs;
Packets_dropped_sa = out_sa.packets_dropped;
cot_sa = zeros(sim_time,1);
rapd_sa = zeros(sim_time,1);

disp('Slotted ALOHA completed');

%% simulation

%whittle_idx_nsw_all = zeros(n,sim_time);
for q=1:sim_time  
    % nos1
    whittle_idx_nos1 = EquationSolve_NSP(X_nos1,whittle_idx_nos1,G,M,PSI,C,l,1);
    active_id_nos1 = Active_ID(G,whittle_idx_nos1);
    %whittle_idx_nsw_all(:,q) = whittle_idx_nsw;
    
    % nos2
    whittle_idx_nos2 = EquationSolve_NSP(X_nos2,whittle_idx_nos2,G,M,PSI,C,l,2);
    active_id_nos2 = Active_ID(G,whittle_idx_nos2);
    
    % nes1
    whittle_idx_nes1 = EquationSolve_SP(X_nes1,whittle_idx_nes1,G,M,PSI,C,l,1);
    active_id_nes1 = Active_ID(G,whittle_idx_nes1);
    
    % nes2
    whittle_idx_nes2 = EquationSolve_SP(X_nes2,whittle_idx_nes2,G,M,PSI,C,l,2);
    active_id_nes2 = Active_ID(G,whittle_idx_nes2);

    for m=1:n
        Cost_nos1(q) = Cost_nos1(q) + C*X_nos1(m) + active_id_nos1(m)*f(min(X_nos1(m),PSI(m)));
        Transmissions_nos1(q) = Transmissions_nos1(q) + active_id_nos1(m)*(min(X_nos1(m),PSI(m)));
        Packets_dropped_nos1(q) = Packets_dropped_nos1(q) + max(0,X_nos1(m) - active_id_nos1(m)*(min(X_nos1(m),PSI(m))) + input_sa(q,m) - M(m));
        X_nos1(m) = min(X_nos1(m) - active_id_nos1(m)*(min(X_nos1(m),PSI(m))) + input_sa(q,m),M(m));

        Cost_nos2(q) = Cost_nos2(q) + C*X_nos2(m) + active_id_nos2(m)*f(min(X_nos2(m),PSI(m)));
        Transmissions_nos2(q) = Transmissions_nos2(q) + active_id_nos2(m)*(min(X_nos2(m),PSI(m)));
        Packets_dropped_nos2(q) = Packets_dropped_nos2(q) + max(0,X_nos2(m) - active_id_nos2(m)*(min(X_nos2(m),PSI(m))) + input_sa(q,m) - M(m));
        X_nos2(m) = min(X_nos2(m) - active_id_nos2(m)*(min(X_nos2(m),PSI(m))) + input_sa(q,m),M(m));
        
        Cost_nes1(q) = Cost_nes1(q) + C*X_nes1(m) + active_id_nes1(m)*f(min(X_nes1(m),PSI(m)));
        Transmissions_nes1(q) = Transmissions_nes1(q) + active_id_nes1(m)*(min(X_nes1(m),PSI(m)));
        Packets_dropped_nes1(q) = Packets_dropped_nes1(q) + max(0,X_nes1(m) - active_id_nes1(m)*(min(X_nes1(m),PSI(m))) + input_sa(q,m) - M(m));
        X_nes1(m) = min(X_nes1(m) - active_id_nes1(m)*(min(X_nes1(m),PSI(m))) + input_sa(q,m),M(m));
        
        Cost_nes2(q) = Cost_nes2(q) + C*X_nes2(m) + active_id_nes2(m)*f(min(X_nes2(m),PSI(m)));
        Transmissions_nes2(q) = Transmissions_nes2(q) + active_id_nes2(m)*(min(X_nes2(m),PSI(m)));
        Packets_dropped_nes2(q) = Packets_dropped_nes2(q) + max(0,X_nes2(m) - active_id_nes2(m)*(min(X_nes2(m),PSI(m))) + input_sa(q,m) - M(m));
        X_nes2(m) = min(X_nes2(m) - active_id_nes2(m)*(min(X_nes2(m),PSI(m))) + input_sa(q,m),M(m));
    end
   
    %Compute the average cost
    if q == 1
        cot_nos1(q) = Cost_nos1(q);
        cot_nos2(q) = Cost_nos2(q);
        cot_nes1(q) = Cost_nes1(q);
        cot_nes2(q) = Cost_nes2(q);
        cot_sa(q) = Cost_sa(q);
        
        rapd_nos1(q) = Packets_dropped_nos1(q);
        rapd_nos2(q) = Packets_dropped_nos2(q);
        rapd_nes1(q) = Packets_dropped_nes1(q);
        rapd_nes2(q) = Packets_dropped_nes2(q);
        rapd_sa(q) = Packets_dropped_sa(q);
    else
        cot_nos1(q) = (cot_nos1(q-1)*(q-1) + Cost_nos1(q))/q;
        cot_nos2(q) = (cot_nos2(q-1)*(q-1) + Cost_nos2(q))/q;
        cot_nes1(q) = (cot_nes1(q-1)*(q-1) + Cost_nes1(q))/q;
        cot_nes2(q) = (cot_nes2(q-1)*(q-1) + Cost_nes2(q))/q;
        cot_sa(q) = (cot_sa(q-1)*(q-1) + Cost_sa(q))/q;
        
        rapd_nos1(q) = (rapd_nos1(q-1)*(q-1) + Packets_dropped_nos1(q))/q;
        rapd_nos2(q) = (rapd_nos2(q-1)*(q-1) + Packets_dropped_nos2(q))/q;
        rapd_nes1(q) = (rapd_nes1(q-1)*(q-1) + Packets_dropped_nes1(q))/q;
        rapd_nes2(q) = (rapd_nes2(q-1)*(q-1) + Packets_dropped_nes2(q))/q;
        rapd_sa(q) = (rapd_sa(q-1)*(q-1) + Packets_dropped_sa(q))/q;
    end
end
disp("All schemes done");
toc;

%% Display/store results

figure('Position', get(0, 'Screensize'));
X_axis = 1:sim_time;
plot(X_axis,cot_nos1,'bo',X_axis,cot_nos2,'gx',X_axis,cot_nes1,'r+',X_axis,cot_nes2,'c*');
legend('Non stationary type-1','Non stationary type-2','New stationary type-1','New stationary type-2','Location','southeast','fontsize',16);
ylabel("Average Cost",'fontweight','bold','fontsize',16);
xlabel("Time",'fontweight','bold','fontsize',16);
title(arrival_size+" arrivals, "+transmission+"\_M\_range="+num2str(a)+"\_"+num2str(b)+" case",'fontsize',16);
saveas(gcf,"../../results/FourFormulations/AC"+"_"+num2str(arrival_size)+"_"+num2str(transmission)+"_M_range="+num2str(a)+"_"+num2str(b)+".png");

figure('Position', get(0, 'Screensize'));
X_axis = 1:sim_time;
plot(X_axis,rapd_nos1,'bo',X_axis,rapd_nos2,'gx',X_axis,rapd_nes1,'r+',X_axis,rapd_nes2,'c*');
legend('Non stationary type-1','Non stationary type-2','New stationary type-1','New stationary type-2','Location','southeast','fontsize',16);
ylabel("Packets Dropped",'fontweight','bold','fontsize',16);
xlabel("Time",'fontweight','bold','fontsize',16);
title(arrival_size+" arrivals, "+transmission+"\_M\_range="+num2str(a)+"\_"+num2str(b)+" case",'fontsize',16);
saveas(gcf,"../../results/FourFormulations/PD"+"_"+num2str(arrival_size)+"_"+num2str(transmission)+"_M_range="+num2str(a)+"_"+num2str(b)+".png");
end
%save("matfiles/"+arrival_size+"_"+transmission+".mat");