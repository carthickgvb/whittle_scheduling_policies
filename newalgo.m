% non-stationary whittle policy
clear;
tic;
%% setup

arrival_size = "small";
transmission = "unrestricted";

load('../mat files/n_M100_fixed_'+arrival_size+'_arrival.mat','n','G','M','psi','C','l','mis','sim_time');
sim_time = 10000;

% restricted or unrestricted
if transmission == "unrestricted"
    PSI = realmax*ones(n,1);
else
    PSI = psi;
end

% jobs in queue
X_nsw = zeros(n,1);


% Average cost: cot, Transmissions, packets dropped, rapd

whittle_idx_nsw = zeros(n,1);
Cost_nsw = zeros(sim_time,1);
Transmissions_nsw = zeros(sim_time,1);
Packets_dropped_nsw = zeros(sim_time,1);
cot_nsw = zeros(sim_time,1);
rapd_nsw = zeros(sim_time,1);

%% slotted aloha

out_sa = slottedAloha(n,G,X_nsw,M,psi,l,C,sim_time);
Cost_sa = out_sa.cost;
Transmissions_sa = out_sa.transmissions;
input_sa = out_sa.inputs;
Packets_dropped_sa = out_sa.packets_dropped;
cot_sa = zeros(sim_time,1);
rapd_sa = zeros(sim_time,1);

disp('Slotted ALOHA completed');

%% whittle based: non-stationary
gamma=[0.00001,0.00005,0.0001,0.0005,0.001,0.005,0.01];
%gamma=[1,2];
cot_nsw_all = zeros(sim_time,length(gamma));
for i=1:length(gamma)
    X_nsw = zeros(n,1);
    % Average cost: cot, Transmissions, packets dropped, rapd

    whittle_idx_nsw = zeros(n,1);
    Cost_nsw = zeros(sim_time,1);
    Transmissions_nsw = zeros(sim_time,1);
    Packets_dropped_nsw = zeros(sim_time,1);
    cot_nsw = zeros(sim_time,1);
    rapd_nsw = zeros(sim_time,1);

for q=1:sim_time
    % whittle-ns
    whittle_idx_nsw = EquationSolve_NSP(X_nsw,whittle_idx_nsw,G,M,PSI,C,l,gamma(i));
    active_id_nsw = Active_ID(G,whittle_idx_nsw);

    for m=1:n
        Cost_nsw(q) = Cost_nsw(q) + C*X_nsw(m) + active_id_nsw(m)*f(min(X_nsw(m),PSI(m)));
        Transmissions_nsw(q) = Transmissions_nsw(q) + active_id_nsw(m)*(min(X_nsw(m),PSI(m)));
        Packets_dropped_nsw(q) = Packets_dropped_nsw(q) + max(0,X_nsw(m) - active_id_nsw(m)*(min(X_nsw(m),PSI(m))) + input_sa(q,m) - M(m));
        X_nsw(m) = min(X_nsw(m) - active_id_nsw(m)*(min(X_nsw(m),PSI(m))) + input_sa(q,m),M(m));
    end
   
    %Compute the average cost
    if q == 1
        cot_nsw(q) = Cost_nsw(q);
        cot_sa(q) = Cost_sa(q);
     
        rapd_nsw(q) = Packets_dropped_nsw(q);
        rapd_sa(q) = Packets_dropped_sa(q);
    else
        cot_nsw(q) = (cot_nsw(q-1)*(q-1) + Cost_nsw(q))/q;
        cot_sa(q) = (cot_sa(q-1)*(q-1) + Cost_sa(q))/q;
        
        rapd_nsw(q) = (rapd_nsw(q-1)*(q-1) + Packets_dropped_nsw(q))/q;
        rapd_sa(q) = (rapd_sa(q-1)*(q-1) + Packets_dropped_sa(q))/q;
    end
end
cot_nsw_all(:,i) = cot_nsw;
end
toc;
%%
X_axis=1:sim_time;
figure('Position', get(0, 'Screensize'));
plot(X_axis,cot_nsw_all(:,1),'b*',X_axis,cot_nsw_all(:,2),'g*',X_axis,cot_nsw_all(:,3),'r*',X_axis,cot_nsw_all(:,4),'c*',X_axis,cot_nsw_all(:,5),'m*',X_axis,cot_nsw_all(:,6),'y*',X_axis,cot_nsw_all(:,7),'k*');
legend("\gamma="+num2str(gamma(1)),"\gamma="+num2str(gamma(2)),"\gamma="+num2str(gamma(3)),"\gamma="+num2str(gamma(4)),"\gamma="+num2str(gamma(5)),"\gamma="+num2str(gamma(6)),"\gamma="+num2str(gamma(7)),'Location','south','fontsize',16);
ylabel("Average Cost",'fontweight','bold','fontsize',16);
xlabel("Time",'fontweight','bold','fontsize',16);
title(arrival_size+" arrivals, "+transmission+" case","fontsize",16);
    
saveas(gcf,"../../results/final_results/gamma_var"+"_"+num2str(arrival_size)+"_"+num2str(transmission)+".png");
%%
% X_axis=1:sim_time;
% figure('Position', get(0, 'Screensize'));
% plot(X_axis,cot_nsw_all(:,1),'b*',X_axis,cot_nsw_all(:,2),'g*');
% legend("\lambda iterate without other \lambda s","\lambda iterate with other \lambda s","fontsize",16,"location","south");
% ylabel("Average Cost",'fontweight','bold','fontsize',16);
% xlabel("Time",'fontweight','bold','fontsize',16);
% title(arrival_size+" arrivals, "+transmission+" case","fontsize",16);
% % 
% saveas(gcf,"../../results/final_results/comparision"+"_"+num2str(arrival_size)+"_"+num2str(transmission)+".png");
%% Display/store results

% figure();
% X_axis = 1:sim_time;
% plot(X_axis,cot_nsw,'r+',X_axis,cot_sa,'ms');
% legend('Non-stationary Whittle Policy','slotted ALOHA','Location','southeast','fontsize',16);
% ylabel("Average Cost",'fontweight','bold','fontsize',16);
% xlabel("Time",'fontweight','bold','fontsize',16);
% title(arrival_size+" arrivals, "+transmission+" case");
% 
% figure();
% X_axis = 1:sim_time;
% plot(X_axis,rapd_nsw,'r+',X_axis,rapd_sa,'ms');
% legend('Non-stationary Whittle Policy','slotted ALOHA','Location','southeast','fontsize',16);
% ylabel("Packets Dropped",'fontweight','bold','fontsize',16);
% xlabel("Time",'fontweight','bold','fontsize',16);
% title(arrival_size+" arrivals, "+transmission+" case");