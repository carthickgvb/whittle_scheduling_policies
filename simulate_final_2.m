% Running all schemes:
clear;
tic;
%% setup
%size_vec = ["large","large","small","small"];
%trans_vec = ["unrestricted","restricted","unrestricted","restricted"];

%for i=1:4

arrival_size = "small";
transmission = "restricted";

%arrival_size = size_vec(i);
%transmission = trans_vec(i);

load('../mat files/n_M100_fixed_'+arrival_size+'_arrival.mat','n','G','M','psi','C','l','mis','sim_time');
sim_time=10000;

% restricted or unrestricted
if transmission == "unrestricted"
    PSI = realmax*ones(n,1);
else
    PSI = psi;
end

% order
% nos1 - 1, nos2 - 2, nes1 - 3, nes2 - 4, wc - 5, wg - 6, cella - 7,
% slotted aloha - 8, MWS - 9, Lyapunov - 10
% 10 policies - 6 whittle based, 6+1 need jobs variable

% jobs in queue
X = zeros(n,7);

% Average cost: cot, Transmissions, packets dropped, rapd
whittle_idx = zeros(n,6);
active_id = zeros(n,7);
Cost = zeros(sim_time,10);
Transmissions = zeros(sim_time,10);
Packets_dropped = zeros(sim_time,10);
cot = zeros(sim_time,10);
rapd = zeros(sim_time,10);

%% whittle indices for clique, graph

w_idx_c = Calculate_whittleIndices(M,PSI,C,l,mis,'clique');
w_idx_g = Calculate_whittleIndices(M,PSI,C,l,mis,'graph');

%% slotted aloha

out_sa = slottedAloha(n,G,zeros(n,1),M,psi,l,C,sim_time);
Cost(:,8) = out_sa.cost;
Transmissions(:,8) = out_sa.transmissions;
input_sa = out_sa.inputs;
Packets_dropped(:,8) = out_sa.packets_dropped;
cot(:,8) = zeros(sim_time,1);
rapd(:,8) = zeros(sim_time,1);

disp('Slotted ALOHA completed');

%% MWS algorithm

out_mwis = MWIS_algo(n,G,M,psi,l,C,sim_time,input_sa);
Cost(:,9) = out_mwis.cost;
Transmissions(:,9) = out_mwis.transmissions;
Packets_dropped(:,9) = out_mwis.packets_dropped;
%cpt_mwis = Cost_mwis./Transmissions_mwis;
cot(:,9) = zeros(sim_time,1);
rapd(:,9) = zeros(sim_time,1);

disp('MWIS completed');

%% Lyapunov algorithm

out_neelay = MWIS_algo3(n,G,M,psi,l,C,sim_time,input_sa);
Cost(:,10) = out_neelay.cost;
Transmissions(:,10)= out_neelay.transmissions;
Packets_dropped(:,10) = out_neelay.packets_dropped;
%cpt_neelay = Cost_neelay./Transmissions_neelay;
cot(:,10) = zeros(sim_time,1);
rapd(:,10) = zeros(sim_time,1);

disp('Lyapunov Drift algorithm completed');

%% whittle based: clique, graph, non-stationary and Cella

% for cella
[a,b] = findMIS_new(G);
rr=size(b,1);

for q=1:sim_time 
    % whittle-nos1
    whittle_idx(:,1) = EquationSolve_NSP(X(:,1),whittle_idx(:,1),G,M,PSI,C,l,1);
    active_id(:,1) = Active_ID(G,whittle_idx(:,1));
    
    % whittle-nos2
    whittle_idx(:,2) = EquationSolve_NSP(X(:,2),whittle_idx(:,2),G,M,PSI,C,l,2);
    active_id(:,2) = Active_ID(G,whittle_idx(:,2));
    
    % whittle-nes1
    whittle_idx(:,3) = EquationSolve_SP(X(:,3),whittle_idx(:,3),G,M,PSI,C,l,1);
    active_id(:,3) = Active_ID(G,whittle_idx(:,3));
    
    % whittle-nes2
    whittle_idx(:,4) = EquationSolve_SP(X(:,4),whittle_idx(:,4),G,M,PSI,C,l,2);
    active_id(:,4) = Active_ID(G,whittle_idx(:,4));
    
    % whittle-clique
    for m=1:n
        whittle_idx(m,5) = w_idx_c(m,X(m,5)+1);
    end
    active_id(:,5) = Active_ID(G,whittle_idx(:,5));
    
    % whittle-graph
    for m=1:n
        whittle_idx(m,6) = w_idx_g(m,X(m,6)+1);
    end
    active_id(:,6) = Active_ID(G,whittle_idx(:,6));
    
    % Cella
    active_id(:,7)=(b(mod(q-1,rr)+1,:))';
    
    for m=1:n
        for i=1:7
            Cost(q,i) = Cost(q,i) + C*X(m,i) + active_id(m,i)*f(min(X(m,i),PSI(m)));
            Transmissions(q,i) = Transmissions(q,i) + active_id(m,i)*(min(X(m,i),PSI(m)));
            Packets_dropped(q,i) = Packets_dropped(q,i) + max(0,X(m,i) - active_id(m,i)*(min(X(m,i),PSI(m))) + input_sa(q,m) - M(m));
            X(m,i) = min(X(m,i) - active_id(m,i)*(min(X(m,i),PSI(m))) + input_sa(q,m),M(m));
        end
    end
   
    %Compute the average cost
    if q == 1
        for i=1:10
            cot(q,i) = Cost(q,i);
            rapd(q,i) = Packets_dropped(q,i);
        end 
    else
        for i=1:10
            cot(q,i) = (cot(q-1,i)*(q-1) + Cost(q,i))/q;
            rapd(q,i) = (rapd(q-1,i)*(q-1) + Packets_dropped(q,i))/q;
        end
    end
end
disp("All schemes done");

%save("new_matfiles/"+arrival_size+"_"+transmission+".mat");
%end
toc;
temp1 = cot;
temp2 = rapd;

%% Display/store results
X_axis = 1:500:sim_time;
%X_axis = 1:sim_time;

cot=temp1(1:500:end,:);
figure();
% style = ['bo','go','ro','co','mo','yo','ko','wo','b.','g.'];
plot(X_axis,cot(:,1),'ro',X_axis,cot(:,2),'rx',X_axis,cot(:,3),'r.-',X_axis,cot(:,4),'r*',X_axis,cot(:,5),'ks',X_axis,cot(:,6),'kd',X_axis,cot(:,7),'kv',X_axis,cot(:,8),'k^',X_axis,cot(:,9),'kp',X_axis,cot(:,10),'k+','MarkerSize',7);
legend('Non stationary type-1','Non stationary type-2','New stationary type-1','New stationary type-2','Clique Whittle Policy','Graphical Whittle Policy','Cella algorithm','slotted ALOHA','MWS algorithm','Lyapunov drift algorithm','NumColumns',4,'Location','southeast','fontsize',16);
ylabel("Average Cost",'fontweight','bold','fontsize',16);
xlabel("Time",'fontweight','bold','fontsize',16);
title(arrival_size+" arrivals, "+transmission+" case",'fontsize',16);
%ax = gca;
%ax.XLim = [0 10000];
%ax.YLim = [7700 10800];

rapd = temp2(1:500:end,:);
figure();
% style = ['bo','go','ro','co','mo','yo','ko','wo','b.','g.'];
plot(X_axis,rapd(:,1),'ro',X_axis,rapd(:,2),'rx',X_axis,rapd(:,3),'r.-',X_axis,rapd(:,4),'r*',X_axis,rapd(:,5),'ks',X_axis,rapd(:,6),'kd',X_axis,rapd(:,7),'kv',X_axis,rapd(:,8),'k^',X_axis,rapd(:,9),'kp',X_axis,rapd(:,10),'k+','MarkerSize',7);
legend('Non stationary type-1','Non stationary type-2','New stationary type-1','New stationary type-2','Clique Whittle Policy','Graphical Whittle Policy','Cella algorithm','slotted ALOHA','MWS algorithm','Lyapunov drift algorithm','NumColumns',4,'Location','southeast','fontsize',16);
ylabel("Packets Dropped",'fontweight','bold','fontsize',16);
xlabel("Time",'fontweight','bold','fontsize',16);
title(arrival_size+" arrivals, "+transmission+" case",'fontsize',16);
%ax = gca;
%ax.XLim = [9995.7 10000];
%ax.YLim = [17 40];

% figure('Position', get(0, 'Screensize'));
% X_axis = 1:sim_time;
% plot(X_axis,cot(:,1),'bo',X_axis,cot(:,1),'gx',X_axis,cot(:,2),'r+',X_axis,cot(:,3),'c*',X_axis,cot(:,4),'ms',X_axis,cot(:,5),'yd',X_axis,cot_neelay,'kv');
% legend('Clique Whittle Policy','Graphical Whittle Policy','New stationary Whittle Policy','Cella algorithm','slotted ALOHA','MWS algorithm','Lyapunov drift algorithm','Location','southeast','fontsize',16);
% ylabel("Average Cost",'fontweight','bold','fontsize',16);
% xlabel("Time",'fontweight','bold','fontsize',16);
% title(arrival_size+" arrivals, "+transmission+" case",'fontsize',16);
% saveas(gcf,"../../results/new_stationary_1/AC"+"_"+num2str(arrival_size)+"_"+num2str(transmission)+".png");
% 
% figure('Position', get(0, 'Screensize'));
% X_axis = 1:sim_time;
% plot(X_axis,rapd_wc,'bo',X_axis,rapd_wg,'gx',X_axis,rapd_nsw,'r+',X_axis,rapd_cella,'c*',X_axis,rapd_sa,'ms',X_axis,rapd_mwis,'yd',X_axis,rapd_neelay,'kv');
% legend('Clique Whittle Policy','Graphical Whittle Policy','New stationary Whittle Policy','Cella algorithm','slotted ALOHA','MWS algorithm','Lyapunov drift algorithm','Location','southeast','fontsize',16);
% ylabel("Packets Dropped",'fontweight','bold','fontsize',16);
% xlabel("Time",'fontweight','bold','fontsize',16);
% title(arrival_size+" arrivals, "+transmission+" case",'fontsize',16);
% saveas(gcf,"../../results/new_stationary_1/PD"+"_"+num2str(arrival_size)+"_"+num2str(transmission)+".png");
% %save("matfiles/"+arrival_size+"_"+transmission+".mat");