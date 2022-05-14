clc;
D=load('Varying_connectivity_unlimited_new.mat');
X_axis = 1:10000;
ax1 = subplot(2,1,1);
plot(X_axis(100:400:end),D.cot_w_g4(100:400:end),'o',X_axis(100:400:end),D.cot_w_g6(100:400:end),'+',X_axis(100:400:end),D.cot_w_g8(100:400:end),'*',X_axis(100:400:end),D.cot_w_g10(100:400:end),'^',X_axis(100:400:end),D.cot_w_g12(100:400:end),'s',X_axis(100:400:end),D.cot_w_g15(100:400:end),'d')
hLegend=legend('0.4','0.6','0.8','1','1.2','1.5');
hLegend.Title.String = 'd_{thresh-computaion}';
xlabel('Time');
ylabel('Average Cost');
%title('M100 fixedG small arrival unlimited');
ax2 = subplot(2,1,2);
plot(X_axis(100:400:end),D.rapd_w_g4(100:400:end),'o',X_axis(100:400:end),D.rapd_w_g6(100:400:end),'+',X_axis(100:400:end),D.rapd_w_g8(100:400:end),'*',X_axis(100:400:end),D.rapd_w_g10(100:400:end),'^',X_axis(100:400:end),D.rapd_w_g12(100:400:end),'s',X_axis(100:400:end),D.rapd_w_g15(100:400:end),'d')
SLegend=legend('0.4','0.6','0.8','1','1.2','1.5');
SLegend.Title.String = 'd_{thresh-computaion}';
%legend('Packets dropped whittle','Packets dropped graphical','Packets dropped SA','Packets dropped MWIS');
xlabel('Time');
ylabel('Packets Dropped');