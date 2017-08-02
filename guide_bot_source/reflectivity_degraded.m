% showcasing degraded reflectivity curve
clear all;close all;

set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 24)

h=figure
R0=0.99;
qc=0.0217;
alpha=6.07;
m=5;
W=0.003;

q=0:0.001:0.13;

for i=1:length(q)
    if q(i)<qc
        R1(i)=R0;
    else
        R1(i)=0.5*R0.*(1-tanh((q(i)-m.*qc)./W)).*(1-alpha.*(q(i)-qc));
        
        %R(i)=0;
    end
end

R0=0.99;
qc=0.0217;
alpha=6.07*1.4;
m=5*0.8;
W=0.003;

for i=1:length(q)
    if q(i)<qc
        R2(i)=R0;
    else
        R2(i)=0.5*R0.*(1-tanh((q(i)-m.*qc)./W)).*(1-alpha.*(q(i)-qc));
        
        %R(i)=0;
    end
end


plot(q/qc,R1,'k','linewidth',2)
hold on
plot(q/qc,R2,'r','linewidth',2)
set(gca,'XTick',[0 1 2 3 4 5 6])

xlabel('$q$ $[q_c]$','interpreter','latex')
ylabel('Reflectivity','interpreter','latex')
%title('Reflectivity curve comparison')

print(h,'-dpng','-r250',['reflectivity_curve.png'])
print(h,'-depsc','reflectivity_degraded.eps')

%%

% showcasing degraded reflectivity curve
clear all;close all;

set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 24)

h=figure
set(h, 'Position', [0 0 464 400])
axes('fontsize',24)
box on
R0=0.99;
qc=0.0217;
alpha=6.07;%4.35
m=2:7;
W=0.003;

q0=0;
q1=0.2;
q=linspace(0,q1,1000);

for j=1:length(m)
for i=1:length(q)
    if q(i)<qc
        R1(i,j)=R0;
    else
        R1(i,j)=0.5*R0.*(1-tanh((q(i)-m(j).*qc)./W)).*(1-alpha.*(q(i)-qc));
    end
end
end
color = {'k' 'b' 'r' 'm' 'g' [0.0 0.4 0.0]};
%set(0,'DefaultTextfontsize',40);
hold on
for i = length(m):-1:1
handle(i) = plot(q,R1(:,i),'color',color{i},'linewidth',2)
end
ylim([0 1.1])
xlim([0 0.22])
set(gca,'linewidth',2)


L = legend(handle,'\quad $m=2$\quad','\quad $m=3$ ','\quad $m=4$ ','\quad $m=5$ ','\quad $m=6$ ','\quad $m=7$ ','false');
set(L,'fontsize',20,'interpreter','latex')
set(L,'linewidth',1)


%set(L,'PlotBoxAspectRatioMode','manual');
%set(L,'PlotBoxAspectRatio',[1 1 1]);
%set(L,'outerposition',[0.6455 0.5047 0.2572 0.5060]);

% hc = get(L,'children')
% 
% yd = zeros(length(hc)/3,1);
% for n = 1:length(yd);
%   yd(n) = get(hc(1+(n-1)*3),'YData');
% end
% 
% yd2 = yd.*1.2-0.2;
% 
% % markers
% for n = 1:length(yd2)
%    set(hc(1+(n-1)*3),'YData',yd2(n));
% end
% 
% % lines
% for n = 1:length(yd2)
%    set(hc(2+(n-1)*3),'YData',[yd2(n),yd2(n)]);
% end
% 
% % text
% for n = 1:length(yd2);
%     pos = get(hc(3+(n-1)*3),'Position');
%     pos(2) = yd2(n);
%     set(hc(3+(n-1)*3),'Position',pos);
% end

%set(L,'PlotBoxAspectRatioMode','manual');
%set(L,'PlotBoxAspectRatio',[1 1 1]);
xlabel('$q$ [\AA$^{-1}$]','fontsize',26,'interpreter','latex')
ylabel('Reflectivity','fontsize',26,'interpreter','latex')

print(h,'-dpng','-r350',['reflectivity_curve.png'])
print(h,'-depsc','reflectivity_curves2.eps')
