et1_kgs=-arrayfun(@(y) y.qin(1),bcof);
area1=0.5D-6; % hard coding
et1_mmday=et1_kgs/area1*86400;
time_day=[bcof.tout]/3600/24;

% saturation profile
%d(1).terms{p_idx}(1:inp.nn1)

%dp0_idx=min(0-nod(1).terms{2}(1:inp.nn1)); 

% write x and y coordinates in matrix form.
x_matrix=reshape(nod(1).terms{x_idx},[inp.nn1,inp.nn2]);
y_matrix=reshape(nod(1).terms{y_idx},[inp.nn1,inp.nn2]);


dp0_idx=inp.nn1;
dp1_idx=inp.nn1-5;

p_top=arrayfun(@(y) y.terms{p_idx}(inp.nn1),nod);
sw_top=arrayfun(@(y) y.terms{s_idx}(inp.nn1),nod);


fig_pos.left=0.09;
fig_pos.bottom=0.6;
fig_pos.length=0.25;
fig_pos.height=0.3;
fig_pos.hori_gap=0.05

%nt=10;


% this is copied form Saltmarsh/Program/Paper1/Fujimaki/Toyoura/5.0
load('Toyoura_paper1.mat');
a.fs=12;
a.lw=2;
a.cz=8;
fs = 5; % sampling frequency
a.fig=figure;
%% --------- plot sat profiles  ------------------
a.sub2=subplot('position'...
     ,[fig_pos.left...
     ,fig_pos.bottom...
     ,fig_pos.length,fig_pos.height]);
iv1=2;iv2=3;
plot(a1(5,(f3(5)+2):iv1:(2*(f3(5)+1)),j(2)),cl-a1(2,(f3(5)+2):iv1:(2*(f3(5)+1)),j(2))*scl,'rv-','LineWidth',a.lw);hold on
plot(a1(5,(f3(5)+2):iv2:(2*(f3(5)+1)),j(3)),cl-a3(2,(f3(5)+2):iv2:(2*(f3(5)+1)),j(3))*scl,'go-','LineWidth',a.lw);hold on
plot(slab(2,1:8,1),cl-slab(1,1:8,1)*scl,'kd','MarkerFaceColor','k','MarkerSize',cz);hold on;
set(gca,'FontSize',15,'FontWeight','bold','linewidth',a.lw)
title('(a)','fontweight','b','Units', 'normalized','Position', [0 1],'HorizontalAlignment', 'left','FontSize',12)
xlabel('Saturation (-)','FontSize',a.fs)
ylabel('Depth (mm)','fontsize',a.fs)
ax1 = gca;
set(ax1,'YAxisLocation','Left','LineWidth',lw,'YDir','reverse','XAxisLocation','bottom')
axis([-0.05 0.65 a1(2,(f3(5)+2),j(1)) a1(2,(2*(f3(5)+1)),j(1))*scl])

%% --------- plot conc profiles  ------------------
a.sub2=subplot('position'...
     ,[fig_pos.left+fig_pos.hori_gap+fig_pos.length...
     ,fig_pos.bottom...
     ,fig_pos.length,fig_pos.height]);
plot(nod(end).terms{c_idx}(1:inp.nn1)...
     ,50-nod(end).terms{y_idx}(1:inp.nn1)*1000 ...
     ,'rv-','linewidth',a.lw) ;hold on
plot(a1(4,[(f3(5)+2):af(3):(2*(f3(5)-int)),(2*(f3(5)-int+1)):(2*(f3(5)+1))],j(3))...
	,50-a1(2,[(f3(5)+2):af(3):(2*(f3(5)-int)),(2*(f3(5)-int+1)):(2*(f3(5)+1))],j(3))*scl,'go-','LineWidth',a.lw);hold on 
%plot(slab(2,:,1),slab(1,:,1),'rd',slab(2,:,2),...
%    slab(1,:,2),'go',slab(2,:,3),slab(1,:,3),'cx','linewidth',a.lw);hold on;
plot(clab(2,:,3),50-clab(1,:,3)*1000,'kd','linewidth',a.lw,'markerfacecolor','k');hold off;
ax1=gca;
title('(b)','fontweight','b','Units', 'normalized','Position', [0 1],'HorizontalAlignment', 'left','FontSize',12)
xlabel('Concentration (kg kg^{-1})','FontSize',fs)
set(ax1,'YTickLabel','','LineWidth',lw,'YDir','reverse','XAxisLocation','bottom')

%% --------- plot y velocity  ------------------
a.sub2=subplot('position'...
     ,[fig_pos.left+2*fig_pos.length,...
     fig_pos.bottom-3*fig_pos.height-0.25,...
      fig_pos.length-0.09,fig_pos.height+0.1]);
plot(ele(nt).terms{vy_idx}(1:inp.nn1-1)*86400 ...
     ,ele(nt).terms{yele_idx}(1:inp.nn1-1)...
     ,'linewidth',a.lw) 
%axis([-0.05, 0.28,0,0.05])
get(gca,'xtick');
set(gca,'fontsize',12);
xlabel('vy (mm/day)','FontSize',a.fs);
ylabel('y (m)','FontSize',a.fs);










for nt=2:10:length(nod)-1
%% -------------  sub 1 ET over time  --------------
a.sub1=subplot('position'...
     ,[fig_pos.left,fig_pos.bottom,...
      fig_pos.length,fig_pos.height]);



a.plot1=plot(time_day(1:nt),et1_mmday(1:nt),...
         'k-','linewidth',a.lw);hold on
a.plot1=plot(eslab(1,:),eslab(2,:),'cx','linewidth',a.lw);
get(gca,'xtick');
set(gca,'fontsize',12);
xlabel('Time (day)','FontSize',a.fs);
ylabel('Evt (mm/day)','FontSize',a.fs);
axis([-1 time_day(end) -1 21])
%title('ead profile');
%legend('show','Location','East')



%% -------------  sub 2 sat over time  --------------
a.sub2=subplot('position'...
     ,[fig_pos.left+0.5,fig_pos.bottom,...
      fig_pos.length,fig_pos.height]);
a.plot2=plot(time_day(1:nt),...
       arrayfun(@(y) y.terms{s_idx}(inp.nn1),nod(1:nt)),...
         'r-','linewidth',a.lw);hold on
a.plot2=plot(time_day(1:nt),...
       arrayfun(@(y) y.terms{s_idx}(inp.nn1-2),nod(1:nt)),...
         'g-','linewidth',a.lw);hold on
a.plot2=plot(time_day(1:nt),...
       arrayfun(@(y) y.terms{s_idx}(inp.nn1-10),nod(1:nt)),...
         'b-','linewidth',a.lw);hold on
get(gca,'xtick');
set(gca,'fontsize',12);
xlabel('Time (day)','FontSize',a.fs);
ylabel('sat (-)','FontSize',a.fs);
axis([-1 time_day(end) -0.1 1.1])


%% -------------  sub 3 conc over time  --------------
a.sub2=subplot('position'...
     ,[fig_pos.left,fig_pos.bottom-fig_pos.height,...
      fig_pos.length,fig_pos.height]);
a.plot2=plot(time_day(1:nt),...
       arrayfun(@(y) y.terms{c_idx}(inp.nn1),nod(1:nt)),...
         'r-','linewidth',a.lw);hold on
a.plot2=plot(time_day(1:nt),...
       arrayfun(@(y) y.terms{c_idx}(inp.nn1-2),nod(1:nt)),...
         'g-','linewidth',a.lw);hold on
a.plot2=plot(time_day(1:nt),...
       arrayfun(@(y) y.terms{c_idx}(inp.nn1-10),nod(1:nt)),...
         'b-','linewidth',a.lw);hold on
%a.plot2=plot(time_day,sw_top,...
%         'k-','linewidth',a.lw);hold on
get(gca,'xtick');
set(gca,'fontsize',12);
xlabel('Time (day)','FontSize',a.fs);
ylabel('Conc','FontSize',a.fs);
axis([-1 time_day(end) -0.005 0.3])


%% -------------  sub 4 sm over time  --------------

salt_dep=arrayfun(@(y) y.terms{sm_idx}(inp.nn1),nod(1:nt))/2165/area1*1000;
a.sub2=subplot('position'...
     ,[fig_pos.left+0.5,fig_pos.bottom-fig_pos.height,...
      fig_pos.length,fig_pos.height]);
a.plot2=plot(time_day(1:nt),...
       salt_dep, ...  
         'r-','linewidth',a.lw);hold on
get(gca,'xtick');
set(gca,'fontsize',12);
xlabel('Time (day)','FontSize',a.fs);
ylabel('solid (kg)','FontSize',a.fs);
%axis([-1 time_day(end) -0.005 0.3])


%% -------- contour plot on saturation ---------
a.sub2=subplot('position'...
     ,[fig_pos.left,fig_pos.bottom-2*fig_pos.height-0.07,...
      fig_pos.length+0.02,fig_pos.height]);
% write pressure and conc in matrix form.
s_matrix  = reshape(nod(nt).terms{s_idx},[inp.nn1,inp.nn2]);

contourf(x_matrix,y_matrix,s_matrix);
colorbar
get(gca,'xtick');
set(gca,'fontsize',12);
xlabel('x (m)','FontSize',a.fs);
ylabel('y (m)','FontSize',a.fs);


%% -------- contour plot on concentration ---------
a.sub2=subplot('position'...
     ,[fig_pos.left+0.5,...
     fig_pos.bottom-2*fig_pos.height-0.07,...
      fig_pos.length+0.02,fig_pos.height]);

% write pressure and conc in matrix form.
c_matrix = reshape(nod(nt).terms{c_idx},[inp.nn1,inp.nn2]);

contourf(x_matrix,y_matrix,c_matrix);
colorbar
get(gca,'xtick');
set(gca,'fontsize',12);
xlabel('x (m)','FontSize',a.fs);
ylabel('y (m)','FontSize',a.fs);
%axis([10, 40,9,10])



F = getframe(gcf); % save the current figure
writeVideo(mov,F);% add it as the next frame of the movie
end
savefig(a.fig,'Calibration.fig')
close(mov);
