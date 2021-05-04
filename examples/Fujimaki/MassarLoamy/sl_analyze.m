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
fig_pos.bottom=0.8;
fig_pos.length=0.3;
fig_pos.height=0.15;

%nt=10;
%a.fig=figure;
a.fs=12;
a.lw=2;
a.cz=8;
  fs = 5; % sampling frequency
  % Creating the movie : quality = 100%, no compression is used and the
  % timing is adjusted to the sampling frequency of the time interval
  qt=100;
  %A number from 0 through 100. Higher quality numbers result in higher video quality and larger file sizes
  a.fig=figure;
  %set(gcf,'Units','normalized', 'WindowStyle','docked','OuterPosition',[0 0 1 1]);  % maximize the plotting figure
  mov =  VideoWriter('linux.avi');% avifile('pvc1.avi','quality',qt,'compression','indeo5','fps',fs);
  mov.FrameRate = 5;mov.Quality=qt;
  open(mov);

for nt=2:1:length(nod)-1

    fprintf('Plotting the %i out of %i results\n', nt, length(nod));

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
    
    
    %% --------- plot sat profiles  ------------------
    a.sub2=subplot('position'...
         ,[fig_pos.left,...
         fig_pos.bottom-3*fig_pos.height-0.25,...
          fig_pos.length-0.09,fig_pos.height+0.1]);
    plot(nod(nt).terms{s_idx}(1:inp.nn1)...
         ,nod(nt).terms{y_idx}(1:inp.nn1)...
         ,'linewidth',a.lw) ;hold on
    plot(slab(2,:,1),slab(1,:,1),'rd',slab(2,:,2),...
        slab(1,:,2),'go',slab(2,:,3),slab(1,:,3),'cx','linewidth',a.lw);hold off;
    axis([-0.05, 1.05,0,0.05])
    get(gca,'xtick');
    set(gca,'fontsize',12);
    xlabel('sat (-)','FontSize',a.fs);
    ylabel('y (m)','FontSize',a.fs);
    
    %% --------- plot conc profiles  ------------------
    a.sub2=subplot('position'...
         ,[fig_pos.left+fig_pos.length,...
         fig_pos.bottom-3*fig_pos.height-0.25,...
          fig_pos.length-0.09,fig_pos.height+0.1]);
    plot(nod(nt).terms{c_idx}(1:inp.nn1)...
         ,nod(nt).terms{y_idx}(1:inp.nn1)...
         ,'linewidth',a.lw) ;hold on
    plot(clab(2,:,1),clab(1,:,1),'rd',clab(2,:,2),clab(1,:,2)...
    ,'go',clab(2,:,3),clab(1,:,3),'cx','MarkerSize',a.cz,'linewidth',a.lw);hold off
    axis([-0.05, 0.28,0,0.05])
    get(gca,'xtick');
    set(gca,'fontsize',12);
    xlabel('conc (-)','FontSize',a.fs);
    ylabel('y (m)','FontSize',a.fs);
    
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
    
    F = getframe(gcf); % save the current figure
    writeVideo(mov,F);% add it as the next frame of the movie
end

close(mov);
