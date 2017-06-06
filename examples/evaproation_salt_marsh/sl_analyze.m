%% a check about when the smass and others becomes negative
smass_neg=struct;
i=0;
for n=1:length(nod)
  neg=find( nod(n).terms{smass_idx}<0   );
  if ~isempty(neg)
    i=i+1;
    smass_neg(i).nodeno=neg;
    smass_neg(i).outputno=n;
  end
end

sm_neg=struct;
i=0;
for n=1:length(nod)
  neg=find( nod(n).terms{sm_idx}<0);
  if ~isempty(neg)
    i=i+1;
    sm_neg(i).nodeno=neg;
    sm_neg(i).outputno=n;
  end
end



wmass_neg=struct;
i=0;
for n=1:length(nod)
  neg=find( nod(n).terms{wmass_idx}<0);
  if ~isempty(neg)
    i=i+1;
    wmass_neg(i).nodeno=neg;
    wmass_neg(i).outputno=n;
  end
end

x_nod_mtx=reshape(nod(1).terms{xnod_idx},[inp.nn1,inp.nn2])';
y_nod_mtx=reshape(nod(1).terms{ynod_idx},[inp.nn1,inp.nn2])';

x_ele_mtx=reshape(ele(1).terms{xele_idx},[inp.nn1-1,inp.nn2-1])';
yele_mtx=reshape(ele(1).terms{yele_idx},[inp.nn1-1,inp.nn2-1])';

et_x_ay=nod(1).terms{xnod_idx}(bcof(1).i);
et_y_ay=nod(1).terms{ynod_idx}(bcof(1).i);

% nod
%nodeindx ind_surface_ay=nod(1).terms(

%% now we analyze the problem
a.iv=12;
a.ivy=2;
a.fs=15;

a.left =  0.08;
a.bot  =  0.80;
a.width= 0.9;
a.width2=0.823;
a.height= 0.15;
a.h_interval=0.15;
a.lw=2;
%
%
n=8000;

fs = 5; % sampling frequency
% Creating the movie : quality = 100%, no compression is used and the
% timing is adjusted to the sampling frequency of the time interval
qt=100;
%A number from 0 through 100. Higher quality numbers result in higher video quality and larger file sizes
a.fig=figure;
set(gcf,'Units','normalized', 'WindowStyle','docked','OuterPosition',[0 0 1 1]);  % maximize the plotting figure
mov =  VideoWriter('linux.avi');% avifile('pvc1.avi','quality',qt,'compression','indeo5','fps',fs);
mov.FrameRate = 5;mov.Quality=qt;
open(mov);

for n=400:1:length(nod)-1
  s_mtx=reshape(nod(n).terms{s_nod_idx},[inp.nn1,inp.nn2])';
  c_mtx=reshape(nod(n).terms{c_nod_idx},[inp.nn1,inp.nn2])';
  p_mtx=reshape(nod(n).terms{p_nod_idx},[inp.nn1,inp.nn2])';
  vx_mtx=reshape(ele(n+1).terms{vx_ele_idx},[inp.nn1-1,inp.nn2-1])';
  vy_mtx=reshape(ele(n+1).terms{vy_ele_idx},[inp.nn1-1,inp.nn2-1])';
  
   
  % mask matrix for water table
  p_mtx_watertable_mask=zeros(size(p_mtx))+1;
  p_mtx_watertable_mask(p_mtx<0)=0;
  % http://stackoverflow.com/questions/15716898/matlab-find-row-indice-of-first-occurrence-for-each-column-of-matrix-without-u
  [~,idx] = max(p_mtx_watertable_mask(:,sum(p_mtx_watertable_mask)>0));
  %http://stackoverflow.com/questions/32941314/matlab-get-data-from-a-matrix-with-data-row-and-column-indeces-stored-in-arrays/32941538#32941538
  p_mtx_watertable_mask = sub2ind(size(p_mtx_watertable_mask),idx,[1:size(p_mtx,2)]);
  watertable_ay=p_mtx(p_mtx_watertable_mask)/9800+y_nod_mtx(p_mtx_watertable_mask);
  
  
  
  n_fig=0;
  %% sat over x axis 
  subplot('position',[a.left,a.bot-n_fig*a.h_interval,a.width2,a.height])
  plot(et_x_ay,s_mtx(1,:),'k-','linewidth',a.lw)
  get(gca,'xtick');
  set(gca,'fontsize',15);
  xlabel('Time (day)','FontSize',a.fs);
  ylab=ylabel('Saturation','FontSize',a.fs);
  set(ylab, 'Units', 'Normalized', 'Position', [-0.06, 0.5, 0]);
  ax1 = gca;
  set(ax1,'XTickLabel','')
  axis([et_x_ay(1) et_x_ay(end) -0.4 1.1])
  title(sprintf('Simulation %f/%f  days',bcop(n).tout/3600/24,bcop(end).tout/3600/24))
  
  n_fig=n_fig+1;
  %%  plot saturation over time
  subplot('position',[a.left,a.bot-n_fig*a.h_interval,a.width,a.height])
  contourf ( x_nod_mtx,y_nod_mtx,s_mtx);hold on
  quiver ( x_ele_mtx(1:a.iv:end,1:a.ivy:end),yele_mtx(1:a.iv:end,1:a.ivy:end)...
        ,vx_mtx(1:a.iv:end,1:a.ivy:end),vy_mtx(1:a.iv:end,1:a.ivy:end)...
	,'w-','linewidth',a.lw);hold on
  plot(et_x_ay,watertable_ay,'g-','linewidth',a.lw)
  plot(et_x_ay,watertable_ay,'g-','linewidth',a.lw)
  colorbar
  set(gca,'fontsize',15);
  xlabel('','FontSize',a.fs);
  ylab=ylabel('Y (m)','FontSize',a.fs);
  set(ylab, 'Units', 'Normalized', 'Position', [-0.06, 0.5, 0]);
  %set(gca,'XTick',[])
  ax1 = gca;
  set(ax1,'XTickLabel','')
  axis([et_x_ay(1) et_x_ay(end) 0.0 8.5])
  
  n_fig=n_fig+1;
  %% concentration over x 
  subplot('position',[a.left,a.bot-n_fig*a.h_interval,a.width2,a.height])
  plot(et_x_ay,c_mtx(1,:),'k-','linewidth',a.lw)
  get(gca,'xtick');
  set(gca,'fontsize',15);
  xlabel('');
  ylab=ylabel('C (kg/kg)','FontSize',a.fs);
  set(ylab, 'Units', 'Normalized', 'Position', [-0.06, 0.5, 0]);
  ax1 = gca;
  set(ax1,'XTickLabel','')
  axis([et_x_ay(1) et_x_ay(end) -0.05 0.39])
  
  n_fig=n_fig+1;
  %% plot concentration contour
  subplot('position',[a.left,a.bot-n_fig*a.h_interval,a.width,a.height])
  contourf ( x_nod_mtx,y_nod_mtx,c_mtx);hold on
  colorbar
  plot(et_x_ay,watertable_ay,'g-','linewidth',a.lw)
  set(gca,'fontsize',15);
  xlabel('');
  %set(gca,'XTick',[])
  ax1 = gca;
  set(ax1,'XTickLabel','')
  ylab=ylabel('C (kg/kg)','FontSize',a.fs);
  set(ylab, 'Units', 'Normalized', 'Position', [-0.06, 0.5, 0]);
  axis([et_x_ay(1) et_x_ay(end) 0.0 8.5])


  n_fig=n_fig+1;
  %% ET over TIME
  subplot('position',[a.left,a.bot-n_fig*a.h_interval,a.width2,a.height])
  et_mmday=-bcof(n).qin./inp.qinc'*3600*24;
  set(gca,'fontsize',15);
  set(gca,'XTick',[])
  plot(et_x_ay,et_mmday,'k-','linewidth',a.lw);hold off
  %xlabel('');
  ax1 = gca;
  set(ax1,'XTickLabel','')
  ylab=ylabel('Evp.(mm/day)','FontSize',a.fs);
  set(ylab, 'Units', 'Normalized', 'Position', [-0.06, 0.5, 0]);
  axis([et_x_ay(1) et_x_ay(end) -3 18])
  
  n_fig=n_fig+1;
  %% plot head over time
  subplot('position',[a.left,a.bot-n_fig*a.h_interval,a.width2,a.height])
  exchange_flux_mmday=bcop(n).qpl./inp.qinc'*3600*24;
  plot(et_x_ay,exchange_flux_mmday,'k-','linewidth',a.lw)
  set(gca,'fontsize',15);
  xlabel('X (m)');
  ylab=ylabel('flux (mm/day)','FontSize',a.fs);
  set(ylab, 'Units', 'Normalized', 'Position', [-0.06, 0.5, 0]);
  axis([et_x_ay(1) et_x_ay(end) -90 1050])
  
  F = getframe(gcf); % save the current figure
  writeVideo(mov,F);% add it as the next frame of the movie
end
close(mov);
  
  
  
  
  
  
saveas(a.fig,'result.fig')  ;
  
  
%   [et_x_ay et_y_ay bcop(n).pbc bcop(n).pvec bcop(n).qpl bcof(n).qin]
%   tide(bcop(n).itout,3)
