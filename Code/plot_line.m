function plot_line(x,I,sty,lim)

plot(x,I,sty,'linewidth',2);
xlabel('x (\mum)','FontSize',16);ylabel('Normalised intensity','FontSize',16);
xlim([-lim*10^-6 lim*10^-6]);
set(gca, 'FontSize', 16);
set(gca,'XTick',-lim*10^-6:(lim/2)*10^-6:lim*10^-6); 
set(gca,'XTickLabel',{-lim,-lim/2,0,lim/2,lim});
set(gca,'YTick',0:0.2:1); 
set(gca,'YTickLabel',{'0','0.2','0.4','0.6','0.8','1'});
legend('z = 0','z = z_R','z = 2z_R');

end