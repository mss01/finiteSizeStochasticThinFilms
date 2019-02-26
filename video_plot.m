function M = video_plot(L_flat, x,Y, t_new)

area(x',Y)
ylim([0 1.4])
xlim([-1.5*L_flat 1.5*L_flat])
xlabel('x [-]','Fontsize',16)
ylabel('h [-]','Fontsize',16)
set(gca,'FontSize',18)
set(gca,'YTick',[0 0.2 0.4 0.6 0.8 1.0 1.2 1.4])
legend(t_new)
M = print('-RGBImage',sprintf('-r%d',150)); %% 150 = resolution, normally gives a good animation, and also occupies less space
end