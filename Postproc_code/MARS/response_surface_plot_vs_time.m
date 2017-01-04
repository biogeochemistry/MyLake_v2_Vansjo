
figure (2)

x = [1:length(post_big_run(1,:))];
plot_id = 1;

for i = 1:20
    
    
    subplot(4,5,i)
    
    plot(x,post_big_run(plot_id:plot_id+19,:))
    plot_id =  plot_id+20;
end

