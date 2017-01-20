function graphs( A,x,m,t,years,t_resolution,filename)
% this part only for better visualization of graphs %

	fig_counter = 0;
	for i=0:size(A,1)-1
		[specie_vis,t_vis]    = visual_graphs(A{i+1},t,years,m,t_resolution);

		if rem(i,4) == 0
			fig_counter=fig_counter + 1;
		    h(fig_counter) = figure('units','normalized','outerposition',[0 0 1 1]);
		end
		s = subplot(2,2,rem(i,4)+1); mesh(t_vis,x,specie_vis); view(160, 15);xlabel('Time (Years)');ylabel('Depth(cm)');zlabel(A(i+1,2));
		set(gca,'xdir','reverse')
		title(s, A(i+1,2))
	end	

	hgsave(h,strcat(filename,'.fig'))
end

function [ C_vis, t_vis ] = visual_graphs(C,t,years_of_simulations,time_steps,yearly_resolution)
%VISUAL_GRAPHS this is function for better representing of the graphs
% 
% This function removes overloaded data from concentration matrices.
% For instance, as result of our computation we may have the matrix of 
% concentrations with size 511x100000 and this is very hard to represent 
% graphically. This function removes overloaded data and return the concentration 
% with needed resolution, for instance 36 lines per 1 year

lines_of_waterfall = years_of_simulations*yearly_resolution; %resolution of the visualization%
skip_lines = floor(time_steps/lines_of_waterfall);           %

for i=2:time_steps
  if rem(i,skip_lines) == 0
    C_vis(:,i/skip_lines) = C(:,i);
    t_vis(:,i/skip_lines) = t(:,i);
  end
end

end
