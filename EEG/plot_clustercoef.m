%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % plots cluster coefficient, from unweighted graph
  % analysis using phase lag index
  
  % VBarone Nov, 2022

function h = plot_clustercoef(ch_list,toplot,pres)


    %channels location
    load('Standard_10-20_81ch.mat', 'locations');
    ch_list = upper(ch_list);
    idx = NaN(length(ch_list),1);
        for ch = 1:length(ch_list)
            if isempty(find(strcmp(locations.labels,ch_list{ch})))
                warning('[plot_topography] Cannot find the %s electrode.',ch_list{ch});
                ch_list{ch} = [];
                idx(ch)     = [];
            else
                idx(ch) = find(strcmp(locations.labels,ch_list{ch}));
            end
        end
     % Creating the rectangle grid (-1,1)
    [ch_x, ch_y] = pol2cart((pi/180).*((-1).*locations.theta(idx)+90), ...
                                locations.radius(idx));     % X, Y channel coords
    INTERP_POINTS = 1000;
    % Points out of the head to reach more natural interpolation
    r_ext_points = 1.2;
    [add_x, add_y] = pol2cart(0:pi/4:7*pi/4,r_ext_points*ones(1,8));
    linear_grid = linspace(-r_ext_points,r_ext_points,INTERP_POINTS);         % Linear grid (-1,1)
    [interp_x, interp_y] = meshgrid(linear_grid, linear_grid);
    
    
% Global parameters
    %   Note: Head radius should be set as 0.6388 if the 10-20 system is used.
    %   This number was calculated taking into account that the distance from Fpz
    %   to Oz is d=2*0.511. Thus, if the circle head must cross the nasion and
    %   the inion, it should be set at 5d/8 = 0.6388.
    %   Note2: When the number of interpolation points rises, the plots become
    %   smoother and more accurate, however, computational time also rises.
    HEAD_RADIUS     = 5*2*0.6388/8;  % 1/2  of the nasion-inion distance
    HEAD_EXTRA      = 1*2*0.6388/8;  % 1/10 of the nasion-inion distance
    k = 4;        


% Plotting the head limits as a circle         
    head_rho    = HEAD_RADIUS;                      % Head radius
    head_theta  = linspace(0,2*pi,INTERP_POINTS);   % From 0 to 360????
    head_x      = head_rho.*cos(head_theta);        % Cartesian X of the head
    head_y      = head_rho.*sin(head_theta);        % Cartesian Y of the head
    plot(head_x, head_y, 'Color', 'k', 'LineWidth',4);
    hold on;

    % Plotting the nose
    nt = 0.15;      % Half-nose width (in percentage of pi/2)
    nr = 0.22;      % Nose length (in radius units)
    nose_rho   = [head_rho, head_rho+head_rho*nr, head_rho];
    nose_theta = [(pi/2)+(nt*pi/2), pi/2, (pi/2)-(nt*pi/2)];
    nose_x     = nose_rho.*cos(nose_theta);
    nose_y     = nose_rho.*sin(nose_theta);
    plot(nose_x, nose_y, 'Color', 'k', 'LineWidth',4);
    hold on;

    % Plotting the ears as ellipses
    ellipse_a = 0.08;                               % Horizontal exentricity
    ellipse_b = 0.16;                               % Vertical exentricity
    ear_angle = 0.9*pi/8;                           % Mask angle
    offset    = 0.05*HEAD_RADIUS;                   % Ear offset
    ear_rho   = @(ear_theta) 1./(sqrt(((cos(ear_theta).^2)./(ellipse_a^2)) ...
        +((sin(ear_theta).^2)./(ellipse_b^2))));    % Ellipse formula in polar coords
    ear_theta_right = linspace(-pi/2-ear_angle,pi/2+ear_angle,INTERP_POINTS);
    ear_theta_left  = linspace(pi/2-ear_angle,3*pi/2+ear_angle,INTERP_POINTS);
    ear_x_right = ear_rho(ear_theta_right).*cos(ear_theta_right);          
    ear_y_right = ear_rho(ear_theta_right).*sin(ear_theta_right); 
    ear_x_left  = ear_rho(ear_theta_left).*cos(ear_theta_left);         
    ear_y_left  = ear_rho(ear_theta_left).*sin(ear_theta_left); 
    plot(ear_x_right+head_rho+offset, ear_y_right, 'Color', 'k', 'LineWidth',4); hold on;
    plot(ear_x_left-head_rho-offset, ear_y_left, 'Color', 'k', 'LineWidth',4); hold on;
    hold on
    scatter(ch_x, ch_y, 60,'k', 'LineWidth',1.5);
    hold on
    
    %plot connections between nodes
    for l = 1:19
      for m = 1:19
        if toplot(l,m) ==1
           plot([ch_x(l) ch_x(m)],[ch_y(l) ch_y(m)], 'k-', 'LineWidth', 1)
           hold on   
        end
      end  
    end
    axis off
    if contains(pres, 'pres')
       title('Preserved', 'Units', 'normalized', 'Position', [0.5, 1.04, 0],'FontName','times','FontSize',14)
    elseif contains(pres, 'unp')  
       title('Unpreserved', 'Units', 'normalized', 'Position', [0.5, 1.04, 0],'FontName','times','FontSize',14)

    end   

hold off
end