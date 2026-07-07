function [g2]=plot_mundi(r,env,kk,swi,app_,rsun,sat,v,attitude,cases)
%kk is the number of the plot
%swi= if you want to define the ecplipse period, use swi=1

%figure(kk) %ijk

if isempty(app_)
    app_ = axes();
end

C1 = repmat([0,0,0],numel(r(r(:,5)==1,1)),1);
C2 = repmat([0,0.70,1],numel(r(r(:,5)==0,1)),1);
S1 = repmat(5,numel(r(r(:,5)==1,1)),1);
S2 = repmat(5,numel(r(r(:,5)==0,1)),1);

switch swi
    case 0
        g2=plot3(app_,r(:,1),r(:,2),r(:,3),'.','color',[0 0.8 1]);
    case 1
%         g2.ecl=plot3(r(r(:,5)==1,1),r(r(:,5)==1,2),r(r(:,5)==1,3),'.','color',[0 0.8 1]);
%         g2.sun=plot3(r(r(:,5)==0,1),r(r(:,5)==0,2),r(r(:,5)==0,3),'.','color','k');
        g2.ecl=scatter3(app_,r(r(:,5)==1,1),r(r(:,5)==1,2),r(r(:,5)==1,3),S1,C1,'MarkerFaceColor','flat');
        hold(app_,'on')
        g2.sun=scatter3(app_,r(r(:,5)==0,1),r(r(:,5)==0,2),r(r(:,5)==0,3),S2,C2,'MarkerFaceColor','flat');
end
%plot Planet
%space_color = 'k';
npanels = 180;   % Number of globe panels around the equator deg/panel = 360/npanels
alpha   = 1; % globe transparency level, 1 = opaque, through 0 = invisible
%GMST0 = []; % Don't set up rotatable globe (ECEF)
GMST0 = 4.89496121282306; % Set up a rotatable globe at J2000.0

% Earth texture image
% Anything imread() will handle, but needs to be a 2:1 unprojected globe
% image.

image_file = 'Earth.jpg';

% Mean spherical earth

erad    = env.Rp; % equatorial radius (meters)
prad    = env.Rp; % polar radius (meters)
erot    = 7.2921158553e-5; % earth rotation rate (radians/sec)

%% Create figure

set(app_, 'Visible','off');
axis(app_,'equal');

% Set initial view
view(app_,0,30);
axis(app_,'vis3d');

%% Create wireframe globe

% Create a 3D meshgrid of the sphere points using the ellipsoid function

[x, y, z] = ellipsoid(0, 0, 0, erad, erad, prad, npanels);
hold(app_,'on')
globe = surf(app_,x, y, -z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);

if ~isempty(GMST0)
    hgx = hgtransform('Parent',app_);
    set(hgx,'Matrix', makehgtform('zrotate',GMST0));
    set(globe,'Parent',hgx);
end

%% Texturemap the globe

% Load Earth image for texture map

cdata = imread(image_file);

% Set image as color data (cdata) property, and set face color to indicate
% a texturemap, which Matlab expects to be in cdata. Turn off the mesh edges.

set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');
%sun vector
sunvect=[zeros(3,1) [rsun(1,1)/norm(rsun)*(env.Rp+5000) ...
     rsun(1,2)/norm(rsun)*(env.Rp+5000) rsun(1,3)/norm(rsun)*(env.Rp+5000)]'];         %Vector to plot sun vector
sun.pl=plot3(app_,sunvect(1,:),sunvect(2,:),sunvect(3,:),'-','linewidth',3,'color',[1 0.65 0]);
x_ax=[zeros(3,1) [1*(env.Rp+5000) 0 0]'];         %x_vector
y_ax=[zeros(3,1) [0 1*(env.Rp+5000) 0]'];         %y_vector
z_ax=[zeros(3,1) [0 0 1*(env.Rp+5000)]'];         %z_vector
xyz(1).pl=plot3(app_,x_ax(1,:),x_ax(2,:),x_ax(3,:),'-','linewidth',3,'color',[1 0 0]);
xyz(2).pl=plot3(app_,y_ax(1,:),y_ax(2,:),y_ax(3,:),'-','linewidth',3,'color',[0 1 0]);
xyz(3).pl=plot3(app_,z_ax(1,:),z_ax(2,:),z_ax(3,:),'-','linewidth',3,'color',[0 0 1]);
legend(app_,[xyz(1).pl,xyz(2).pl,xyz(3).pl,sun.pl],{'X ECI','Y ECI','Z ECI','Sun Vector'});

%% Satellite orientation

if exist('sat','var')
    n_pos = 20;
    x_ax_sc=[zeros(3,1) [1*1000 0 0]'];         %x_vector
    y_ax_sc=[zeros(3,1) [0 1*1000 0]'];         %y_vector
    z_ax_sc=[zeros(3,1) [0 0 1*1000]'];         %z_vector
    switch cases
        case 1
            tex_at='info1';
            tex_sp='info';
        case 2
            tex_at='info1_cold';
            tex_sp='info_cold';
    end
    if isfield(attitude,tex_at) % if attitude contains the field info1 or info1_cold -> Nadir-Velocity pointing
        cont_color = 0;
        for j=round(1:(size(r,1)/n_pos):size(r,1)) % for every temporal instant...
            cont_color = cont_color + 1;
            mat=Eci2body(sat,attitude,r(j,1:3),v(j,1:3),cases); %  BODY -> ECI
            x_sc = mat'*x_ax_sc + r(j,1:3)';
            y_sc = mat'*y_ax_sc + r(j,1:3)';
            z_sc = mat'*z_ax_sc + r(j,1:3)';
            xyz(4).pl = plot3(app_,x_sc(1,:),x_sc(2,:),x_sc(3,:),'-','linewidth',1.5,'color',[1 0 0]);
            xyz(5).pl = plot3(app_,y_sc(1,:),y_sc(2,:),y_sc(3,:),'-','linewidth',1.5,'color',[0 1 0]);
            xyz(6).pl = plot3(app_,z_sc(1,:),z_sc(2,:),z_sc(3,:),'-','linewidth',1.5,'color',[0 0 1]);
            if cont_color == 1
                xyz(7).pl = plot3(app_,r(j,1),r(j,2),r(j,3),'o','color',[0.5 1 1],'MarkerSize',5,'LineWidth',5);
            end
        end
    elseif isfield(attitude,tex_sp)  % if attitude contains the field info or info_cold -> random spin
        a0=[0 0 0];                                                            %initial angle
        rota=0.3*rand(size(r,1),3);                                         %angular velocity random (maximum 5deg/s)    
        a=zeros(size(r,1),3);                                               %rotation angle for each direction
        a_=rota.*r(:,4);
        cont_color = 0;
        for j=round(1:(size(r,1)/n_pos):size(r,1)) % for every temporal instant...
            cont_color = cont_color + 1;
            a(j,:)=a0+a_(j,:);
            a0=a(j,:);
            mat=Eci2body2(r,a(j,:)); %  ECI -> BODY
            x_sc = mat*x_ax_sc + r(j,1:3)';
            y_sc = mat*y_ax_sc + r(j,1:3)';
            z_sc = mat*z_ax_sc + r(j,1:3)';
            xyz(4).pl = plot3(app_,x_sc(1,:),x_sc(2,:),x_sc(3,:),'-','linewidth',0.5,'color',[1 0 0]);
            xyz(5).pl = plot3(app_,y_sc(1,:),y_sc(2,:),y_sc(3,:),'-','linewidth',0.5,'color',[0 1 0]);
            xyz(6).pl = plot3(app_,z_sc(1,:),z_sc(2,:),z_sc(3,:),'-','linewidth',0.5,'color',[0 0 1]);
            if cont_color == 1
                xyz(7).pl = plot3(app_,r(j,1),r(j,2),r(j,3),'o','color',[0.5 1 1],'MarkerSize',5,'LineWidth',5);
            end
        end
    end
    legend(app_,[xyz(1).pl,xyz(2).pl,xyz(3).pl,sun.pl,xyz(4).pl,xyz(5).pl,xyz(6).pl,xyz(7).pl],{'X ECI','Y ECI','Z ECI','Sun Vector','X S/C','Y S/C','Z S/C','Time = 0 s'});
end
hold(app_,'off');
end