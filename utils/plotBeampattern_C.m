function [P P_max] =  plotBeampattern_C(pos, k_min, k_max)
% function plotBeampattern_C(pos, k_min, k_max)
% plots the modulo of the spectrum P(k,phi)
% P the complex spectrum for given wavenumbers and direction of arrival (DOA)
% K_vec vector of wavenumbers
% Theta_vec vector of azimuths [rad]

sideview=false;

K_vec=linspace(k_min, k_max, 200);
Theta_vec=linspace(0, 2*pi, 360);

N_sensors=size(pos,1);


Kx=transpose(K_vec)*cos(Theta_vec);
Ky=transpose(K_vec)*sin(Theta_vec);

P=zeros(size(Kx));
for nn=1:N_sensors
  for mm=nn+1:N_sensors
	  P=P+cos( 2*pi*( Kx*(pos(nn,1)-pos(mm,1)) + Ky*(pos(nn,2)-pos(mm,2)) ) );
  end
end
P=2*P/(N_sensors^2)+1/N_sensors;

P_max = max(P(:));


K_vec=linspace(0, k_max*1.5, 200);
Theta_vec=linspace(0, 2*pi, 360);

N_sensors=size(pos,1);


Kx=transpose(K_vec)*cos(Theta_vec);
Ky=transpose(K_vec)*sin(Theta_vec);

P=zeros(size(Kx));
for nn=1:N_sensors
  for mm=nn+1:N_sensors
	  P=P+cos( 2*pi*( Kx*(pos(nn,1)-pos(mm,1)) + Ky*(pos(nn,2)-pos(mm,2)) ) );
  end
end
P=2*P/(N_sensors^2)+1/N_sensors;




[r th]=meshgrid(K_vec,Theta_vec);
[x y]=pol2cart(th,r);

if(~isreal(P))
    P=abs(P);
end

z=transpose(P);

if(sideview)
    subplot(1,2,1)
end
surf(x,y,z+0.1) % we add this small value (+0.1) for the purpose of visualization. Otherwise values =0 are rendered poorly.

axis square
view(2)
xlim([min(x(:)) max(x(:))])
ylim([min(y(:)) max(y(:))])
shading flat
shading interp
% set(gcf,'Renderer','Zbuffer');
set(gcf, 'Renderer', 'opengl');

colormap(jet)
% colorbar B/W
% cmap = [linspace(1,0,256).^2';];
% colormap([cmap cmap cmap])
% set(gcf,'Renderer','Zbuffer')

if(sideview)
    subplot(1,2,2)
    surf(x,y,z)
    axis square
    view(90,0)
    xlim([min(x(:)) max(x(:))])
    ylim([min(y(:)) max(y(:))])
%     colorbar
    shading flat
    shading interp
    set(gcf, 'Renderer', 'opengl');
%     set(gcf,'Renderer','Zbuffer')
end
xlabel('Wavenumber x [1/m]')
ylabel('Wavenumber y [1/m]')
return
% see more code here
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/243095
