function plotWavenumberLimits(k_min, k_max)

hold on
plot3(k_min*cos(linspace(0,2*pi,100)) , k_min*sin(linspace(0,2*pi,100)), ones(100,1),'k','linewidth',4)
plot3(k_max*cos(linspace(0,2*pi,100)) , k_max*sin(linspace(0,2*pi,100)), ones(100,1),'k','linewidth',4)
plot3(k_min*cos(linspace(0,2*pi,100)) , k_min*sin(linspace(0,2*pi,100)), ones(100,1),'m','linewidth',2)
plot3(k_max*cos(linspace(0,2*pi,100)) , k_max*sin(linspace(0,2*pi,100)), ones(100,1),'m','linewidth',2)
xlim(1.5*k_max*[-1 1 ])
ylim(1.5*k_max*[-1 1 ])
hold off