mex reedsconnect.c
mex reedsdistance.c
clear, cla
delta=5*pi/180;
n=1000;
L=0.68;
r=L/tan(0.4);
hold on;
hi=quiver(0,0,1,1,'color','b','LineWidth',1);
he=quiver(0,0,1,1,'color','r','LineWidth',1);
reeds = robotics.ReedsSheppConnection('MinTurningRadius',r);
t1=tic;
for i = 1:n
    
    xi=[10*rand(1),10*rand(1),2*pi*rand(1)];
    xe=[10*rand(1),10*rand(1),2*pi*rand(1)];
    set(hi,'XData',xi(1),'YData',xi(2),'UData',L*cos(xi(3)),'VData',L*sin(xi(3)));
    set(he,'XData',xe(1),'YData',xe(2),'UData',L*cos(xe(3)),'VData',L*sin(xe(3)));
    h1=plot(path(:,1),path(:,2),'r-');
    quiver(x(1),x(2),L/3*cos(x(3)),L/3*sin(x(3)),'color','b','LineWidth',2);
    
%    [pose,l(i)] = connect(reeds,xi,xe);	%% MATLAB's own implementation
%    path = interpolate(pose{1},[0:0.3:l(i)]);  %% works only with licensed robotics toolbox with 							%% R2018b

    [l(i),t,u,v,num] = reedsdistance(xi(1),xi(2),xi(3),xe(1),xe(2),xe(3),r);
    [x,y,theta]=reedsconnect(t,u,v,num,xi(1),xi(2),xi(3),r,delta,1000);
    h=plot(x,y,'k--');
     axis equal
    drawnow update 
    pause(0.1)
    delete(h)
end
toc(t1)
