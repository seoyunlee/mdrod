%% Molecular Dynamics Simulator
% Seoyun Soy Lee

% 2/26/15
% torque implemented and no longer self imploding
% forces are additive

%% Experimental Parameters
N = 9; % number of rods
R = 1; % radius of circles/distance away from line
L = 4; % length of center line
K = 10; % spring constant
M = 3; % mass of particles
I=M*(L/2)^2; % moment of inertia

Lx=10*2;  % size of box
Ly=10*2;  

%% Initial Conditions
[x, y]=ndgrid(L/2+R:L+2*R:Lx-L/2-R,L/2+R:L+2*R:Lx-L/2-R);  
[~, ii]=sort(rand(1,numel(x)));
x=x(ii(1:N)); % random center x
y=y(ii(1:N)); % random center y
theta=2*pi*rand(1,N); % random angle
vx=randn(1,N)/10;
vy=randn(1,N)/10;
w=zeros(1,N); % beginning angular velocity

ax_old=0*x;
ay_old=0*y;
alpha_old=0*x;

total=zeros(1,100000);

%% Display Parameters
plotit=true;
Nplotskip=400;

dt=0.25e-1;

%% Setup Plotting
clf;
h=zeros(1,2*N);
hl=zeros(1,2*N);
for np=1:N
  h(np)=rectangle('Position',[x(np)+cos(theta(np))*L/2-R y(np)+sin(theta(np))*L/2-R 2*R 2*R],'Curvature',[1 1],'edgecolor','b');
  h(np+N)=rectangle('Position',[x(np)-cos(theta(np))*L/2-R y(np)-sin(theta(np))*L/2-R 2*R 2*R],'Curvature',[1 1],'edgecolor','b');
  hl(np)=line([x(np)+L/2*cos(theta(np))+R*cos(theta(np)-pi/2),x(np)-L/2*cos(theta(np))+R*cos(theta(np)-pi/2)],...
              [y(np)+L/2*sin(theta(np))+R*sin(theta(np)-pi/2),y(np)-L/2*sin(theta(np))+R*sin(theta(np)-pi/2)],'Color','b');
  hl(np+N)=line([x(np)+L/2*cos(theta(np))-R*cos(theta(np)-pi/2),x(np)-L/2*cos(theta(np))-R*cos(theta(np)-pi/2)],...
                [y(np)+L/2*sin(theta(np))-R*sin(theta(np)-pi/2),y(np)-L/2*sin(theta(np))-R*sin(theta(np)-pi/2)],'Color','b');
end
axis('equal');
axis([0 Lx 0 Ly]);

%% Main Loop

for t=1:100000
    
    if(rem(t,Nplotskip)==0)
        for np=1:N
            %plot circles
            set(h(np),'Position',[x(np)+cos(theta(np))*L/2-R, y(np)+sin(theta(np))*L/2-R, 2*R, 2*R]);
            set(h(np+N),'Position',[x(np)-cos(theta(np))*L/2-R, y(np)-sin(theta(np))*L/2-R, 2*R, 2*R]);
            set(hl(np),'XData',[x(np)+L/2*cos(theta(np))+R*cos(theta(np)-pi/2),x(np)-L/2*cos(theta(np))+R*cos(theta(np)-pi/2)],...
                       'YData',[y(np)+L/2*sin(theta(np))+R*sin(theta(np)-pi/2),y(np)-L/2*sin(theta(np))+R*sin(theta(np)-pi/2)]);
            set(hl(np+N),'XData',[x(np)+L/2*cos(theta(np))-R*cos(theta(np)-pi/2),x(np)-L/2*cos(theta(np))-R*cos(theta(np)-pi/2)],...
                         'YData',[y(np)+L/2*sin(theta(np))-R*sin(theta(np)-pi/2),y(np)-L/2*sin(theta(np))-R*sin(theta(np)-pi/2)]);
        end
        drawnow;
    end

    % KE and PE calcs
    PE=0;
    
    x=x+vx*dt+ax_old.*dt.^2/2;  % first step in Verlet integration
    y=y+vy*dt+ay_old.*dt.^2/2;
    theta=theta+w*dt+alpha_old.*dt.^2/2;
    
    % Interaction detector and Force Law
    Fx=zeros(1,N);
    Fy=zeros(1,N);
    torque=zeros(1,N);
    
    x1=x+cos(theta).*L/2;
    y1=y+sin(theta).*L/2;
    x2=x-cos(theta).*L/2;
    y2=y-sin(theta).*L/2;
    
    % check to see which distance is shortest between each particle
    for nn=1:N-1
        for mm=nn+1:N
            % End of rod mm hits somewhere on rod nn
            dot=(x1(mm)-x(nn))*(x1(nn)-x(nn))+(y1(mm)-y(nn))*(y1(nn)-y(nn));
            nnm=dot/(L/2);
            if(nnm < L/2 && nnm > -L/2)
                tx=-(-(x1(mm)-x(nn))+nnm*cos(theta(nn)));
                ty=-(-(y1(mm)-y(nn))+nnm*sin(theta(nn)));
                if(sqrt(tx.^2+ty.^2)<2*R)
                    F=-K*(2*R/sqrt(tx.^2+ty.^2)-1);
                    Fx(nn)=Fx(nn)+F.*tx;  % particle-particle Force Law
                    Fx(mm)=Fx(mm)-F.*tx;
                    Fy(nn)=Fy(nn)+F.*ty;  % particle-particle Force Law
                    Fy(mm)=Fy(mm)-F.*ty;
                    torque(nn)=torque(nn)+F*nnm*(cos(theta(nn))*ty-sin(theta(nn))*tx);
                    torque(mm)=torque(mm)-F*L/2*(cos(theta(mm))*ty-sin(theta(mm))*tx);
                    PE=PE+K/2*((2*R-sqrt(tx.^2+ty.^2))^2);
                end
            end
            
            dot=(x2(mm)-x(nn))*(x1(nn)-x(nn))+(y2(mm)-y(nn))*(y1(nn)-y(nn));
            nnm=dot/(L/2);
            if(nnm < L/2 && nnm > -L/2)
                tx=-(-(x2(mm)-x(nn))+nnm*cos(theta(nn)));
                ty=-(-(y2(mm)-y(nn))+nnm*sin(theta(nn)));
                if(sqrt(tx.^2+ty.^2)<2*R)
                    F=-K*(2*R/sqrt(tx.^2+ty.^2)-1);
                    Fx(nn)=Fx(nn)+F.*tx;  % particle-particle Force Law
                    Fx(mm)=Fx(mm)-F.*tx;
                    Fy(nn)=Fy(nn)+F.*ty;  % particle-particle Force Law
                    Fy(mm)=Fy(mm)-F.*ty;
                    torque(nn)=torque(nn)+F*nnm*(cos(theta(nn))*ty-sin(theta(nn))*tx);
                    torque(mm)=torque(mm)-F*(-L/2)*(cos(theta(mm))*ty-sin(theta(mm))*tx);
                    PE=PE+K/2*((2*R-sqrt(tx.^2+ty.^2))^2);
                end
            end
            
            dot=(x1(nn)-x(mm))*(x1(mm)-x(mm))+(y1(nn)-y(mm))*(y1(mm)-y(mm));
            mmn=dot/(L/2);
            if(mmn < L/2 && mmn > -L/2)
                tx=-((x1(nn)-x(mm))-mmn*cos(theta(mm)));
                ty=-((y1(nn)-y(mm))-mmn*sin(theta(mm)));
                if(sqrt(tx.^2+ty.^2)<2*R)
                    F=-K*(2*R/sqrt(tx.^2+ty.^2)-1);
                    Fx(nn)=Fx(nn)+F.*tx;  % particle-particle Force Law
                    Fx(mm)=Fx(mm)-F.*tx;
                    Fy(nn)=Fy(nn)+F.*ty;  % particle-particle Force Law
                    Fy(mm)=Fy(mm)-F.*ty;
                    torque(nn)=torque(nn)+F*L/2*(cos(theta(nn))*ty-sin(theta(nn))*tx);
                    torque(mm)=torque(mm)-F*mmn*(cos(theta(mm))*ty-sin(theta(mm))*tx);
                    PE=PE+K/2*((2*R-sqrt(tx.^2+ty.^2))^2);
                end    
            end
            
            dot=(x2(nn)-x(mm))*(x1(mm)-x(mm))+(y2(nn)-y(mm))*(y1(mm)-y(mm));
            mmn=dot/(L/2);
            if(mmn < L/2 && mmn > -L/2)
                tx=-((x2(nn)-x(mm))-mmn*cos(theta(mm)));
                ty=-((y2(nn)-y(mm))-mmn*sin(theta(mm)));
                if(sqrt(tx.^2+ty.^2)<2*R)
                    F=-K*(2*R/sqrt(tx.^2+ty.^2)-1);
                    Fx(nn)=Fx(nn)+F.*tx;  % particle-particle Force Law
                    Fx(mm)=Fx(mm)-F.*tx;
                    Fy(nn)=Fy(nn)+F.*ty;  % particle-particle Force Law
                    Fy(mm)=Fy(mm)-F.*ty;
                    torque(nn)=torque(nn)+F*(-L/2)*(cos(theta(nn))*ty-sin(theta(nn))*tx);
                    torque(mm)=torque(mm)-F*mmn*(cos(theta(mm))*ty-sin(theta(mm))*tx);
                    PE=PE+K/2*((2*R-sqrt(tx.^2+ty.^2))^2);
                end
            end

            tx=x1(mm)-x1(nn);
            ty=y1(mm)-y1(nn);
            if(sqrt(tx.^2+ty.^2)<2*R)
                F=-K*(2*R/sqrt(tx.^2+ty.^2)-1);
                Fx(nn)=Fx(nn)+F.*tx;  % particle-particle Force Law
                Fx(mm)=Fx(mm)-F.*tx;
                Fy(nn)=Fy(nn)+F.*ty;  % particle-particle Force Law
                Fy(mm)=Fy(mm)-F.*ty;
                torque(nn)=torque(nn)+F*L/2*(cos(theta(nn))*ty-sin(theta(nn))*tx);
                torque(mm)=torque(mm)-F*L/2*(cos(theta(mm))*ty-sin(theta(mm))*tx);
                PE=PE+K/2*((2*R-sqrt(tx.^2+ty.^2))^2);
            end
            tx=x1(mm)-x2(nn);
            ty=y1(mm)-y2(nn);
            if(sqrt(tx.^2+ty.^2)<2*R)
                F=-K*(2*R/sqrt(tx.^2+ty.^2)-1);
                Fx(nn)=Fx(nn)+F.*tx;  % particle-particle Force Law
                Fx(mm)=Fx(mm)-F.*tx;
                Fy(nn)=Fy(nn)+F.*ty;  % particle-particle Force Law
                Fy(mm)=Fy(mm)-F.*ty;
                torque(nn)=torque(nn)+F*-L/2*(cos(theta(nn))*ty-sin(theta(nn))*tx);
                torque(mm)=torque(mm)-F*L/2*(cos(theta(mm))*ty-sin(theta(mm))*tx);
                PE=PE+K/2*((2*R-sqrt(tx.^2+ty.^2))^2);
            end

            tx=x2(mm)-x2(nn);
            ty=y2(mm)-y2(nn);
            if(sqrt(tx.^2+ty.^2)<2*R)
                F=-K*(2*R/sqrt(tx.^2+ty.^2)-1);
                Fx(nn)=Fx(nn)+F.*tx;  % particle-particle Force Law
                Fx(mm)=Fx(mm)-F.*tx;
                Fy(nn)=Fy(nn)+F.*ty;  % particle-particle Force Law
                Fy(mm)=Fy(mm)-F.*ty;
                torque(nn)=torque(nn)+F*-L/2*(cos(theta(nn))*ty-sin(theta(nn))*tx);
                torque(mm)=torque(mm)-F*-L/2*(cos(theta(mm))*ty-sin(theta(mm))*tx);
                PE=PE+K/2*((2*R-sqrt(tx.^2+ty.^2))^2);
            end

            tx=x2(mm)-x1(nn);
            ty=y2(mm)-y1(nn);
            if(sqrt(tx.^2+ty.^2)<2*R)
                F=-K*(2*R/sqrt(tx.^2+ty.^2)-1);
                Fx(nn)=Fx(nn)+F.*tx;  % particle-particle Force Law
                Fx(mm)=Fx(mm)-F.*tx;
                Fy(nn)=Fy(nn)+F.*ty;  % particle-particle Force Law
                Fy(mm)=Fy(mm)-F.*ty;
                torque(nn)=torque(nn)+F*L/2*(cos(theta(nn))*ty-sin(theta(nn))*tx);
                torque(mm)=torque(mm)-F*-L/2*(cos(theta(mm))*ty-sin(theta(mm))*tx);
                PE=PE+K/2*((2*R-sqrt(tx.^2+ty.^2))^2);
            end
        end
    end
    
    PEw=0;
    
    Fx=Fx-K*(x1-R).*(x1<R);  % Left wall
    PEw=PEw+K/2*sum(((x1-R).*(x1<R)).^2);
    Fx=Fx-K*(x2-R).*(x2<R);
    PEw=PEw+K/2*sum(((x2-R).*(x2<R)).^2);
    
    Fy=Fy-K*(y1-R).*(y1<R);  % Bottom wall
    PEw=PEw+K/2*sum(((y1-R).*(y1<R)).^2);
    Fy=Fy-K*(y2-R).*(y2<R);
    PEw=PEw+K/2*sum(((y2-R).*(y2<R)).^2);
    
    Fx=Fx-K*(x1-(Lx-R)).*(x1>Lx-R);  % Right wall
    PEw=PEw+K/2*sum(((x1-(Lx-R)).*(x1>Lx-R)).^2);
    Fx=Fx-K*(x2-(Lx-R)).*(x2>Lx-R);
    PEw=PEw+K/2*sum(((x2-(Lx-R)).*(x2>Lx-R)).^2);
    
    Fy=Fy-K*(y1-(Ly-R)).*(y1>Ly-R);
    PEw=PEw+K/2*sum(((y1-(Ly-R)).*(y1>Ly-R)).^2);
    Fy=Fy-K*(y2-(Ly-R)).*(y2>Ly-R); % Top wall
    PEw=PEw+K/2*sum(((y2-(Ly-R)).*(y2>Ly-R)).^2);
    
    for i=1:N
        % Left wall
        torque(i)=torque(i)+L/2*sin(theta(i))*K*(x1(i)-R).*(x1(i)<R);
        torque(i)=torque(i)+(-L/2)*sin(theta(i))*K*(x2(i)-R).*(x2(i)<R);
        % Bottom wall
        torque(i)=torque(i)-L/2*cos(theta(i))*K*(y1(i)-R).*(y1(i)<R);
        torque(i)=torque(i)-(-L/2)*cos(theta(i))*K*(y2(i)-R).*(y2(i)<R);
        % Right wall
        torque(i)=torque(i)+L/2*sin(theta(i))*K*(x1(i)-(Lx-R)).*(x1(i)>Lx-R);
        torque(i)=torque(i)+(-L/2)*sin(theta(i))*K*(x2(i)-(Lx-R)).*(x2(i)>Lx-R);
        % Top wall
        torque(i)=torque(i)-L/2*cos(theta(i))*K*(y1(i)-(Ly-R)).*(y1(i)>(Ly-R));
        torque(i)=torque(i)-(-L/2)*cos(theta(i))*K*(y2(i)-(Ly-R)).*(y2(i)>(Ly-R));
    end
    
    ax=Fx./M;
    ay=Fy./M;
    alpha=torque./I;
    
    vx=vx+(ax_old+ax).*dt/2;  % second step in Verlet integration
    vy=vy+(ay_old+ay).*dt/2;
    w=w+(alpha_old+alpha).*dt/2;
    
    
    ax_old=ax;
    ay_old=ay;
    alpha_old=alpha;
    
     % translational KE
    KE=M/2*sum(vx.^2+vy.^2);

    % rotational KE
    KR=I/2*sum(w.^2);
    
    total(t)=KE+KR+PE+PEw;
    
    if(t>1 && total(t)-total(t-1)>1e-4)
        clf;
        h=zeros(1,2*N);
        hl=zeros(1,2*N);
        for np=1:N
            h(np)=rectangle('Position',[x(np)+cos(theta(np))*L/2-R y(np)+sin(theta(np))*L/2-R 2*R 2*R],'Curvature',[1 1],'edgecolor','b');
            h(np+N)=rectangle('Position',[x(np)-cos(theta(np))*L/2-R y(np)-sin(theta(np))*L/2-R 2*R 2*R],'Curvature',[1 1],'edgecolor','b');
            hl(np)=line([x(np)+L/2*cos(theta(np))+R*cos(theta(np)-pi/2),x(np)-L/2*cos(theta(np))+R*cos(theta(np)-pi/2)],...
                [y(np)+L/2*sin(theta(np))+R*sin(theta(np)-pi/2),y(np)-L/2*sin(theta(np))+R*sin(theta(np)-pi/2)],'Color','b');
            hl(np+N)=line([x(np)+L/2*cos(theta(np))-R*cos(theta(np)-pi/2),x(np)-L/2*cos(theta(np))-R*cos(theta(np)-pi/2)],...
                [y(np)+L/2*sin(theta(np))-R*sin(theta(np)-pi/2),y(np)-L/2*sin(theta(np))-R*sin(theta(np)-pi/2)],'Color','b');
        end
        axis('equal');
        axis([0 Lx 0 Ly]);
        pause;
    end


end    

plot(total);