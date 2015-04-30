%% Molecular Dynamics Simulator

% 4/29/15
% implement for rods of two different lengths/radii
% 4/29/15
% working model for the one contact method
% all rods in this are one length/radius

%% Experimental Parameters
N1 = 5; % number of rods of size 1
N2 = 4; % number of rods of size 2
N = N1+N2; % total number of rods
R1 = 1; % radius of circles/distance away from line
R2 = 0.5;
R = [repmat(R1,1,N1), repmat(R2,1,N2)];
ratio = R/(R1+R2);
L1 = 4; % length of center line
L2 = 3;
L = [repmat(L1,1,N1), repmat(L2,1,N2)];
g1=L1/R1;
g2=L2/R2;
g = [repmat(g1,1,N1), repmat(g2,1,N2)];
K = 10; % spring constant
M = 3; % mass of particles
I = M*(L/2).^2;
B=1; % drag

Lx=20;  % size of box
Ly=20; 

% for the packing
Pthresh=1e-5;
Pstep=1e-4;
Ptol=1e-6;

% Different size particles
% place random points and find the minimum distance and set that equal to R
%% Initial Conditions
Rmax = max([R1 R2]);
Lmax = max([L1 L2]);
[x, y]=ndgrid(Lmax/2+Rmax:Lmax+2*Rmax:Lx-Lmax/2-Rmax,L/2+Rmax:Lmax+2*Rmax:Lx-Lmax/2-Rmax);  
[~, ii]=sort(rand(1,numel(x)));
x=x(ii(1:N)); % random center x
y=y(ii(1:N)); % random center y
theta=2*pi*rand(1,N); % random angle
vx=randn(1,N)/30;
vy=randn(1,N)/30;
w=zeros(1,N); % beginning angular velocity

ax_old=0*x;
ay_old=0*y;
alpha_old=0*x;

%% Simulation Parameters
dt=0.25e-1; % time step
Nt=10000; % max time

%% Display Parameters
plotit=true;
Nplotskip=400;

%% Save energy calculations
total=zeros(1,Nt);
PEs=zeros(1,Nt);
KEs=zeros(1,Nt);
KRs=zeros(1,Nt);
PEws=zeros(1,Nt);

%% Setup Plotting
clf;
h=zeros(1,2*N);
hl=zeros(1,2*N);
% size 1 particle
for np=1:N
  h(np)=rectangle('Position',[x(np)+cos(theta(np))*L(np)/2-R(np) y(np)+sin(theta(np))*L(np)/2-R(np) 2*R(np) 2*R(np)],'Curvature',[1 1],'edgecolor','b');
  h(np+N)=rectangle('Position',[x(np)-cos(theta(np))*L(np)/2-R(np) y(np)-sin(theta(np))*L(np)/2-R(np) 2*R(np) 2*R(np)],'Curvature',[1 1],'edgecolor','b');
  hl(np)=line([x(np)+L(np)/2*cos(theta(np))+R(np)*cos(theta(np)-pi/2),x(np)-L(np)/2*cos(theta(np))+R(np)*cos(theta(np)-pi/2)],...
              [y(np)+L(np)/2*sin(theta(np))+R(np)*sin(theta(np)-pi/2),y(np)-L(np)/2*sin(theta(np))+R(np)*sin(theta(np)-pi/2)],'Color','b');
  hl(np+N)=line([x(np)+L(np)/2*cos(theta(np))-R(np)*cos(theta(np)-pi/2),x(np)-L(np)/2*cos(theta(np))-R(np)*cos(theta(np)-pi/2)],...
                [y(np)+L(np)/2*sin(theta(np))-R(np)*sin(theta(np)-pi/2),y(np)-L(np)/2*sin(theta(np))-R(np)*sin(theta(np)-pi/2)],'Color','b');
end
axis('equal');
axis([0 Lx 0 Ly]);

%% Set Up Movie
% writerObj=VideoWriter('simrods.avi');
% open(writerObj);
% 
% axis('equal');
% axis([0 Lx 0 Ly]);
% set(gca,'nextplot','replacechildren');
% set(gcf,'Renderer','zbuffer');

%% Main Loop
Dh=-1;
Dl=-1;

for nt=1:Nt
    
    if(rem(nt,Nplotskip)==0)
        for np=1:N
            %plot circles
            set(h(np),'Position',[x(np)+cos(theta(np))*L(np)/2-R(np), y(np)+sin(theta(np))*L(np)/2-R(np), 2*R(np), 2*R(np)]);
            set(h(np+N),'Position',[x(np)-cos(theta(np))*L(np)/2-R(np), y(np)-sin(theta(np))*L(np)/2-R(np), 2*R(np), 2*R(np)]);
            set(hl(np),'XData',[x(np)+L(np)/2*cos(theta(np))+R(np)*cos(theta(np)-pi/2),x(np)-L(np)/2*cos(theta(np))+R(np)*cos(theta(np)-pi/2)],...
                       'YData',[y(np)+L(np)/2*sin(theta(np))+R(np)*sin(theta(np)-pi/2),y(np)-L(np)/2*sin(theta(np))+R(np)*sin(theta(np)-pi/2)]);
            set(hl(np+N),'XData',[x(np)+L(np)/2*cos(theta(np))-R(np)*cos(theta(np)-pi/2),x(np)-L(np)/2*cos(theta(np))-R(np)*cos(theta(np)-pi/2)],...
                         'YData',[y(np)+L(np)/2*sin(theta(np))-R(np)*sin(theta(np)-pi/2),y(np)-L(np)/2*sin(theta(np))-R(np)*sin(theta(np)-pi/2)]);
        end
        drawnow;
%         frame=getframe;
%         writeVideo(writerObj,frame);
    end

    % KE and PE calcs
    PE=0;
    
    x=x+vx*dt+ax_old.*dt.^2/2;  % first step in Verlet integration
    y=y+vy*dt+ay_old.*dt.^2/2;
    theta=theta+w*dt+alpha_old*dt^2/2;
    
    % Interaction detector and Force Law
    Fx=zeros(1,N);
    Fy=zeros(1,N);
    torque=zeros(1,N);
    
    x1=x+cos(theta).*L/2;
    y1=y+sin(theta).*L/2;
    x2=x-cos(theta).*L/2;
    y2=y-sin(theta).*L/2;
    
    C=0; % counts
    % check to see which distance is shortest between each particle
    for nn=1:N-1
        for mm=nn+1:N
            Dnm=(R(nn)+R(mm));
            dnmfinal=Dnm;
            txfinal=0;
            tyfinal=0;
            nnarm=0;
            mmarm=0;
            
            % End of rod mm hits somewhere on rod nn
            dot=(x1(mm)-x(nn))*(x1(nn)-x(nn))+(y1(mm)-y(nn))*(y1(nn)-y(nn));
            proj=dot/(L(nn)/2);
            if(proj < L(nn)/2 && proj > -L(nn)/2)
                tx=x1(mm)-x(nn)-proj*cos(theta(nn));
                ty=y1(mm)-y(nn)-proj*sin(theta(nn));
                dnm=tx^2+ty^2;
                if(dnm<Dnm^2 && dnm<dnmfinal^2)
                    dnmfinal=sqrt(dnm);
                    txfinal=tx;
                    tyfinal=ty;
                    nnarm=proj;
                    mmarm=L(mm)/2;
                end
            end
            
            dot=(x2(mm)-x(nn))*(x1(nn)-x(nn))+(y2(mm)-y(nn))*(y1(nn)-y(nn));
            proj=dot/(L(nn)/2);
            if(proj < L(nn)/2 && proj > -L(nn)/2)
                tx=x2(mm)-x(nn)-proj*cos(theta(nn));
                ty=y2(mm)-y(nn)-proj*sin(theta(nn));
                dnm=tx^2+ty^2;
                if(dnm<Dnm^2 && dnm<dnmfinal^2)
                    dnmfinal=sqrt(dnm);
                    txfinal=tx;
                    tyfinal=ty;
                    nnarm=proj;
                    mmarm=-L(mm)/2;
                end
            end
            
            dot=(x1(nn)-x(mm))*(x1(mm)-x(mm))+(y1(nn)-y(mm))*(y1(mm)-y(mm));
            proj=dot/(L(mm)/2);
            if(proj < L(mm)/2 && proj > -L(mm)/2)
                tx=-((x1(nn)-x(mm))-proj*cos(theta(mm)));
                ty=-((y1(nn)-y(mm))-proj*sin(theta(mm)));
                dnm=tx^2+ty^2;
                if(dnm<Dnm^2 && dnm<dnmfinal^2)
                    dnmfinal=sqrt(dnm);
                    txfinal=tx;
                    tyfinal=ty;
                    nnarm=L(nn)/2;
                    mmarm=proj;
                end
            end
            
            dot=(x2(nn)-x(mm))*(x1(mm)-x(mm))+(y2(nn)-y(mm))*(y1(mm)-y(mm));
            proj=dot/(L(mm)/2);
            if(proj < L(mm)/2 && proj > -L(mm)/2)
                tx=-((x2(nn)-x(mm))-proj*cos(theta(mm)));
                ty=-((y2(nn)-y(mm))-proj*sin(theta(mm)));
                dnm=tx^2+ty^2;
                if(dnm<Dnm^2 && dnm<dnmfinal^2)
                    dnmfinal=sqrt(dnm);
                    txfinal=tx;
                    tyfinal=ty;
                    nnarm=-L(nn)/2;
                    mmarm=proj;
                end
            end
            
            tx=x1(mm)-x1(nn);
            ty=y1(mm)-y1(nn);
            dnm=tx^2+ty^2;
            if(dnm<Dnm^2 && dnm<dnmfinal^2)
                dnmfinal=sqrt(dnm);
                txfinal=tx;
                tyfinal=ty;
                nnarm=L(nn)/2;
                mmarm=L(mm)/2;
            end
            tx=x1(mm)-x2(nn);
            ty=y1(mm)-y2(nn);
            dnm=tx^2+ty^2;
            if(dnm<Dnm^2 && dnm<dnmfinal^2)
                dnmfinal=sqrt(dnm);
                txfinal=tx;
                tyfinal=ty;
                nnarm=-L(nn)/2;
                mmarm=L(mm)/2;
            end
            
            tx=x2(mm)-x2(nn);
            ty=y2(mm)-y2(nn);
            dnm=tx^2+ty^2;
            if(dnm<Dnm^2 && dnm<dnmfinal^2)
                dnmfinal=sqrt(dnm);
                txfinal=tx;
                tyfinal=ty;
                nnarm=-L(nn)/2;
                mmarm=-L(mm)/2;
            end
            
            tx=x2(mm)-x1(nn);
            ty=y2(mm)-y1(nn);
            dnm=tx^2+ty^2;
            if(dnm<Dnm^2 && dnm<dnmfinal^2)
                dnmfinal=sqrt(dnm);
                txfinal=tx;
                tyfinal=ty;
                nnarm=L(nn)/2;
                mmarm=-L(mm)/2;
            end
            
            if(dnmfinal<Dnm)
                C=C+1;
                F=-K*(Dnm/dnmfinal-1);
                Fx(nn)=Fx(nn)+F.*txfinal;  % particle-particle Force Law
                Fx(mm)=Fx(mm)-F.*txfinal;
                Fy(nn)=Fy(nn)+F.*tyfinal;  % particle-particle Force Law
                Fy(mm)=Fy(mm)-F.*tyfinal;
                torque(nn)=torque(nn)+F*nnarm*(cos(theta(nn))*tyfinal-sin(theta(nn))*txfinal);
                torque(mm)=torque(mm)-F*mmarm*(cos(theta(mm))*tyfinal-sin(theta(mm))*txfinal);
                PE=PE+K/2*((Dnm-dnmfinal)^2);
            end

        end
    end
    
    PEw=0;
    
    Fx=Fx-K*(x1-R).*(x1<R).*(x1<x2);  % Left wall
    PEw=PEw+K/2*sum(((x1-R).*(x1<R).*(x1<x2)).^2);
    Fx=Fx-K*(x2-R).*(x2<R).*(x2<x1);
    PEw=PEw+K/2*sum(((x2-R).*(x2<R).*(x2<x1)).^2);
    
    C=C+sum(x1<R | x2<R);
    
    Fy=Fy-K*(y1-R).*(y1<R).*(y1<y2);  % Bottom wall
    PEw=PEw+K/2*sum(((y1-R).*(y1<R).*(y1<y2)).^2);
    Fy=Fy-K*(y2-R).*(y2<R).*(y2<y1);
    PEw=PEw+K/2*sum(((y2-R).*(y2<R).*(y2<y1)).^2);
    
    C=C+sum(y1<R | y2<R);
    
    Fx=Fx-K*(x1-(Lx-R)).*(x1>Lx-R).*(x1>x2);  % Right wall
    PEw=PEw+K/2*sum(((x1-(Lx-R)).*(x1>Lx-R).*(x1>x2)).^2);
    Fx=Fx-K*(x2-(Lx-R)).*(x2>Lx-R).*(x2>x1);
    PEw=PEw+K/2*sum(((x2-(Lx-R)).*(x2>Lx-R).*(x2>x1)).^2);
    
    C=C+sum((x1>Lx-R) | (x2>Lx-R));
    
    Fy=Fy-K*(y1-(Ly-R)).*(y1>Ly-R).*(y1>y2);
    PEw=PEw+K/2*sum(((y1-(Ly-R)).*(y1>Ly-R).*(y1>y2)).^2);
    Fy=Fy-K*(y2-(Ly-R)).*(y2>Ly-R).*(y2>y1); % Top wall
    PEw=PEw+K/2*sum(((y2-(Ly-R)).*(y2>Ly-R).*(y2>y1)).^2);
    
    C=C+sum((y1>Ly-R) | (y2>Ly-R));
    
    for i=1:N
        % Left wall
        torque(i)=torque(i)+L(i)/2*sin(theta(i))*K*(x1(i)-R(i)).*(x1(i)<R(i)).*(x1(i)<x2(i));
        torque(i)=torque(i)+(-L(i)/2)*sin(theta(i))*K*(x2(i)-R(i)).*(x2(i)<R(i)).*(x2(i)<x1(i));
        % Bottom wall
        torque(i)=torque(i)-L(i)/2*cos(theta(i))*K*(y1(i)-R(i)).*(y1(i)<R(i)).*(y1(i)<y2(i));
        torque(i)=torque(i)-(-L(i)/2)*cos(theta(i))*K*(y2(i)-R(i)).*(y2(i)<R(i)).*(y2(i)<y1(i));
        % Right wall
        torque(i)=torque(i)+L(i)/2*sin(theta(i))*K*(x1(i)-(Lx-R(i))).*(x1(i)>Lx-R(i)).*(x1(i)>x2(i));
        torque(i)=torque(i)+(-L(i)/2)*sin(theta(i))*K*(x2(i)-(Lx-R(i))).*(x2(i)>Lx-R(i)).*(x2(i)>x1(i));
        % Top wall
        torque(i)=torque(i)-L(i)/2*cos(theta(i))*K*(y1(i)-(Ly-R(i))).*(y1(i)>(Ly-R(i))).*(y1(i)>y2(i));
        torque(i)=torque(i)-(-L(i)/2)*cos(theta(i))*K*(y2(i)-(Ly-R(i))).*(y2(i)>(Ly-R(i))).*(y2(i)>y1(i));
    end
    
    ax=Fx./M;
    ay=Fy./M;
    alpha=torque./I;
       
     % translational KE
    KE=M/2*sum(vx.^2+vy.^2);

    % rotational KE
    KR=sum(I.*(w.^2)/2);
    
%     if(PE+PEw<Pthresh*N && Dh<0)
%         R=R+sqrt(Pstep/K)*ratio;
%         L=g.*R;
%         I=M*(L/2).^2;
%     elseif(PE+PEw<Ptol*N/100 && Dh>0)
%         Dl=R;
%         R=(Dl+Dh)/2;
%         L=g.*R;
%         I=M*(L/2).^2;
%     elseif(KE+KR<1e-15*N && C>N/2 && PE+PEw>Ptol*N)
%         if(Dh<0)
%             Dh=R;
%         end
%         if(Dl>0)
%             Dh=R;
%             R=(Dl+Dh)/2;
%             L=g.*R;
%             I=M*(L/2).^2;
%         else
%             R=R-sqrt(Pstep/K)*ratio;
%             L=g.*R;
%             I=M*(L/2).^2;
%         end
%     elseif(KE+KR<1e-25*N && C>N/2)
%         break; % packing done
%     end
    
    vx=vx+(ax_old+ax).*dt/2;  % second step in Verlet integration
    vy=vy+(ay_old+ay).*dt/2;
    w=w+(alpha_old+alpha).*dt/2;
    
    
    ax_old=ax;
    ay_old=ay;
    alpha_old=alpha;
    
    total(nt)=KE+KR+PE+PEw;
    KRs(nt)=KR;
    KEs(nt)=KE;
    PEs(nt)=PE;
    PEws(nt)=PEw;
    
 
%     if(nt>1 && total(nt)-total(nt-1)>1e-4)
%         clf;
%         h=zeros(1,2*N);
%         hl=zeros(1,2*N);
%         for np=1:N
%           h(np)=rectangle('Position',[x(np)+cos(theta(np))*L/2-R y(np)+sin(theta(np))*L/2-R Dnm Dnm],'Curvature',[1 1],'edgecolor','b');
%           h(np+N)=rectangle('Position',[x(np)-cos(theta(np))*L/2-R y(np)-sin(theta(np))*L/2-R Dnm Dnm],'Curvature',[1 1],'edgecolor','b');
%           hl(np)=line([x(np)+L/2*cos(theta(np))+R*cos(theta(np)-pi/2),x(np)-L/2*cos(theta(np))+R*cos(theta(np)-pi/2)],...
%             [y(np)+L/2*sin(theta(np))+R*sin(theta(np)-pi/2),y(np)-L/2*sin(theta(np))+R*sin(theta(np)-pi/2)],'Color','b');
%           hl(np+N)=line([x(np)+L/2*cos(theta(np))-R*cos(theta(np)-pi/2),x(np)-L/2*cos(theta(np))-R*cos(theta(np)-pi/2)],...
%             [y(np)+L/2*sin(theta(np))-R*sin(theta(np)-pi/2),y(np)-L/2*sin(theta(np))-R*sin(theta(np)-pi/2)],'Color','b');
%         end
%         axis('equal');
%         axis([0 Lx 0 Ly]);
%         text(x,y,num2str((1:N)'));
%         text(x+cos(theta).*L/2,y+sin(theta).*L/2,num2str(1))
%         text(x-cos(theta).*L/2,y-sin(theta).*L/2,num2str(2))
%         pause;
%     end
end    
