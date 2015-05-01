%% Molecular Dynamics Simulator
% Seoyun Soy Lee

% 4/30/15
% working to generate packings of two different sized rods.
% 4/29/15
% implement for rods of two different lengths/radii

% Forces are calculated by taking the projection if possible
% and only checking end-to-end if none of the projections work

%% Experimental Parameters
N1 = 10; % number of rods of size 1
N2 = 10; % number of rods of size 2
N = N1+N2; % total number of rods
R1 = .01; % radius of circles/distance away from line
R2 = .016;
R = [repmat(R1,1,N1), repmat(R2,1,N2)];
R0 = min([R1, R2]);
ratio = R/(R0);
L1 = .03; % length of center line
L2 = .048;
L = [repmat(L1,1,N1), repmat(L2,1,N2)];
g1=L1/R1;
g2=L2/R2;
g = [repmat(g1,1,N1), repmat(g2,1,N2)];
K = 10; % spring constant
M = 3; % mass of particles
I = M*(L/2).^2;
B=1; % drag

Lx=1;  % size of box
Ly=1; 

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
Nt=1000000; % max time

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
    
    ct=cos(theta);
    st=sin(theta);
    
    lx=L/2.*ct;
    ly=L/2.*st;
    
    C=0; % counts
    % check to see which distance is shortest between each particle
    for nn=1:N-1
        for mm=nn+1:N
          Dnm=R(nn)+R(mm);
            % End of rod mm hits somewhere on rod nn
          p1mmsnn=false;
          nnm=(x1(mm)-x(nn))*ct(nn)+(y1(mm)-y(nn))*st(nn);
          if(abs(nnm) < L(nn)/2)
            tx=-(-(x1(mm)-x(nn))+nnm*ct(nn));
            ty=-(-(y1(mm)-y(nn))+nnm*st(nn));
            dnm=tx.^2+ty.^2;
            if(dnm<Dnm^2)
              p1mmsnn=true;
              dnm=sqrt(dnm);
              F=-K*(Dnm/dnm-1);
              Fx(nn)=Fx(nn)+F.*tx;  % particle-particle Force Law
              Fx(mm)=Fx(mm)-F.*tx;
              Fy(nn)=Fy(nn)+F.*ty;  % particle-particle Force Law
              Fy(mm)=Fy(mm)-F.*ty;
              torque(nn)=torque(nn)+F*nnm*(ct(nn)*ty-st(nn)*tx);
              torque(mm)=torque(mm)-F*L(mm)/2*(ct(mm)*ty-st(mm)*tx);
              PE=PE+K/2*((Dnm-dnm)^2);
            end
          end

          p2mmsnn=false;
          nnm=(x2(mm)-x(nn))*ct(nn)+(y2(mm)-y(nn))*st(nn);
          if(abs(nnm) < L(nn)/2)
            tx=-(-(x2(mm)-x(nn))+nnm*ct(nn));
            ty=-(-(y2(mm)-y(nn))+nnm*st(nn));
            dnm=tx.^2+ty.^2;
            if(dnm<Dnm^2)
              p2mmsnn=true;
              dnm=sqrt(dnm);
              F=-K*(Dnm/dnm-1);
              Fx(nn)=Fx(nn)+F.*tx;  % particle-particle Force Law
              Fx(mm)=Fx(mm)-F.*tx;
              Fy(nn)=Fy(nn)+F.*ty;  % particle-particle Force Law
              Fy(mm)=Fy(mm)-F.*ty;
              torque(nn)=torque(nn)+F*nnm*(ct(nn)*ty-st(nn)*tx);
              torque(mm)=torque(mm)-F*(-L(mm)/2)*(ct(mm)*ty-st(mm)*tx);
              PE=PE+K/2*((Dnm-dnm)^2);
            end
          end

          p1nnsmm=false;
          mmn=(x1(nn)-x(mm))*ct(mm)+(y1(nn)-y(mm))*st(mm);
          if(abs(mmn) < L(mm)/2)
            tx=-((x1(nn)-x(mm))-mmn*ct(mm));
            ty=-((y1(nn)-y(mm))-mmn*st(mm));
            dnm=tx.^2+ty.^2;
            if(dnm<Dnm^2)
              p1nnsmm=true;
              dnm=sqrt(dnm);
              F=-K*(Dnm/dnm-1);
              Fx(nn)=Fx(nn)+F.*tx;  % particle-particle Force Law
              Fx(mm)=Fx(mm)-F.*tx;
              Fy(nn)=Fy(nn)+F.*ty;  % particle-particle Force Law
              Fy(mm)=Fy(mm)-F.*ty;
              torque(nn)=torque(nn)+F*L(nn)/2*(ct(nn)*ty-st(nn)*tx);
              torque(mm)=torque(mm)-F*mmn*(ct(mm)*ty-st(mm)*tx);
              PE=PE+K/2*((Dnm-dnm)^2);
            end
          end

          p2nnsmm=false;
          mmn=(x2(nn)-x(mm))*ct(mm)+(y2(nn)-y(mm))*st(mm);
          if(mmn < L(mm)/2 && mmn > -L(mm)/2)
            tx=-((x2(nn)-x(mm))-mmn*ct(mm));
            ty=-((y2(nn)-y(mm))-mmn*st(mm));
            dnm=tx.^2+ty.^2;
            if(dnm<Dnm^2)
              p2nnsmm=true;
              dnm=sqrt(dnm);
              F=-K*(Dnm/dnm-1);
              Fx(nn)=Fx(nn)+F.*tx;  % particle-particle Force Law
              Fx(mm)=Fx(mm)-F.*tx;
              Fy(nn)=Fy(nn)+F.*ty;  % particle-particle Force Law
              Fy(mm)=Fy(mm)-F.*ty;
              torque(nn)=torque(nn)+F*(-L(nn)/2)*(ct(nn)*ty-st(nn)*tx);
              torque(mm)=torque(mm)-F*mmn*(ct(mm)*ty-st(mm)*tx);
              PE=PE+K/2*((Dnm-dnm)^2);
            end
          end

          C=C+(p1mmsnn || p1nnsmm || p2nnsmm || p2mmsnn);

          if(~p1mmsnn && ~p1nnsmm)
            tx=x1(mm)-x1(nn);
            ty=y1(mm)-y1(nn);
            dnm=tx.^2+ty.^2;
            if(dnm<Dnm^2)
              dnm=sqrt(dnm);
              F=-K*(Dnm/dnm-1);
              Fx(nn)=Fx(nn)+F.*tx;  % particle-particle Force Law
              Fx(mm)=Fx(mm)-F.*tx;
              Fy(nn)=Fy(nn)+F.*ty;  % particle-particle Force Law
              Fy(mm)=Fy(mm)-F.*ty;
              torque(nn)=torque(nn)+F*(lx(nn)*ty-ly(nn)*tx);
              torque(mm)=torque(mm)-F*(lx(mm)*ty-ly(mm)*tx);
              PE=PE+K/2*((Dnm-dnm)^2);
              C=C+1;
            end
          end

          if(~p1mmsnn && ~p2nnsmm)
            tx=x1(mm)-x2(nn);
            ty=y1(mm)-y2(nn);
            dnm=tx.^2+ty.^2;
            if(dnm<Dnm^2)
              dnm=sqrt(dnm);
              F=-K*(Dnm/dnm-1);
              Fx(nn)=Fx(nn)+F.*tx;  % particle-particle Force Law
              Fx(mm)=Fx(mm)-F.*tx;
              Fy(nn)=Fy(nn)+F.*ty;  % particle-particle Force Law
              Fy(mm)=Fy(mm)-F.*ty;
              torque(nn)=torque(nn)-F*(lx(nn)*ty-ly(nn)*tx);
              torque(mm)=torque(mm)-F*(lx(mm)*ty-ly(mm)*tx);
              PE=PE+K/2*((Dnm-dnm)^2);
              C=C+1;
            end
          end

          if(~p2mmsnn && ~p2nnsmm)
            tx=x2(mm)-x2(nn);
            ty=y2(mm)-y2(nn);
            dnm=tx.^2+ty.^2;
            if(dnm<Dnm^2)
              dnm=sqrt(dnm);
              F=-K*(Dnm/dnm-1);
              Fx(nn)=Fx(nn)+F.*tx;  % particle-particle Force Law
              Fx(mm)=Fx(mm)-F.*tx;
              Fy(nn)=Fy(nn)+F.*ty;  % particle-particle Force Law
              Fy(mm)=Fy(mm)-F.*ty;
              torque(nn)=torque(nn)-F*(lx(nn)*ty-ly(nn)*tx);
              torque(mm)=torque(mm)+F*(lx(mm)*ty-ly(mm)*tx);
              PE=PE+K/2*((Dnm-dnm)^2);
              C=C+1;
            end
          end

          if(~p2mmsnn && ~p1nnsmm)
            tx=x2(mm)-x1(nn);
            ty=y2(mm)-y1(nn);
            dnm=tx.^2+ty.^2;
            if(dnm<Dnm^2)
              dnm=sqrt(dnm);
              F=-K*(Dnm/dnm-1);
              Fx(nn)=Fx(nn)+F.*tx;  % particle-particle Force Law
              Fx(mm)=Fx(mm)-F.*tx;
              Fy(nn)=Fy(nn)+F.*ty;  % particle-particle Force Law
              Fy(mm)=Fy(mm)-F.*ty;
              torque(nn)=torque(nn)+F*(lx(nn)*ty-ly(nn)*tx);
              torque(mm)=torque(mm)+F*(lx(mm)*ty-ly(mm)*tx);
              PE=PE+K/2*((Dnm-dnm)^2);
              C=C+1;
            end
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
  
  C=C+sum(x1<R)+sum(y1<R)+sum(x2<R)+sum(y2<R);
  C=C+sum(x1>Lx-R)+sum(x2>Lx-R)+sum(y1>Ly-R)+sum(y2>Ly-R); % total number of contacts
  
  for i=1:N
    % Left wall
    torque(i)=torque(i)+L(i)/2*sin(theta(i))*K*(x1(i)-R(i)).*(x1(i)<R(i));
    torque(i)=torque(i)+(-L(i)/2)*sin(theta(i))*K*(x2(i)-R(i)).*(x2(i)<R(i));
    % Bottom wall
    torque(i)=torque(i)-L(i)/2*cos(theta(i))*K*(y1(i)-R(i)).*(y1(i)<R(i));
    torque(i)=torque(i)-(-L(i)/2)*cos(theta(i))*K*(y2(i)-R(i)).*(y2(i)<R(i));
    % Right wall
    torque(i)=torque(i)+L(i)/2*sin(theta(i))*K*(x1(i)-(Lx-R(i))).*(x1(i)>Lx-R(i));
    torque(i)=torque(i)+(-L(i)/2)*sin(theta(i))*K*(x2(i)-(Lx-R(i))).*(x2(i)>Lx-R(i));
    % Top wall
    torque(i)=torque(i)-L(i)/2*cos(theta(i))*K*(y1(i)-(Ly-R(i))).*(y1(i)>(Ly-R(i)));
    torque(i)=torque(i)-(-L(i)/2)*cos(theta(i))*K*(y2(i)-(Ly-R(i))).*(y2(i)>(Ly-R(i)));
  end
  
    ax=Fx./M-B*vx;
    ay=Fy./M-B*vy;
    alpha=torque./I-B*w;
    
    if(PE+PEw<Pthresh*N && Dh<0)
        R0=R0+sqrt(Pstep/K)/2;
        R=R0*ratio;
        L=g.*R;
        I=M*(L/2).^2;
    elseif(PE+PEw<Ptol*N/100 && Dh>0)
        Dl=R0;
        R0=(Dl+Dh)/2;
        R=R0*ratio;
        L=g.*R;
        I=M*(L/2).^2;
    elseif(KE+KR<1e-15*N && C>N/2 && PE+PEw>Ptol*N)
        if(Dh<0)
            Dh=R0;
        end
        if(Dl>0)
            Dh=R0;
            R0=(Dl+Dh)/2;
            R=R0*ratio;
            L=g.*R;
            I=M*(L/2).^2;
        else
            R0=R0-sqrt(Pstep/K)/2;
            R=R0*ratio;
            L=g.*R;
            I=M*(L/2).^2;
        end
    elseif(KE+KR<1e-25*N && C>N/2)
        break; % packing done
    end
    
     % translational KE
    KE=M/2*sum(vx.^2+vy.^2);

    % rotational KE
    KR=sum(I.*(w.^2)/2);
    
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
end    
