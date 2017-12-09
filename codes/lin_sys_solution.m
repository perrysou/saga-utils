A1 = [-1/5 2 ; -1/2 -1/5 ] ;
A2 = [ -1/5 1/2 ; -2 -1/5 ] ;

delta_t = pi/2;
t0 = 0 ;
tf = 100*delta_t ;
time = t0:delta_t:tf;

Phi1 = expm(A1*delta_t) ;
x1 = zeros(2,length(time));
x1(:,1) = [0;1] ;
for k=2:length(time)
    x1(:,k) = Phi1*x1(:,k-1);
end

figure(1)
plot(x1(1,:),x1(2,:))
grid
axis(1.5*[-1 1 -1 1])
xlabel('x_1(t)')
ylabel('x_2(t)')
title('State trajectory when A = A_1')

Phi2 = expm(A2*delta_t) ;
x2 = zeros(2,length(time));
x2(:,1) = [0;1] ;
for k=2:length(time)
    x2(:,k) = Phi2*x2(:,k-1);
end

figure(2)
plot(x2(1,:),x2(2,:))
grid
axis(1.5*[-1 1 -1 1])
xlabel('x_1(t)')
ylabel('x_2(t)')
title('State trajectory when A = A_2')

x3 = zeros(2,length(time));
x3(:,1) = [0;1] ;
for k=2:length(time)
    t_interval = floor(time(k-1)/(pi/2));
    if (t_interval/2)==floor(t_interval/2)
           x3(:,k) = Phi2*x3(:,k-1);
    else
        x3(:,k) = Phi1*x3(:,k-1);
    end
end

figure(3)
plot(x3(1,:),x3(2,:))
grid
%axis(1.5*[-1 1 -1 1])
xlabel('x_1(t)')
ylabel('x_2(t)')
title('State trajectory when A switches beteween A_1 and A_2')

% G = ss(A1,[0;0],eye(2),[0;0]);
% y = lsim(G,0*time,time,[1;0]) ;
% figure(4)
% plot(y(:,1),y(:,2))
% grid
% axis(1.5*[-1 1 -1 1])
% xlabel('x_1(t)')
% ylabel('x_2(t)')
% title('State trajectory when A = A_1, using LSIM command')
%
% figure(5)
% plot(x1(1,:),x1(2,:),y(:,1),y(:,2))
% grid
% axis(1.5*[-1 1 -1 1])
% xlabel('x_1(t)')
% ylabel('x_2(t)')
% title('State trajectory when A = A_1')
% legend({'x(t) using Phi','x(t) using LSIM'})
