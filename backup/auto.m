clear all;  clc; close all;
global a b1 b2 b3 k omega theta0 c0 M0 M1 M2 mu xi0 z0 L T A nu x0;

%строит зависимость ошибки от шага
    L=120; T=16; 
    M0=-1.48; M1=6.16; k=1.6; 
    
    mu=sqrt(M1/(M0*(M1+6*M0)));             
    omega = k^2-1/12 *(M0*M1)/(M1+6*M0)-1/6 *M0;
    b2=3/(4*M0)*(M1/(M1+6*M0)-2);
    b3=4/(M0*(M1+6*M0));
    c0=2*k; 
    theta0=0; xi0=0; z0=-20;
    
    aa=log((sqrt(4*M0*M1 + M1^2)-2*M0-M1)/(2*M0))/mu+xi0; bb=log((-sqrt(4*M0*M1 + M1^2)-2*M0-M1)/(2*M0))/mu+xi0;
    aaa=aa+(10^-14)/1; bbb=bb-(10^-14)/1;
    middle=(aaa+bbb)/2; alog=-2; blog=1; numpoints=10000;
    XI1=(logspace(alog,blog,numpoints)-10^alog)*(middle-aaa)/(10^blog-10^alog)+aaa;
    XI2=(logspace(blog,alog,numpoints)-10^blog)*(bbb-middle)/(10^alog-10^blog)+middle;
    XI=cat(2,XI1(1:end-1),XI2);
    Z=zofxi(XI);
    [Z,XI]=fix(Z,XI,14,aa,bb);
    YY=yofxi(XI);

H=[1 0.5 0.4 0.25];
Tau=H.^2;
for kk=1:max(size(H))
    h=H(kk); tau=Tau(kk);
    N=L/h; Nx=N+1; Nt=round(T/tau)+1;
    j=linspace(-N/2,N/2,Nx);
    x=j*h; t=linspace(0,T,Nt);
    U=zeros(Nx-1,Nt);
    Utrue=zeros(Nx-1,Nt);
    V=zeros(Nx-1,Nt); 
mun=2*pi/L .*linspace(-N/2,N/2-1,N); %частоты
M=exp(-1i.*mun'.^2 .*tau); %соотношение между фурье коэфтами 
direct=(h/L.*exp(-1i.*mun'*x(1:end-1))); inverse=exp(1i.*mun'*x(1:end-1));

    U(:,1) = arrayfun(@(x) initial(x,Z,XI) ,x(1:end-1)); 
    for i=1:1:(max(size(x))-1)
        Utrue(i,:)=true(x(i),t,Z,XI);
    end
    for m=progress(1:Nt-1)
        V(:,m)=exp(1i*tau* (1*(abs(U(:,m))).^2 + b2*(abs(U(:,m))).^4 + b3*(abs(U(:,m))).^6 ) ).*U(:,m);
        Unf=M.*(direct*V(:,m));
        U(:,m+1)=inverse*Unf;
    end
    err=zeros(1,Nt);
    for m=1:Nt
        err(m)=max(abs(U(:,m))-abs(Utrue(:,m)));
    end
    Err(kk)=max(err);
end

figure('Name','ошибка от сетки'); hold on; grid on;

plot(H,Err,'k^-','LineWidth',1);
set ( gca, 'xdir', 'reverse' )
xlabel('h','FontSize',12); ylabel('maximum relative error','FontSize',12);

% Mc=round(Tc./tau)+1;
% for i=1:max(size(Tc))
%     figure('Name',['что есть при T= ',num2str(Tc(i))]); hold on; grid on;
%     ylim([0 ymax]); xlim([-L/2 L/2]);
%     xlabel('x','FontSize',12); ylabel('|U|','FontSize',12,'Rotation',0);
% 
%     plot(x(1:end-1),abs(U(:,Mc(i))),'k','LineWidth',1);
%     plot(x(1:end-1),abs(Utrue(:,Mc(i))),'r--','LineWidth',2    );
%     legend(join(['numerical solution']),join(['analytical solution']),'FontSize',12);
%     ax=gca; ax.FontSize=12;
%     xticks(linspace(-L/2,L/2,6));
%     yticks(linspace(0,ymax,6));
% end


function res=initial(x,Z,XI)
    res=true(x,0,Z,XI);
    %res=0.5.*(1+0.1.*cos(pi.*x./8));
end
function res=true(x,t,Z,XI)
global M0 M1 k c0 xi0 L omega mu theta0;
    t(t>(L/2+x)/c0)=t(t>(L/2+x)/c0)-(L/c0)*floor(1/2+(t(t>(L/2+x)/c0).*c0-x)/L);
    z=x-c0.*t;
    xi = interp1(Z,XI,z);
    size(xi)
    res=((M0+M1./(1+exp(mu.*(xi-xi0)))-M1./((1+exp(mu.*(xi-xi0))).^2)).^(1/2) ).*exp(1i.*(k.*x-omega.*t-theta0));
end
function res=zofxi(xi)
global M0 M1 mu z0 xi0;
   res=z0+  xi./M0  +(2.*M1)./(mu.*M0.*sqrt(4.*M0.*M1+M1.^2)) .*atanh((2.*M0.*exp(mu.*(xi-xi0))+2.*M0+M1)./(sqrt(4.*M0.*M1+M1.^2)));
end
function res=yofxi(xi)
global M0 M1 mu xi0;
   res=(M0+M1./(1+exp(mu.*(xi-xi0)))-M1./((1+exp(mu.*(xi-xi0))).^2)).^(1/2);
end
function [Zn, Xin]=fix(Z,XI,GP,aa,bb)
global M0 M1 L;
    general_power=GP;
    if (~isreal(Z))
        disp('Z is complex!! fixing.')
        general_power=general_power-1;
        XI=linspace(aa+10^-general_power,bb-10^-general_power,100000);
        Z=zofxi(XI,M0,M1);
        if (~isreal(Z))
            disp('not today')
        end
    end
    if max(Z)==Inf
        disp('max is inf!!')
    end
    if min(Z)==-Inf
        disp('min is -inf!!')
    end
    k1=0; k2=0; rightflag=0;
    while ((L/2>max(Z)) && k1<5 && isreal(Z) &&  ~(max(Z)==Inf || min(Z)==-Inf))
        disp('right true expansion iteration')
        k1=k1+1;
        XI=linspace(aa+10^-general_power,bb-(10^-general_power)/(2^k1),100000);
        Z=zofxi(XI);
        if (max(Z)==Inf || ~isreal(Z))
            disp('right true expansion jumped to +inf or became complex. Restoring and using KTL')
            k1=k1-1;
            XI=linspace(aa+10^-general_power,bb-(10^-general_power)/(2^k1),100000);
            Z=zofxi(XI);
            Z(end+1)=10*L;
            XI(end+1)=bb;  
            rightflag=1;
            continue;
        end
    end
    while ((-L/2<min(Z)) && k2<5 && isreal(Z) &&  ~(max(Z)==Inf || min(Z)==-Inf))
        disp('left true expansion iteration')
        k2=k2+1;
        XI=linspace(aa+(10^-general_power)/(2^k2),bb-(10^-general_power)/(2^k1),100000);
        Z=zofxi(XI);

        if (min(Z)==-Inf || ~isreal(Z))
            disp('left true expansion jumped to -inf or became complex. Restoring and using KTL')
            k2=k2-1;
            XI=linspace(aa+(10^-general_power)/(2^k2),bb-(10^-general_power)/(2^k1),100000);
            Z=zofxi(XI);
            Z(end+1)=-10*L;
            XI(end+1)=aa;  
            continue;
        end
    end
    if rightflag==1
        Z(end+1)=10*L;
        XI(end+1)=bb; 
    end
    Zn=Z;
    Xin=XI;
end