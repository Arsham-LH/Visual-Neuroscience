function [S] = WANG_edit(X1)

% StimDur=.5;
StimDur=1.5;
% endTime=2.5;
endTime=1.5;
 JN11=0.2609;
 JN12 =0.0497;
 dt=.001;
 SDnoise=.02;

thr=X1(1);
coh=X1(2);
I0=X1(3);
u0=X1(4);

JN22=JN11;
JN21=JN12;
Tampa=2e-3;
Ts=100e-3;
gama=.641;
a= 270;
b=108;
d=.154;
Ja=5.2e-4;

S1=0;
S2=0;
Inoise1=0;
Inoise2=0;

RT=0;
deltaS=0;
choice=0;

T1=((Tampa/dt*SDnoise^2)^.5);
I1=Ja*u0*(1+coh/100);
I2=Ja*u0*(1-coh/100);

StimDurS=floor(StimDur/dt);
InoStim=Ja*u0*(1+0);
ti0=endTime/dt;
% MS1=zeros();
% MS2=zeros();
for i0=1:floor(ti0)
    if i0==StimDurS 
        I1=InoStim;
        I2=InoStim;
    end
    
    Inoise1=Inoise1+(dt/Tampa)*( -Inoise1-randn()*T1  );
    Inoise2=Inoise2+(dt/Tampa)*( -Inoise2-randn()*T1 );
    
    x1=JN11*S1-JN12*S2+I0+I1+Inoise1;
    x2=JN22*S2-JN21*S1+I0+I2+Inoise2;
    
    H1=(a*x1-b)/(1-exp(-d*(a*x1-b)) );
    H1(H1 < 0) = 0; 
    S1=S1+(dt)*(-S1/Ts+(1-S1)*gama*H1);
    
    H2=(a*x2-b)/(1-exp(-d*(a*x2-b)) );
    H2(H2 < 0) = 0; 
    S2=S2+(dt)*(-S2/Ts+(1-S2)*gama*H2);
    
    if (S1>thr || S2>thr) && (RT==0)
        RT=i0*dt;
        if S1>S2
           choice=1;
        end
        deltaS=abs(S2-S1);
        break;
    end
% MS1(i0)=S1; 
% MS2(i0)=S2;        

end
if RT==0
    RT=endTime;
    if S1>S2
       choice=1;
    end
    deltaS=abs(S2-S1);
end
% plot(MS1)
% hold on
% plot(MS2)

S=[RT deltaS choice];
end