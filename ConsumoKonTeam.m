clc;
% clear all;
run linearizeTrack;

%%'Constantes'
fr_a=0.479; % [Area frontal m2]
cl=0.592; %[coeficiente de sustentacion []]
cd=0.079; %[coeficiente de arrastre []]
dr=0.508; %[diametro rueda m]
crr=0.0031; %[Coeficiente de rodadua []]
m=99; %[masa [kg]]

%%Caracteristicas Motor
h=0.00068; %[Inductancia H];
ohm=0.268; %[resistencia interna ohm]
pn1=0.450; %[Nm]
pn2=0.0005; %[Nm/rpm]
kv=7.5; %[V/rpm]
kt=60/(2*pi()*kv); %[Nm/A]

%%controlador
i_bat=15; %A
i_pwm=30;
v_bat=48;
tmax=6; %Nm
%% 
%%Sacar pendiente en funcion al largo
track_L=linearizedTrack(:, 1);
A=[track_L;track_L(1)];
A=A(2:length(A));
track_h=linearizedTrack(:, 2);
B=[track_h; track_h(1)];
B=B(2:length(B));
track_p=(track_h - B)./(track_L - A);
track_p=atan(track_p)*180/pi();
track_d=[track_L,track_h,track_p];
track_d=[track_d,subsampledTrack];
scatter(track_d(:,4),track_d(:,5),abs(track_d(:,3))*100,track_d(:,3),'filled')

%% 
dt=0.1;
t=0;
d=0;
v=0.00001;
simudata=[];
a=0.018;
i=0;
v_cru=7;
E_acu=0;
while d<=track_d(length(track_d),1)
    i=1+i;
    t=t+dt; %alm
    v_rpm=v*60/(pi()*dr); % alm
    pen=interp1(track_L,track_p,d); %alm
    v_mot=v_rpm/kv; %alm
    
    %perdidas en fuerzas
    t_mag=(v/dr*0.5)*pn2 + pn1;
    f_mag=t_mag/(dr*0.5);
    f_sust=1.29*0.5*cl*v*v*fr_a;
    f_arr=1.29*0.5*cd*v*v*fr_a;
    f_crr=(f_sust+m*9.81*cos(pen*pi()/180))*crr;
    f_g=m*9.81*sin(pen*pi()/180);
    f_n=f_g+f_crr+f_arr+f_mag;
    i_mot=0;
    p_con=0;
    if f_n<0
        f_a=(f_n)*-1;
    elseif v<(v_cru-0.1)
        f_a=(tmax*2/dr)-f_n;
        p_con=(f_a+f_n)*v;
        i_mot=((f_a+f_n)*0.5*dr)/kt;
    elseif v>(v_cru-0.1) & v<(v_cru+0.1) &  f_n<=(tmax*2/dr)
        f_a=0;
    elseif f_n>(tmax*2/dr)
        f_a=(tmax*2/dr)-f_n;
        p_con=(f_a+f_n)*v;
        i_mot=((f_a+f_n)*0.5*dr)/kt;
    else 
        f_a=(f_n)*-1;
    end
    v_pwm=i_mot*ohm+v_mot;
    if i_mot>0;
        eff_st=1-(i_mot*ohm/v_pwm);
        eff_t=(1-(t_mag/(i_mot*kt)))*eff_st;
        p_bat=p_con/eff_st;
    else
        eff_st=0;
        eff_t=0;
        p_bat=0;
    end
    E_bat=p_bat*dt;
    E_acu=E_acu+E_bat;
    d=v*dt+d;
    v=f_a*dt/(m*1.05)+v;
    simudata=[simudata;d,E_bat,E_acu,f_a,p_bat,v,pen];
end
    
'terminado'