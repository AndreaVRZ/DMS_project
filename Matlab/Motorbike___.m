close all
clear all

load('Motorbike_mkr.mat')
%--------------------------------------------------------------------------
%Dividiamo le matrici appena caricate in sottomatrici
C=R;

MFF=M(1:39,1:39);
MFC=M(1:39,40:47);
MCF=M(40:47,1:39);
MCC=M(40:47,40:47);

KFF=K(1:39,1:39);
KFC=K(1:39,40:47);
KCF=K(40:47,1:39);
KCC=K(40:47,40:47);

CFF=C(1:39,1:39);
CFC=C(1:39,40:47);
CCF=C(40:47,1:39);
CCC=C(40:47,40:47);
%--------------------------------------------------------------------------
%Ricalcoliamo le FREQUENZE NATURALI e i MODI DI VIBRARE come verifica
[modi,autovalori]=eig(MFF\KFF);
freq=sqrt(diag(autovalori))/(2*pi);
modi;
%--------------------------------------------------------------------------
%PUNTO 4 a)
%Definiamo il vettore di frequenze in cui effettuare l'analisi
vett_f=0:0.1:300; %f_min=0Hz f_max=300Hz passo=0.1Hz
i=sqrt(-1);
%Costruiamo il vettore dei vincoli dal 40 al 47, inserendo 1 solo nei 
%vincoli con movimento armonico
y_davanti=1;
y_dietro=1;
x_c=[0; y_dietro; 0; 0; y_davanti; 0; 0; 0];  

for j=1:length(vett_f)
    omega=2*pi*vett_f(j);
    A=-omega^2*MFF+i*omega*CFF+KFF;
    Q0=-(-omega^2*MFC+i*omega*CFC+KFC)*x_c;
    x0=A\Q0;
    
    %Vogliamo l'accelerazione verticale della sella: guardando la IDB
    %matrix risulta il grado di libertà 34 
    y_sella=x0(34);
    y_selladd=-omega^2*y_sella;
    
    modulo_y_selladd(j)=abs(y_selladd);
    fase_y_selladd(j)=(angle(y_selladd)*180)/pi;
end

figure
subplot 211; plot(vett_f,modulo_y_selladd); grid; xlabel('Freq. [Hz]'); ylabel('Ampl. [m/s^2/m]'); title('Accelerazione verticale della sella - nodo 16, g.d.l. 34')
subplot 212; plot(vett_f,fase_y_selladd); grid; xlabel('Freq. [Hz]'); ylabel('Fase [deg]')

%PUNTO 4 b)
vett_f=0:0.1:300;
i=sqrt(-1);
phi=0.25*pi;
y_davanti=1;
y_dietro=exp(i*phi);
x_c=[0; y_dietro; 0; 0; y_davanti; 0; 0; 0];

for j=1:length(vett_f)
    omega=2*pi*vett_f(j);
    A=-omega^2*MFF+i*omega*CFF+KFF;
    Q0=-(-omega^2*MFC+i*omega*CFC+KFC)*x_c;
    x0=A\Q0;
    
    y_sella=x0(34);
    y_selladd=-omega^2*y_sella;
    
    modulo_y_selladd(j)=abs(y_selladd);
    fase_y_selladd(j)=(angle(y_selladd)*180)/pi;
end

figure
subplot 211; plot(vett_f,modulo_y_selladd); grid; xlabel('Freq. [Hz]'); ylabel('Ampl. [m/s^2/m]'); title('Accelerazione verticale della sella - nodo 16, g.d.l. 34')
subplot 212; plot(vett_f,fase_y_selladd); grid; xlabel('Freq. [Hz]'); ylabel('Fase [deg]')
%--------------------------------------------------------------------------
%PUNTO 5
A1=0.13;
A2=0.009;
A3=0.0025;
lambda1=10;
lambda2=1;
lambda3=0.1;
%Trasformiamo lo spazio nel tempo (spazio=velocità*tempo)
V=180/3.6; %conversione velocità km/h --> m/s
t1=lambda1/V;
t2=lambda2/V;
t3=lambda3/V;
t=0:0.001:t1; %vettore dei tempi (t1 è il maggiore dei 3!)
y1=A1*cos(2*pi/lambda1*V*t);
y2=A2*cos(2*pi/lambda2*V*t);
y3=A3*cos(2*pi/lambda3*V*t);
y_davanti=y1+y2+y3;

figure; plot(t,y_davanti); grid; xlabel('Tempo [s]'); ylabel('Spostamento [m]'); title('y\_davanti')

%Stesso procedimento anche per la ruota posteriore tenendo in
%considerazione lo sfasamento temporale
t0=(1.560-0.600)/V;
N=length(t);
t_sfasamento=t0*ones(1,N); %vettore di sfasamenti
y1=A1*cos(2*pi/lambda1*V*(t-t_sfasamento));
y2=A2*cos(2*pi/lambda2*V*(t-t_sfasamento));
y3=A3*cos(2*pi/lambda3*V*(t-t_sfasamento));
y_dietro=y1+y2+y3;

figure; plot(t,y_dietro); grid; xlabel('Tempo [s]'); ylabel('Spostamento [m]'); title('y\_dietro')

%PUNTO 5a
%SPETTRI degli spostamenti impressi alle ruote
df=1/t1;
N=length(t);
fmax=df*(N/2-1); %Shannon
vett_f=0:df:fmax;
fftout=fft(y_davanti);
modulo_davanti(1)=1/N*abs(fftout(1));
modulo_davanti(2:N/2)=2/N*abs(fftout(2:N/2));
fase_davanti(1:N/2)=angle(fftout(1:N/2))*180/pi;

figure
subplot 211; bar(vett_f,modulo_davanti); grid; xlabel('Freq. [Hz]'); ylabel('Modulo'); title('Spettro y\_davanti')
subplot 212; bar(vett_f,fase_davanti); grid; xlabel('Freq. [Hz]'); ylabel('Fase [deg]');
%-----------
df=1/t1;
N=length(t);
fmax=df*(N/2-1); %Shannon
vett_f=0:df:fmax;
fftout=fft(y_dietro);
modulo_dietro(1)=1/N*abs(fftout(1));
modulo_dietro(2:N/2)=2/N*abs(fftout(2:N/2));
fase_dietro(1:N/2)=angle(fftout(1:N/2))*180/pi;

figure
subplot 211; bar(vett_f,modulo_dietro); grid; xlabel('Freq. [Hz]'); ylabel('Modulo'); title('Spettro y\_dietro')
subplot 212; bar(vett_f,fase_dietro); grid; xlabel('Freq. [Hz]'); ylabel('Fase [deg]');

%PUNTO 5b
for k=1:N/2
    omega=2*pi*vett_f(k);
    y_davanti=modulo_davanti(k)*exp(i*fase_davanti(k));
    y_dietro=modulo_dietro(k)*exp(i*fase_dietro(k));
    x_c=[0; y_dietro; 0; 0; y_davanti; 0; 0; 0];
    A=-omega^2*MFF+i*omega*CFF+KFF;
    Q0=-(-omega^2*MFC+i*omega*CFC+KFC)*x_c;
    x0=A\Q0;
    
    y_sella=x0(34);
    
    modulo_sella(k)=abs(y_sella);
    fase_sella(k)=angle(y_sella)*180/pi;
end

figure
subplot 211; plot(vett_f,modulo_sella); grid; xlabel('Freq. [Hz]'); ylabel('Ampl. [m/m]'); title('Spostamento verticale della sella - nodo 16, g.d.l. 34')
subplot 212; plot(vett_f,fase_sella); grid; xlabel('Freq. [Hz]'); ylabel('Fase [deg]')

%PUNTO 5c e 5d
spostamento_verticale_sella=zeros(1,N);
for k=1:N/2
    omega=2*pi*vett_f(k);
    y_davanti=modulo_davanti(k)*exp(i*fase_davanti(k));
    y_dietro=modulo_dietro(k)*exp(i*fase_dietro(k));
    x_c=[0; y_dietro; 0; 0; y_davanti; 0; 0; 0];
    A=-omega^2*MFF+i*omega*CFF+KFF;
    Q0=-(-omega^2*MFC+i*omega*CFC+KFC)*x_c;
    x0=A\Q0;
    
    y_sella=x0(34);
    
    modulo_sella(k)=abs(y_sella);
    fase_sella(k)=angle(y_sella)*180/pi;
    
    spostamento_verticale_sella=spostamento_verticale_sella+modulo_sella(k)*cos(omega*t+fase_sella(k));
end

figure; plot(t,spostamento_verticale_sella); grid; title('Storia temporale dello spostamento verticale della sella')

df=1/t1;
N=length(t);
fmax=df*(N/2-1);
vett_f=0:df:fmax;
fftout=fft(spostamento_verticale_sella);
modulo(1)=1/N*abs(fftout(1));
modulo(2:N/2)=2/N*abs(fftout(2:N/2));
fase(1:N/2)=angle(fftout(1:N/2))*180/pi;

figure
subplot 211; bar(vett_f,modulo); grid; title('Spettro dello spostamento della sella')
subplot 212; bar(vett_f,fase); grid
%----------
f_risonanza_sella=freq(34) %frequenza risonanza della sella
V1=f_risonanza_sella*lambda1
V2=f_risonanza_sella*lambda2
V3=f_risonanza_sella*lambda3
%--------------------------------------------------------------------------
%PUNTO 6
omega_motore=1000*2*pi/60;
f_motore=omega_motore/2/pi

m_motore=150*0.140

x_max=0.001;

k_e_min=(m_motore*9.81)/(2*x_max) %103500 N/m=10,35*10^4 N/m

f_naturale=sqrt((2*k_e_min)/m_motore)/(2*pi)

r=f_motore/f_naturale

k=200000; %valore assegnato nei dati
ome_naturale=sqrt((2*k)/m_motore);
c=2*m_motore*ome_naturale %c_cr=2*m*omega smorzamento critico