clear all
close all
%--------------------------------------------------------------------------

load('Motorbike_mkr.mat')

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
[modi autovalori]=eig(MFF\KFF);
frequenze_naturali=sqrt(diag(autovalori))/2/pi;
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
x_c=[0; y_dietro; 0; 0; y_davanti; 0; 0; 0];  % Costruisco il vettore dei vincoli dal 40 al 47, inserendo 1 solo nei vincoli con movimento armonico
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
subplot 212; plot(vett_f,fase_y_selladd); grid; xlabel('Freq. [Hz]'); ylabel('Fase [rad]')
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
T1=lambda1/V;
T2=lambda2/V;
T3=lambda3/V;
t=0:0.001:T1; %vettore dei tempi (T1 è il maggiore dei 3!) 
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

%PUNTO 5 a)
%SPETTRI degli spostamenti impressi alle ruote
df=1/T1;
N=length(t);
fmax=df*(N/2-1); %Shannon
vett_freq=0:df:fmax;
fftout=fft(y_davanti);
modulo_davanti(1)=1/N*abs(fftout(1));
modulo_davanti(2:N/2)=2/N*abs(fftout(2:N/2));
fase_davanti(1:N/2)=angle(fftout(1:N/2));

figure
subplot 211; bar(vett_freq,modulo_davanti); grid; xlabel('Freq. [Hz]'); ylabel('Modulo'); title('Spettro y\_davanti')
subplot 212; bar(vett_freq,fase_davanti); grid; xlabel('Freq. [Hz]'); ylabel('Fase [rad]');
%-----------
df=1/T1;
N=length(t);
fmax=df*(N/2-1); %Shannon
vett_freq=0:df:fmax;
fftout=fft(y_dietro);
modulo_dietro(1)=1/N*abs(fftout(1));
modulo_dietro(2:N/2)=2/N*abs(fftout(2:N/2));
fase_dietro(1:N/2)=angle(fftout(1:N/2));

figure
subplot 211; bar(vett_freq,modulo_dietro); grid; xlabel('Freq. [Hz]'); ylabel('Modulo'); title('Spettro y\_dietro')
subplot 212; bar(vett_freq,fase_dietro); grid; xlabel('Freq. [Hz]'); ylabel('Fase [rad]');

%PUNTO 5 b)
for k=1:N/2
    omega=2*pi*vett_freq(k);
    y_davanti=modulo_davanti(k)*exp(i*fase_davanti(k));
    y_dietro=modulo_dietro(k)*exp(i*fase_dietro(k));
    x_c=[0; y_dietro; 0; 0; y_davanti; 0; 0; 0];
    A=-omega^2*MFF+i*omega*CFF+KFF;
    Q0=-(-omega^2*MFC+i*omega*CFC+KFC)*x_c;
    x0=A\Q0;
    
    y_sella=x0(34);
    
    modulo_sella(k)=abs(y_sella);
    fase_sella(k)=(angle(y_sella)*180)/pi;
end

figure
subplot 211; plot(vett_freq,modulo_sella); grid; xlabel('Freq. [Hz]'); ylabel('Ampl. [m/m]'); title('Spostamento verticale della sella - nodo 16, g.d.l. 34')
subplot 212; plot(vett_freq,fase_sella); grid; xlabel('Freq. [Hz]'); ylabel('Fase [deg]') 

% Punto 5 c) e d)
spostamento_verticale_sella=zeros(1,N);
for k=1:N/2
    omega=2*pi*vett_freq(k);
    y_davanti=modulo_davanti(k)*exp(i*fase_davanti(k));
    y_dietro=modulo_dietro(k)*exp(i*fase_dietro(k));
    x_c=[0; y_dietro; 0; 0; y_davanti; 0; 0; 0];
    A=-omega^2*MFF+i*omega*CFF+KFF;
    Q0=-(-omega^2*MFC+i*omega*CFC+KFC)*x_c;
    x0=A\Q0;
    
    y_sella=x0(34);
    
    modulo_sella(k)=abs(y_sella);
    fase_sella(k)=angle(y_sella);
    
    spostamento_verticale_sella=spostamento_verticale_sella+modulo_sella(k)*cos(omega*t+fase_sella(k));
end

figure; plot(t,spostamento_verticale_sella); grid; xlabel('Tempo [s]'); ylabel('Spostamento [m]'); title('Storia temporale dello spostamento verticale della sella')

df=1/T1;
N=length(t);
fmax=df*(N/2-1);
vett_freq=0:df:fmax;
fftout=fft(spostamento_verticale_sella);
modulo(1)=1/N*abs(fftout(1));
modulo(2:N/2)=2/N*abs(fftout(2:N/2));
fase(1:N/2)=angle(fftout(1:N/2));

figure
subplot 211; bar(vett_freq,modulo);  grid; xlabel('Freq. [Hz]'); ylabel('Modulo'); title('Spettro dello spostamento della sella')
subplot 212; bar(vett_freq,fase); grid; xlabel('Freq. [Hz]'); ylabel('Fase [rad]'); 

% Ultima domanda del PUNTO 5
% Sappiamo che (dal vettore delle frequenze naturali) la frequenza di
% risonanza verticale della sella (grado di libertà 34) è
% [35.4415418843955] (da frequenze_naturali)
% Vogliamo quindi che lo spostamento imposto alle ruote abbia al suo
% interno una frequenza pari a quella appena citata. Abbiamo tre possibili
% frequenze per lo spostamento verticale delle ruote
f_risonanza_sella=35.4415;
% Imponiamo frequenza_principale=f_risonanza_sella!
V1=f_risonanza_sella*lambda1
V2=f_risonanza_sella*lambda2
V3=f_risonanza_sella*lambda3
%--------------------------------------------------------------------------
%PUNTO 6
omega_motore=1000*2*pi/60;
freq_motore=omega_motore/2/pi
% Quindi ho una frequenza pari a 16,67 Hz sul motore.
% Ci da un limite sulla deflessione statica massima delle sospensioni del
% motore: 1 mm. Vuol dire che il peso del motore, diviso su entrambe le
% molle, deve portare ad una deflessione statica massima di 1 mm. Per cui:
% m/2*g=k_e*x. m_motore=0,14*150=21 kg (calcolato dalla densità della trave
% del motore). Quindi k_e_min=m*g/(2*x_max)=102900 N/m=10,29*10^4

% Per la tecnica di isolamento delle vibrazioni, cerchiamo di rendere il
% coefficiente di trasmissibilità il più piccolo possibile (vedi grafico su
% word). Per fare questo dobbiamo rendere f/fn maggiore di sqrt(2) (più è 
% grande meglio è), quindi dobbiamo rendere la frequenza naturale il più
% piccolo possibile rispetto alla frequenza di eccitazione (pari a 16,67 Hz,
% ricordiamolo). Per rendere la frequenza naturale il più piccolo
% possibile, devo rispettare la richiesta della deflessione massima e
% prendere k_e=k_e_min=102900 N/m. Per cui 
k_e_min=102900;
m_motore=21;
freq_naturale=sqrt(2*k_e_min/m_motore)/2/pi

% Quindi il rapporto tra le frequenze (r) è pari a
r=freq_motore/freq_naturale

% Non possiamo ottenere nessun risultato cambiando il valore della
% rigidezza, quindi manteniamolo al valore precedente e cambiamo solo il
% valore dello smorzamento (ponendo uguale al valore critico) in modo da
% abbassare il più possibile il picco del coefficiente di trasmissibilità T
k=200000; %valore assegnato nei dati
pulsazione_naturale=sqrt(2*k/m_motore)
c=2*m_motore*pulsazione_naturale


    
    
    





