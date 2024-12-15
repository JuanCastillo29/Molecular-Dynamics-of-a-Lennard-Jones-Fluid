#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sys/stat.h>
#include <string.h>

using namespace std;

unsigned int irr[256];
unsigned int ir1;
unsigned char ind_ran,ig1,ig2,ig3;
#define NormRANu (2.3283063671E-10F)
#define Pi 3.14159265
#define Gases 8.314462 //Constante de los gases en J/(K·mol)
#define Avogradro 6.022140 //Constante de Avogradro.

#define NMedidas 131072
#define Termostato 10
#define dt 0.001
#define TTermalizacion 100
#define NPart 128
#define TMeasure 0.05

double TTot = NMedidas*TMeasure;
const double Densidad[4]={0.2, 0.4, 0.6, 0.8};

int MedidasTotal=TTot/dt, Medida = TMeasure/dt, NTermalizacion = TTermalizacion/dt;

#define m 40  //Masa en g/mol.
#define e 0.998 //Energia en kJ/mol
#define sigma 3.4   //Radio molecular en amstrongs.

const double V = sigma*sigma*sigma; //Volumen en Amstrongs^3.
const double t = sqrt(m/e)*sigma/10; //Tiempo en picosegundos
const double Temp = 1000*e/Gases; //Temperatura en Kelvin
const double Pres = 100000*e/Avogradro/V; // Presion en GPa
const double rho = m/Avogradro/V*10; //Densidad en g/cm^3.

//Generador de números de Parisi-Rapuano.
void ini_ran(int SEMILLA)
{
    int INI,FACTOR,SUM,i;
    srand(SEMILLA);
    INI=SEMILLA;
    FACTOR=67397;
    SUM=7364893;
    for(i=0;i<256;i++)
    {
        INI=(INI*FACTOR+SUM);
        irr[i]=INI;
    }
    ind_ran=ig1=ig2=ig3=0;
}

float Random(void)
{
    float r;
    ig1=ind_ran-24;
    ig2=ind_ran-55;
    ig3=ind_ran-61;
    irr[ind_ran]=irr[ig1]+irr[ig2];
    ir1=(irr[ind_ran]^irr[ig3]);
    ind_ran++;
    r=ir1*NormRANu;
    return r;
}

//Generador de los números aleatorios según distribución gaussiana por el método de Box-Muller

double DistrGauss(void){
    double d1, d2;
    //Se generan los numeros planos necesarios
    d1 = Random();
    d2 = Random();
    //Se usa la fórmula para crear un número independiente según una gaussiana
    return -sqrt(-2*log(d1))*cos(2*Pi*d2);
}

double PBC(double x, double L){
if(x>L/2)
    x=x-L;
if(x<-L/2)
    x=x+L;
return x;
}

double IniBCC(double (&r)[NPart][3], double (&v)[NPart][3], double densi){
    int N = int(pow((NPart/2)*1.0001, 1.0/3));
    double l = pow((NPart/2)/densi, 1.0/3), a = l/N;
    for(int i=0; i<NPart/2; i++){
        r[i][2]= (i/(N*N));
        r[i][1]=(i%(N*N))/N;
        r[i][0] = (i)%N;
        for(int j=0; j<3; j++){
            r[i+NPart/2][j] = PBC(a*(r[i][j] + 0.5), l);
            r[i][j] = PBC(a*r[i][j], l);
        }
    }
    for(int i = 0; i<NPart; i++){
        for(int j=0; j<3; j++){
            v[i][j]=0;
        }
    }
return l;
}

void GuardarConfi(double (&r)[NPart][3], int d){
char filename[20];
sprintf(filename, "Difusion/Configuration%d.dat", d);
FILE *f = fopen(filename, "wt");
if(f!=NULL){
    for(int i =0; i<NPart;i++){
        fprintf(f, "%lf\t %lf\t %lf\n", r[i][0], r[i][1], r[i][2]);
    }
    fclose(f);
}
else
    printf("No se ha podido crear el fichero: %s", filename);
}

int InicializarFicheroTermodinamica(char *filename, double densi){
    FILE *f = fopen(filename, "w");
    if(f!=NULL){
        fprintf(f, "Numero de particulas: %d. \n", NPart);
        fprintf(f, "Densidad del sistema: %lf. \n", rho*densi);
        fprintf(f, "Temperatura del baño termico: %lf K. \n", Termostato*Temp);
        fprintf(f, "Tiempo de termalizacion: %lf ps.\n", TTermalizacion*t);
        fprintf(f, "Tiempo total de simulacion (en equilibrio): %lf ps.\n \n", TTot*t);

        fprintf(f, "Tiempo [ps]\t E. Cinetica [kJ/mol]\t E. Potencial [kJ/mol]\t Temperatura [K]\t Presion [GPa]\n");
        fclose(f);
        return 1;
    }
    else{
        printf("No se pudo abrir el fichero: %s.\n", filename);
        return 0;
    }
}

double distancia(double x, double y, double z){
    return x*x+y*y+z*z;
}



double EnergiaPBC(double (&r)[NPart][3], double L){
    double UPot=0, d, d6, d12, cutoff2 = 0.999*L*L/4, cutoff6 = cutoff2*cutoff2*cutoff2, Uoff = -4.0*(1.0/(cutoff6) - 1.0/(cutoff6*cutoff6)), x, y, z;
    for(int i=0;i<NPart; i++){
        for(int j=i+1; j <NPart; j++){
            x=PBC(r[j][0]-r[i][0], L);
            y=PBC(r[j][1]-r[i][1], L);
            z=PBC(r[j][2]-r[i][2], L);
            d = distancia(x,y, z);
            if(d<cutoff2){
                d6 = d*d*d;
                d12=d6*d6;
                UPot=UPot-4.0*(1.0/d6 - 1.0/d12)-Uoff;
            }
        }
    }
    return UPot;
}

void Fuerza(double (&F)[NPart][3],const double (&r)[NPart][3], double  L){
double d2, d8, d14, R[3], cutoff2 = 0.999*L*L*1.0/4, module;
for(int i=0; i<NPart; i++){
    for(int j=0; j<3; j++){
        F[i][j]=0;
    }
}
for(int i=0; i<NPart; i++){
    for(int j=i+1; j<NPart; j++){
        for(int l=0; l<3; l++){
            R[l] = PBC(r[i][l]-r[j][l], L);
        }
        d2=distancia(R[0], R[1], R[2] );
        if(d2 < cutoff2){
            d8 = d2*d2*d2*d2;
            d14 = d8*d2*d2*d2;
            module = 48.0/d14 - 24.0/d8;
            for(int n=0; n<3; n++){
                F[i][n]=F[i][n] + module*R[n];
                F[j][n] = F[j][n] - module*R[n];
            }
        }
    }
}
}

double Kinetic(double (&v)[NPart][3]){
double K = 0;
for(int i=0; i<NPart; i++){
    for(int j=0; j<3; j++)
        K=K+0.5*v[i][j]*v[i][j];
}
return K;
}

void VVS(double (&r)[NPart][3], double (&v)[NPart][3], double (&F)[NPart][3], double L, int (&NVueltas)[NPart][3]){
    //Nuevas posiciones.
    for(int n=0; n<NPart; n++){
        for(int l=0; l<3; l++){
            v[n][l] = v[n][l]+F[n][l]*dt/2;

            if(r[n][l]>L/2){
                NVueltas[n][l]++;
                r[n][l] = PBC(r[n][l] + v[n][l]*dt, L);
            }
            else{
                if(r[n][l]<-L/2){
                    NVueltas[n][l]-=1;
                    r[n][l] = PBC(r[n][l] + v[n][l]*dt, L);
                }
            }
        }
    }

    Fuerza(F, r, L);

    //Nuevas velocidades
    for(int n=0; n<NPart; n++){
        for(int l=0; l<3; l++){
            v[n][l] =v[n][l]+  F[n][l]*dt/2;
        }
    }
}

void Bath(double (&v)[NPart][3], double G){
int n=0;
for(int i=0; i<NPart; i++){
    for(int j=0; j<3; j++){
        if(Random()<0.01){
            v[i][j] = sqrt(G)*DistrGauss();
        }
    }
}
}


void Difusion(int (&NVueltas)[NPart][3], double (&r)[NPart][3], char *filename, double tt, double L){
FILE *fr = fopen("Difusion/Configuration0.dat", "r");
if(fr != NULL){
    FILE *fw = fopen(filename, "a");
    if(fw !=NULL){
        double x, y, z, MSD=0;
        for(int i=0; i<NPart; i++){
            fscanf(fr, "%lf \t %lf\t %lf\n", &x, &y, &z);
            MSD += distancia(r[i][0]-x + L*NVueltas[i][0],r[i][1]-y + L*NVueltas[i][1],r[i][2]-z + + L*NVueltas[i][2]);
        }
        fprintf(fw, "%lf \t %lf \n", t*tt, MSD*sigma*sigma/NPart);
        fclose(fw);
    }
    else
        printf("No se ha podido abrir el fichero: %s.\n", filename);
    fclose(fr);
}
else
    printf("No se ha podido abrir el fichero: Difusion/Configuration0.dat");

}


void InicializarSistema(double (&r)[NPart][3], double (&v)[NPart][3], double L){
    double vmod = 10;
    for(int i=0; i< NPart; i++){
        if(Random()<0.5)
            vmod = -vmod;
        for(int j=0; j<3; j++){
            v[i][j]=vmod;
        }
    }
}

void Termalizacion(double (&r)[NPart][3], double (&v)[NPart][3], double (&F)[NPart][3], double L, double T){
int NVueltas[NPart][3];
for(int i=0; i<NTermalizacion; i++){
    VVS(r, v, F, L, NVueltas);
    Bath(v, T);
}
}



double Presion(double (&r)[NPart][3], double (&F)[NPart][3], double T, double L, double densi){
double P = 0;
for(int i=0; i< NPart; i++){
    for(int j=0; j< 3; j++){
        P=P+r[i][j]*F[i][j];
    }
}
return P*1.0/3/(L*L*L) + densi*T;
}

int Flag(double (&r)[NPart][3], double L){
    for(int n=0; n<NPart; n++){
        for(int j=0; j<3; j++){
            if(fabs(r[n][j])>L/2){
                printf("Ha habido un error en la simulacion.\n");
                return 1;
            }
        }
    }
    return 0;
}

void RadialDistribution(double (&r)[NPart][3], double L, char *filename){
FILE *f = fopen(filename, "a");
if(f!=NULL){
    for(int i=0; i<NPart-1; i++){
        for(int j=i+1; j<NPart; j++)
            fprintf(f,"%lf\n", sigma*sqrt(distancia(r[i][0]-r[j][0],r[i][1]-r[j][1], r[i][2]-r[j][2] )));
    }
    fclose(f);
}
else
    printf("No se pudo abrir el fichero: %s", filename);
}

void Simulacion(char *filename,char *filename2, char *filename3, double densi){

    FILE *fkin=fopen(filename, "a");
    if(fkin!=NULL ){
        double K, T,F[NPart][3], G, r[NPart][3], v[NPart][3], L=IniBCC(r, v, densi), T0 = 2.0/(3*NPart-3);
        int NVueltas[NPart][3];

        //Inicializo el valor de las fuerzas.
        Fuerza(F, r, L);

        //Desordeno el sistema para crear una configuracion liquida.
        printf("Desordenando el sistema...\n");
        Termalizacion(r, v,F, L, 100);

        //Termalizo el sistema para llevarlo al equilibrio.
        printf("Sistema desordenado.\nTermalizando...\n");
        Termalizacion(r, v, F, L, Termostato);

        printf("Sistema termalizado.\nSimulando...\n");
        GuardarConfi(r, 0);

        for(int i=0; i< NPart; i++){
            for(int j=0; j<3; j++)
                NVueltas[i][j]=0;
        }

        for(int i=0; i<MedidasTotal ; i++){
            if(i%Medida==0){
                K =  Kinetic(v);
                T=T0*K;
                Difusion(NVueltas, r, filename2, i*dt, L);
                RadialDistribution(r, L, filename3);
                fprintf(fkin, "%lf\t%lf\t%lf\t%lf\t%lf\n", t*i*dt, K*e/NPart, e*EnergiaPBC(r,L)/NPart, T*Temp, Pres*Presion(r, F, T, L, densi) );
            }
            VVS(r, v, F, L, NVueltas);
            Bath(v, Termostato);
    }
    fclose(fkin);
}

else
    printf("No se pudo abrir el fichero: %s.\n", filename);
}

int main(){
    ini_ran(time(NULL));
    char filename[50], filename2[50], filename3[50];
    double densi;
    for(int i=0; i<5; i++){
        densi = Densidad[i];
        sprintf(filename, "Thermodynamics/density = %lf.dat", densi);
        sprintf(filename2, "Difusion/density = %lf.dat", densi);
        sprintf(filename3,"Distribucion/density = %lf.dat", densi );
        if(InicializarFicheroTermodinamica(filename, densi)==0)
            return 0;
        Simulacion(filename, filename2, filename3, densi);
        printf("He terminado la simulacion con densidad %lf.\n", densi);
    }
    return 0;
}
