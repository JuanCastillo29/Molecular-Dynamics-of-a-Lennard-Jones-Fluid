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

#define densi 0.8
#define NPart 125

#define TTimeStep 10000000
#define dt 0.00001
#define TTermalizacion 100000

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
    //Se usa la fórmula para crear dos números independientes según una gaussiana
    return -sqrt(-2*log(d1))*cos(2*Pi*d2);
}


double IniSC(double (&r)[NPart][3]){
    int N = int(pow(NPart*1.0001, 1.0/3));
    double l = pow(NPart/densi, 1.0/3), a=l/N;
    for(int i=0; i<NPart; i++){
        r[i][2]= (i/(N*N));
        r[i][1]=(i%(N*N))/N;
        r[i][0] = (i)%N;
        for(int j=0; j<3; j++){
            r[i][j] = a*r[i][j];
        }
    }
return l;
}

double IniBCC(double (&r)[NPart][3]){
    int N = int(pow((NPart/2)*1.0001, 1.0/3));
    double l = pow((NPart/2)/densi, 1.0/3), a = l/N;
    for(int i=0; i<NPart/2; i++){
        r[i][2]= (i/(N*N));
        r[i][1]=(i%(N*N))/N;
        r[i][0] = (i)%N;
        for(int j=0; j<3; j++){
            r[i+NPart/2][j] = a*(r[i][j] + 0.5);
            r[i][j] = a*r[i][j];
        }
    }
return l;
}

void GuardarConfi(double (&r)[NPart][3], int d){
char filename[20];
sprintf(filename, "Trayectoria/Configuration%d.xyz", d);
FILE *f = fopen(filename, "wt");
if(f!=NULL){
    fprintf(f, "%d\n\n", NPart);
    for(int i =0; i<NPart;i++){
        fprintf(f, "%d\t%lf\t %lf\t %lf\n", i, r[i][0], r[i][1], r[i][2]);
    }
    fclose(f);
}
else
    printf("No se ha podido crear el fichero: %s", filename);
}

double distancia(double x, double y, double z){
    return x*x+y*y+z*z;
}

double Energia(double (&r)[NPart][3], double cutoff2){
    double UPot=0, d, d6, d12, cutoff6 = cutoff2*cutoff2*cutoff2, Uoff = -4.0*(1.0/(cutoff6) - 1.0/(cutoff6*cutoff6));
    for(int i=0;i<NPart; i++){
        for(int j=i+1; j <NPart; j++){
            d = distancia(r[j][0]-r[i][0],r[j][1]-r[i][1], r[j][2]-r[i][2] );
            if(d<cutoff2){
                d6 = d*d*d;
                d12=d6*d6;
                UPot=UPot-4.0*(1.0/d6 - 1.0/d12)-Uoff;
            }
        }

    }
    return UPot;
}

double PBC(double x, double L){
if(x>L/2)
    x=x-L;
if(x<-L/2)
    x=x+L;
return x;
}

double EnergiaPBC(double (&r)[NPart][3], double L){
    double UPot=0, d, d6, d12, cutoff2 = L*L/4, cutoff6 = cutoff2*cutoff2*cutoff2, Uoff = -4.0*(1.0/(cutoff6) - 1.0/(cutoff6*cutoff6)), x, y, z;
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
double d2, d8, d14, R[3], cutoff2 = L*L*1.0/4, module;
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

void CopyArray(double Conservar[NPart][3], double Borrar[NPart][3]){
    for(int i=0; i<NPart; i++){
        for(int j=0; j<3; j++){
            Borrar[i][j] = Conservar[i][j];
        }
    }
}

void Euler(double (&r)[NPart][3], double (&v)[NPart][3], double L){
FILE *fkin=fopen("Euler/Energia.txt", "wt");
if(fkin!=NULL ){
    double K, P, F[NPart][3];
    for(int i=0; i<TTimeStep; i++){

        GuardarConfi(r, i);
        K=Kinetic(v);
        P=Energia(r,L*L/4);
        fprintf(fkin, "%lf\t%lf\t%lf\t%lf\n", i*dt,K, P, P+K );

        Fuerza(F, r, L);

        //Nuevas velocidades
        for(int n=0; n<NPart; n++){
            for(int l=0; l<3; l++){
                r[n][l] = r[n][l] + v[n][l]*dt;
                v[n][l] =v[n][l]+  (F[n][l])*dt;
            }
        }
    }
}

else{
    printf("Remember that you have to create a directory called Euler.\n");
}
}

void Verlet(double (&r)[NPart][3], double (&r0)[NPart][3], double L){
FILE *fkin=fopen("Verlet/Energia.txt", "wt");
if(fkin!=NULL ){
    double K, P, v[NPart][3], aux[NPart][3], F[NPart][3];
    for(int i=0; i<TTimeStep; i++){

        CopyArray(r, aux);
        GuardarConfi(r, i);
        K=Kinetic(v);
        P=Energia(r,L*L/4);
        fprintf(fkin, "%lf\t%lf\t%lf\t%lf\n", i*dt,K, P, P+K );

        Fuerza(F, r, L);
        for(int n=0; n<NPart; n++){
            for(int l=0; l<3; l++){
                r[n][l] = PBC(2*aux[n][l] - r0[n][l] + F[n][l]*dt*dt, L);
                v[n][l]=(r[n][l]-r0[n][l])/(2*dt);
            }
        }
        CopyArray(aux, r0);

    }
}

else{
    printf("Remember that you have to create a directory called Verlet.\n");
}
}

void VVS(double (&r)[NPart][3], double (&v)[NPart][3], double (&F)[NPart][3], double L){
    //Nuevas posiciones.
    for(int n=0; n<NPart; n++){
        for(int l=0; l<3; l++){
            v[n][l] = v[n][l]+F[n][l]*dt/2;
            r[n][l] = PBC(r[n][l] + v[n][l]*dt, L);
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
        if(Random()<0.1){
            v[i][j] = sqrt(G)*DistrGauss();
        }
    }
}
}

void VerletVel(double (&r)[NPart][3], double (&v)[NPart][3], double L){
FILE *fkin=fopen("Thermodynamics.txt", "wt");
if(fkin!=NULL ){
    double K, P, F[NPart][3], G;
    Fuerza(F, r, L);
    G=3;
    for(int i=0; i<TTimeStep; i++){
        GuardarConfi(r, i);
        K=Kinetic(v);
        P=EnergiaPBC(r,L);
        fprintf(fkin, "%lf\t%lf\t%lf\t%lf\n", i*dt,K, P, 2*K/(3*NPart-3));
        VVS(r, v, F, L);
        //Bath(v, G);
    }
    fclose(fkin);
}

else{
    printf("Remember that you have to create a directory called VerletVel.\n");
}
}

void Termalizacion(double (&r)[NPart][3], double (&v)[NPart][3], double L){
double F[NPart][3], G=100;
Fuerza(F, r, L);
for(int i=0; i<TTermalizacion; i++){
    VVS(r, v, F, L);
    Bath(v, G);
}
GuardarConfi(r, 0);
}

void InicializarSistema( double (&v)[NPart][3]){
    double vmod = 10;
    for(int i=0; i< NPart; i++){
        if(Random()<0.5)
            vmod = -vmod;
        for(int j=0; j<3; j++){
            v[i][j]=vmod;
        }
    }
}

void guardarvel(char (&filename)[50],double (&v)[NPart][3]){
FILE *f = fopen(filename, "wt");
if(f!=NULL){
    for(int i=0; i< NPart;i++){
        for(int j=0; j< 3; j++)
            fprintf(f, "%lf\n", v[i][j]);
    }
}
}

double Presion(double (&r)[NPart][3], double (&F)[NPart][3], double T){
double P = 0;
for(int i=0; i< NPart; i++){
    for(int j=0; j< 3; j++){
        P=P+r[i][j]*F[i][j];
    }
}
return P*1.0/3 + densi*T;
}

void Cinetica(char *filename, int i, const double (&v)[NPart][3], const double E){
char mode[10];
if(i==0)
    sprintf(mode, "wt");
else
    sprintf(mode, "at");

FILE *f = fopen(filename, mode);
if( f!=NULL){
    double P[3] = {0, 0 , 0};
    for(int n=0; n<NPart; n++){
        for(int j=0; j<3; j++){
            P[j] += v[n][j];
        }
    }
    fprintf(f, "%lf\t %lf \t %lf \t %lf \t %lf\n", i*dt, P[0], P[1], P[2], E);
    fclose(f);
}
else
    printf("No se pudo abrir el fichero %s", filename);
}

void Simulacion(void){
char directoryname[150], filename[50], filename2[50];
sprintf(filename, "Verlet/dt=%lf", dt);
sprintf(filename2, "Verlet/dt=%lf", dt);
getcwd(directoryname, sizeof(directoryname));
strcat(directoryname, "/");
strcat(directoryname,filename);
mkdir(directoryname);

strcat(filename2, "/Kinetics.txt");
strcat(filename,"/Thermodynamics.txt" );
FILE *fkin=fopen(filename, "wt");
if(fkin!=NULL ){
    double K, P, T,F[NPart][3], G, r[NPart][3], v[NPart][3], L=IniSC(r);
    Termalizacion(r, v, L);
    sprintf(filename, "Verlet/dt=%lf/VelInicial.txt", dt);
    InicializarSistema(v);
    guardarvel(filename, v);
    Fuerza(F, r, L);
    //G=3;
    for(int i=0; i<TTimeStep; i++){
        //GuardarConfi(r, i);
        if(i%10000==0){
            K=Kinetic(v);
            P=EnergiaPBC(r,L);
            T=2*K/(3*NPart-3);
            fprintf(fkin, "%lf\t%lf\t%lf\t%lf\t%lf\n", i*dt,K, P, T, Presion(r, F, T) );
            Cinetica(filename2, i, v, K+P);
        }
        VVS(r, v, F, L);
        //Bath(v, G);
    }
      sprintf(filename, "Verlet/dt=%lf/VelFinal.txt", dt);
    guardarvel( filename, v);
    fclose(fkin);
}

else{
    printf("Remember that you have to create a directory called VerletVel.\n");
    fclose(fkin);
}
}

int main(){
    ini_ran(time(NULL));
    Simulacion();

return 0;
}
