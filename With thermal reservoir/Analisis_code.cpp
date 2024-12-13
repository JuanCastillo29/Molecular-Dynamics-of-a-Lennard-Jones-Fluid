#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

using namespace std;

const int N[10] = { 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024}, NDatos = 131072, NLineas = 7, NColumnas = 5;
char FicheroAnalisis[50] = {"Thermodynamics/density = 0.400000.dat"};
char AnalisisPresion[50] = {"Analisis/Presion.dat"};
char AnalisisCinetica[50] = {"Analisis/Cinetica.dat"};
char AnalisisPotencial[50] = {"Analisis/Potencial.dat"};
char AnalisisTemperatura[50] = {"Analisis/Temperatura.dat"};

void Binding1(double *data, double &media, double &Varm, const char *X, int c){
    FILE *f = fopen(X, "r");
    if(f!=NULL){
        int n=0;
        double aux;
        Varm = 0;

        // Saltar las primeras N líneas
        for (int i = 0; i < NLineas; i++) {
            if (fscanf(f, "%*[^\n]\n") == EOF) {  // Leer y descartar la línea
                printf("El archivo tiene menos de %d líneas.\n", N);
                fclose(f);
                return;
            }
        }

        // Escaneo del fichero y cálculo de la media
        for (int i = 0; i < NDatos-1; i++) {
            for (int j = 1; j <= NColumnas; j++) {
                if (fscanf(f, "%lf", &aux) != 1) { // Validar la lectura
                    if (feof(f)) { // Fin del archivo alcanzado
                    printf("Fin del archivo alcanzado antes de completar los datos.\n");
                } else { // Error de lectura
                    printf("Error al leer datos en la línea %d, columna %d.\n", NLineas + i + 1, j);
                    }
                fclose(f);
                return;
                }

                if (j == c) { // Solo procesar la columna de interés
                    *(data + i) = aux;
                    media += aux / NDatos;
                }
            }
        }

        //Calculo de la varianza.
        for(int i=0; i< NDatos; i++){
            Varm += (*(data+i) - media)*(*(data+i) - media);
        }
        Varm = sqrt(Varm)/(NDatos);

        fclose(f);
    }
    else
        printf("No se pudo abrir el archivo: %s.\n", X);
}

double BindingN(double *data, int Ldata, const double media){
    double Varianza=0, aux;
    for(int i=0; i<Ldata; i++){
        *(data + i) = (*(data + 2*i) + *(data + 2*i +1))/2;
        aux= *(data+i) - media;
        Varianza += aux*aux;
    }
    return sqrt(Varianza/(Ldata)/(Ldata-1));
}

void Binding(char (&X)[50], char (&AnalisisX)[50], int c){
    FILE *f = fopen(AnalisisX, "w");
    if(f!=NULL){
        double data[NDatos], media=0, Varm = 0;
        Binding1(data, media, Varm, X, c);
        fprintf(f, "La media es: %lf.\n", media);
        fprintf(f, "%d\t %lf\n", 1, Varm);
        for(int n=0; n<sizeof(N)/sizeof(const int); n++){
            fprintf(f, "%d \t %lf\n", N[n], BindingN(data, NDatos/N[n], media));
        }
        fclose(f);
    }
    else
        printf("No se pudo abrir el fichero: %s.\n", AnalisisX);
}

int main(){
    Binding(FicheroAnalisis, AnalisisCinetica, 2);
    Binding(FicheroAnalisis, AnalisisPotencial, 3);
    Binding(FicheroAnalisis, AnalisisTemperatura, 4);
    Binding(FicheroAnalisis, AnalisisPresion, 5);
    return 0;
}
