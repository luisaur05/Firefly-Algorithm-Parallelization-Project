#include <iostream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <time.h>
#include <algorithm>
#include <omp.h>

using namespace std;

// Definición de la función objetivo
double f_0(double x, double y) {
    return (abs(x) + abs(y)) * exp(-0.0625 * (x * x + y * y));
}

//Funciones 
double f_1(double x, double y){
    double a = sin(5*x);
    double b = sin(5*y);
    double fun;

    fun = x*x + y*y + 3*sqrt( a*a + b*b ) + 0.1;
    return (fun);
}

double f_2(double x, double y){
    double a = (4 - 2.1*x*x + (1/2)*x*x*x*x)*x*x;
    double b = x*y;
    double c = 4*(y*y-1)*y*y;
    double resultado;

    resultado=  a + b + c;
    return (resultado);
}

// Encuentra los límites de las coordenadas dentro del rango especificado
void findrange(vector<double>& xn, vector<double>& yn, const double range[]) {
    int length = yn.size();

    #pragma omp parallel for
    for (int i = 0; i < length; i++) {
        if (xn[i] <= range[0])
            xn[i] = range[0];
        if (xn[i] >= range[1])
            xn[i] = range[1];
        if (yn[i] <= range[2])
            yn[i] = range[2];
        if (yn[i] >= range[3])
            yn[i] = range[3];
    }
}

// Inicializa las luciérnagas
void init_ffa(int n, double range[], vector<double>& xn, vector<double>& yn, vector<double>& Lightn) {
    double xrange = range[1] - range[0];
    double yrange = range[3] - range[2];
    xn.resize(n);
    yn.resize(n);
    Lightn.resize(n);

    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        xn[i] = (rand() / (double)RAND_MAX) * xrange + range[0];
        yn[i] = (rand() / (double)RAND_MAX) * yrange + range[2];
        Lightn[i] = 0.0;
    }
}

// Mueve las luciérnagas
void ffa_move(vector<double>& xn, vector<double>& yn, vector<double>& Lightn, const vector<double>& xo, const vector<double>& yo, const vector<double>& Lighto, double alpha, double gamma, const double range[]) {
    int ni = yn.size();
    int nj = yo.size();

    #pragma omp parallel for collapse(2)
    for (int i = 0; i < ni; i++) {
        for (int j = 0; j < nj; j++) {
            double r = sqrt(pow(xn[i] - xo[j], 2) + pow(yn[i] - yo[j], 2));
            if (Lightn[i] < Lighto[j]) {
                double beta0 = 1;
                double beta = beta0 * exp(-gamma * pow(r, 2));
                xn[i] = xn[i] * (1 - beta) + xo[j] * beta + alpha * ((double)rand() / RAND_MAX - 0.5);
                yn[i] = yn[i] * (1 - beta) + yo[j] * beta + alpha * ((double)rand() / RAND_MAX - 0.5);
            }
        }
    }

    findrange(xn, yn, range);
}

int main() {
    int n = 6000;
    int MaxGeneration = 100;
    double range[] = { -5, 5, -5, 5 };
    double alpha = 0.2;
    double gamma = 1.0;

    // Comenzamos a contar el tiempo
    double start = clock();
    srand(time(NULL));

    // Los valores de la cuadrícula se utilizan solo para mostrarlos
    int Ndiv = 100;
    double dx = (range[1] - range[0]) / Ndiv;
    double dy = (range[3] - range[2]) / Ndiv;
    vector<vector<double>> x(Ndiv + 1, vector<double>(Ndiv + 1));
    vector<vector<double>> y(Ndiv + 1, vector<double>(Ndiv + 1));
    vector<vector<double>> z(Ndiv + 1, vector<double>(Ndiv + 1));

    #pragma omp parallel for collapse(2)
    for (int i = 0; i <= Ndiv; i++) {
        for (int j = 0; j <= Ndiv; j++) {
            x[i][j] = range[0] + i * dx;
            y[i][j] = range[2] + j * dy;
            z[i][j] = f_0(x[i][j], y[i][j]); //Aqui tomamos la funcion objetivo que deseemos
        }
    }

    // Generar las ubicaciones iniciales de las luciérnagas
    vector<double> xn, yn, Lightn;
    init_ffa(n, range, xn, yn, Lightn);

    for (int i = 0; i < MaxGeneration; i++) {

        vector<double> zn(n);
        #pragma omp parallel for
        for (int j = 0; j < n; j++) {
            zn[j] = f(xn[j], yn[j]);
        }

        vector<double> Lighto = Lightn;
        vector<double> xo = xn;
        vector<double> yo = yn;

        vector<int> Index(n);
        #pragma omp parallel for
        for (int j = 0; j < n; j++) {
            Index[j] = j;
        }

        // Ordenar los índices según los valores de zn[a] < zn[b]
        sort(Index.begin(), Index.end(), [&](int a, int b) {
            return zn[a] < zn[b];
        });

        #pragma omp parallel for
        for (int h = 0; h < n; h++) { // Acomodar el resto usando sus copias y el orden de índices
            int j = Index[h];
            Lightn[h] = zn[j];
            xn[h] = xo[j];
            yn[h] = yo[j];
            Lighto[h] = Lightn[j];

            xo[h] = xn[h];
            yo[h] = yn[h];
            Lighto[h] = Lightn[h];
        }

        ffa_move(xn, yn, Lightn, xo, yo, Lighto, alpha, gamma, range);
    }

    // Detenemos el tiempo
    double end = clock();
    double total = (end - start) / CLOCKS_PER_SEC;

    // Imprimimos los resultados obtenidos
    for (int j = 0; j < n; j++) {
        cout << "La luciérnaga " << j << " (" << xn[j] << "," << yn[j] << "),\t tiene una intensidad " << Lightn[j] << endl;
    }
    cout << endl;
    cout << "Tiempo = " << total << endl;

    return 0;
}
