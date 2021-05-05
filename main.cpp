/*
 *  Equação de Burger
 *
 *        U_t + U*U_x = v U_xx     , U=U(x,t),  a < x < b,     t>0
 *
 *        U_t = DU/dt   , U_x = DU/dx   e U_xx = D(DU/dx)/dx
 * */

#include <iostream>
#include <cmath>

using namespace std;

typedef double (*FPtr)(double, double);

double alfa = 5.0;
double beta = 4.0;
double v = 0.05;

double x0 = 0.0;
double x1 = 1.0;
double t0 = 0.0;
double t1 = 1.0;

int M = 10;
int N = 2;
double h;
double k;
double sigma;
FPtr Uptr, EsqPtr, DirPtr;

bool debug = true;
double *w1 = nullptr;
double *wAux = nullptr;
double *F = nullptr;
double *S1 = nullptr;
double *S2 = nullptr;
double **w = nullptr;
double **DF = nullptr;
double **DF1 = nullptr;
double **DF2 = nullptr;

double DirFuncaoInicial(double x, double t);

double EsqFuncaoInicial(double x, double t);

double U(double x, double t);

void imprime(char *nome, double x, double t, FPtr Uptr);

void NewMatriz(int nLin, int nCol, double **&M, double m);

void NewVetor(int nCol, double *&V, double v);

//void printArray(char ,double a[][N])
void imprimeMatriz(char *nome, int nlin, int ncol, double **&A);

//void printArray(char ,double a[][N])
void imprimeVetor(char *nome, int ncol, double *&v);

void inicializa(double **w, FPtr Uptr);

void copiaVetor(int ini, int fim, double *&vOrig, double *&vDest);

void copiaColMatrizToVetor(int col, int linIni, int linFim, double **&MOrig, double *&vDest);

void CalculaDF1(double **&DF1);

void CalculaDF2(double **&DF2, double *&w1);

void somaMatrizes(double **&S, double **&A, double **&B);

void escalar_x_matriz(double escalar, double **&A);

void Matriz_x_vetor(double * vR, double ** A, double * v);

void CalculaF(int col, double * F, double ** w, double * wAux);

void CondicaoDeDirichlet(int n, double ** DF, double * F, double * w1, FPtr DirPtr, FPtr EsqPtr);

void copiaVetorToColMatriz(int col, double * v, double ** A);

void ZeraMatriz(double ** A);

int main() {

    if (M > 1) {
        h = (x1 - x0) / M;
        k = (t1 - t0) / N;
    } else {
        h = M;
        k = N;
        M = (x1 - x0) / h;
        N = (t1 - t0) / k;
    }


    sigma = v * k / (h * h);

    Uptr = &U;
    EsqPtr = &EsqFuncaoInicial;
    DirPtr = &DirFuncaoInicial;
    if (debug) {
        imprime("Funcao U(x=0.5; t=1) =", 0.5, 1, Uptr);
    }

    NewMatriz(M + 1, M + 1, DF, 0);
    NewMatriz(M + 1, M + 1, DF1, 0);
    NewMatriz(M + 1, M + 1, DF2, 0);
    NewMatriz(M + 1, N + 1, w, 0);
    NewVetor(M + 1, F, 0);
    NewVetor(M + 1, S1, 1);
    NewVetor(M + 1, S2, 1);
    NewVetor(M + 1, w1, 0);
    NewVetor(M + 1, wAux, 0);


    if (debug) {
        imprimeMatriz("Matriz DF =", M + 1, M + 1, DF);
        imprimeMatriz("Matriz DF1 =", M + 1, M + 1, DF1);
        imprimeMatriz("Matriz DF2 =", M + 1, M + 1, DF2);
        imprimeMatriz("Matriz w =", M + 1, N + 1, w);
        imprimeVetor("Transposta do vetor F =", M + 1, F);
        imprimeVetor("Transposta do vetor w1 =", M + 1, w1);
    }

    inicializa(w, Uptr);
    copiaColMatrizToVetor(0, 0, M, w, w1);

    for (int n = 0; n < N; n++) {
        for (int it = 1; it <= 1; ++it) {
            ZeraMatriz(DF1);
            ZeraMatriz(DF2);
            if (debug) {
                imprimeMatriz("Matriz w =", M + 1, N + 1, w);
                imprimeVetor("Transposta do vetor w1 =", M + 1, w1);
            }

            CalculaDF1(DF1);
            CalculaDF2(DF2, w1);
            somaMatrizes(DF, DF1, DF2);

            if (debug) {
                imprimeMatriz("Matriz DF1 =", M + 1, M + 1, DF1);
                imprimeMatriz("Matriz DF2 =", M + 1, M + 1, DF2);
                imprimeMatriz("Matriz DF =", M + 1, M + 1, DF);
            }

            escalar_x_matriz(0.5, DF2);
            somaMatrizes(DF1, DF1, DF2);
            Matriz_x_vetor(wAux, DF1, w1);
            CalculaF(n, F, w, wAux);
            CondicaoDeDirichlet(n, DF, F, w1, DirPtr, EsqPtr);

            if (debug) {
                imprimeMatriz("Matriz DF2 =", M + 1, M + 1, DF2);
                imprimeMatriz("Matriz DF1 =", M + 1, M + 1, DF1);
                imprimeVetor("Transposta do vetor F =", M + 1, F);
            }

        }
        //w(:,n+1)=w1
        copiaVetorToColMatriz(n + 1, w1, w);
        if (debug) {
            imprimeMatriz("Matriz w =", M + 1, N + 1, w);

        }

    }
    return 0;
}

void inicializa(double **w, FPtr Uptr) {
    for (int i = 0; i <= M; i++)
        w[i][0] = Uptr(x0 + i * h, 0);

}


/*
 *  Entrada:
 *  a1 = n
 *  a2 = v
 *  a3 = x
 */

double DirFuncaoInicial(double x, double t) {
    return 0.0;
}

double EsqFuncaoInicial(double x, double t) {
    return 0.0;
}

double U(double x, double t) {

    return (2 * v * beta * M_PI * sin(M_PI * x) * exp(-v * t * M_PI * M_PI))
           / (alfa + beta * cos(M_PI * x) * exp(-v * t * M_PI * M_PI));
}

void imprime(char *nome, double x, double t, FPtr Uptr) {

    printf("\n  %s \n", nome);
    printf(" %10.6f ", Uptr(x, t));
}

void NewMatriz(int nLin, int nCol, double **&M, double m) {
    if (M == nullptr) {
        M = new double *[nLin];
        for (int i = 0; i < nLin; i++)
            M[i] = new double[nCol];

        for (int i = 0; i < nLin; i++)
            for (int j = 0; j < nCol; j++)
                M[i][j] = m;
    } else
        cout << "Atencao !! Matriz jah foi criada";
}

void NewVetor(int nCol, double *&V, double v) {
    if (V == nullptr) {
        V = new double[nCol];
        for (int i = 0; i < nCol; i++)
            V[i] = v;
    } else
        cout << "Atencao !! Vetor jah foi criado";
}

//void printArray(char ,double a[][N]) {
void imprimeMatriz(char *nome, int nlin, int ncol, double **&A) {
    printf("\n  %s \n", nome);
    // loop through array's rows
    for (int i = 0; i < nlin; ++i) {
        // loop through columns of current row
        for (int j = 0; j < ncol; ++j)
            printf(" %7.4f", A[i][j]);
        cout << endl; // start new line of output
    } // end outer for
} // end function printArray

//void printArray(char ,double a[][N]) {
void imprimeVetor(char *nome, int ncol, double *&v) {
    printf("\n  %s \n", nome);
    // loop through columns of current row
    for (int j = 0; j < ncol; ++j)
        printf(" %6.4f", v[j]);
    cout << endl; // start new line of output
} // end function printArray

void copiaVetor(int ini, int fim, double *&vOrig, double *&vDest) {
    for (int i = ini; i <= fim; i++)
        vDest[i] = vOrig[i];
}

void copiaColMatrizToVetor(int col, int linIni, int linFim, double **&MOrig, double *&vDest) {
    for (int i = linIni; i <= linFim; i++)
        vDest[i] = MOrig[i][col];
}

void copiaVetorToColMatriz(int col, double * v, double ** A) {
    for (int i = 0; i <= M; i++)
        A[i][col] = v[i];
}


void CalculaDF1(double **&DF1) {
    DF1[0][0] = 1; // Condição de Dirichlet
    DF1[M][M] = 1; // Condição de Dirichlet
    DF1[0][1] = -sigma;
    DF1[M][M - 1] = -sigma;
    for (int i = 1; i < M; i++) {
        DF1[i][i] = 1 + 2 * sigma;
        DF1[i][i + 1] = -sigma;
        DF1[i][i - 1] = -sigma;
    }
}

void CalculaDF2(double **&DF2, double *&w1) {

    for (int i = 1; i < M; i++) {
        DF2[i][i] = (w1[i + 1] - w1[i - 1]) * k / (2 * h);
        DF2[i][i + 1] = w1[i] * k / (2 * h);
        DF2[i][i - 1] = -w1[i] * k / (2 * h);
    }
}

void somaMatrizes(double **&S, double **&A, double **&B) {
    for (int i = 0; i <= M; i++)
        for (int j = 0; j <= M; j++)
            S[i][j] = A[i][j] + B[i][j];

}

void escalar_x_matriz(double escalar, double **&A) {
    for (int i = 0; i <= M; i++)
        for (int j = 0; j <= M; j++)
            A[i][j] = A[i][j] * escalar;
}

void Matriz_x_vetor(double *vR, double **A, double *v) {
    double soma;
    for (int i = 0; i <= M; i++) {
        soma = 0.0;
        for (int j = 0; j <= M; j++) {
            soma = soma + A[i][j] * v[j];
        }
        vR[i] = soma;
    }
}

void CalculaF(int j, double *F, double **w, double *wAux) {
    for (int i = 0; i <= M; ++i)
        F[i] = -w[i][j] + wAux[i];

}

void CondicaoDeDirichlet(int n, double **DF, double *F, double *w1, FPtr DirPtr, FPtr EsqPtr) {
    for (int j = 0; j <= M; j++) {
        DF[0][j] = 0.0;
        DF[M][j] = 0.0;
    }
    DF[0][0] = 1.0;
    DF[M][M] = 0.0;

    F[0] = w1[0] - EsqPtr(x0, n * k);
    F[M] = w1[M] - DirPtr(x1, n * k);
}

void ZeraMatriz(double ** A){
    for (int i = 0; i <=M ; ++i) {
        for (int j = 0; j <=M ; ++j) {
            A[i][j]=0.0;

        }

    }
}