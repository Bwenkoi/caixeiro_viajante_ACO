/******************************************************************************
 Pedro Paes Siniscalchi - d2021101303@unifei.edu.br
 Breno de Oliveira Renó - brenooliveirareno@unifei.edu.br
 Maurício Andre de Almeida - d2021101420@unifei.edu.br

 ACO - Otimização por colônia de formigas aplicado ao problema do caixeiro viajante

 // FONTE: https://github.com/diogo-fernan/aco


*******************************************************************************/



#include <cstdio>
#include <iostream>
#include <cstdlib>

#include <cmath>
#include <limits>
#include <climits>


using namespace std;


class Randoms {

private:
    long xpto;
public:
    // Generator seed.
    Randoms (long x) {xpto = -x;}
    // Returns a random Gaussian number.
    double Normal (double avg, double sigma)
    {
        return (avg+sigma*gaussdev(&xpto)) ;
    }
    // Returns a uniform random number between 0 and 1.
    double Uniforme()
    {
        return ran1(&xpto);
    }
    // Returns a random number between -m and m.
    double sorte(int m)
    {
        return (1.0*rand())/(1.0*RAND_MAX)*2.0*m-m;
    }

/*
	Taken from Numerical Recipes in C, Chapter 7.
*/

#define   IA 16807
#define   IM 2147483647
#define   AM (1.0/IM)
#define   IQ 127773
#define   IR 2836
#define   NTAB 32
#define   NDIV (1+(IM-1)/NTAB)
#define   EPS 1.2e-7
#define   RNMX (1.0-EPS)

    float ran1(long *idum)
/*
"Minimal" random number generator of Park and Miller with Bays-Durham shuffle and added
safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint values). Call with idum a negative integer to initialize; thereafter, do not alter idum between successive deviates in a sequence. RNMX should approximate the largest floating value that is less than 1.
*/
    {
        int j;
        long k;
        static long iy=0;
        static long iv[NTAB];
        float temp;
        if (*idum <= 0 || !iy) {           // Initialize.
            if (-(*idum) < 1) *idum=1;     // Be sure to prevent idum = 0.
            else *idum = -(*idum);
            for (j=NTAB+7;j>=0;j--) {      // Load the shuffle table (after 8 warm-ups).
                k=(*idum)/IQ;
                *idum=IA*(*idum-k*IQ)-IR*k;
                if (*idum < 0) *idum += IM;
                if (j < NTAB) iv[j] = *idum;
            }
            iy=iv[0];
        }
        k=(*idum)/IQ;                     // Start here when not initializing.
        *idum=IA*(*idum-k*IQ)-IR*k;       // Compute idum=(IA*idum) % IM without over-
        if (*idum < 0) *idum += IM;       // flows by Schrage's method.
        j=iy/NDIV;                        //  Will be in the range 0..NTAB-1.
        iy=iv[j];                         // Output previously stored value and refill the
        iv[j] = *idum;                    // shuffle table.
        if ((temp=AM*iy) > RNMX)
            return RNMX; 				   // Because users don't expect endpoint values.
        else
            return temp;
    }

    float gaussdev(long *idum)
// Returns a normally distributed deviate with zero mean and unit variance,
// using ran1(idum) as the source of uniform deviates.
    {
//    float ran1(long *idum);

        static int iset=0;
        static float gset;
        float fac,rsq,v1,v2;
        if (*idum < 0) iset=0;      //     Reinitialize.
        if (iset == 0) {            //     We don't have an extra deviate handy, so
            do {
                v1=2.0*ran1(idum)-1.0;    // pick two uniform numbers in the square ex-
                v2=2.0*ran1(idum)-1.0;    // tending from -1 to +1 in each direction,
                rsq=v1*v1+v2*v2;          // see if they are in the unit circle,

            } while (rsq >= 1.0 || rsq == 0.0);  // and if they are not, try again.
            fac=sqrt(-2.0*log(rsq)/rsq);
            // Now make the Box-Muller transformation to get two normal deviates.
            // Return one and save the other for next time.
            gset=v1*fac;
            iset=1;                 //  Set flag.
            return v2*fac;
        } else {                    //   We have an extra deviate handy,
            iset=0;                 //   so unset the flag,
            return gset;            //   and return it.
        }
    }

};

class ACO {
public:
    ACO (int nAnts, int nCities,
         double alpha, double beta, double q, double ro, double taumax,
         int initCity);
    virtual ~ACO ();

    void init ();

    void connectCITIES (int cityi, int cityj, double distance);

    void printPHEROMONES ();
    void printGRAPH ();
    void printRESULTS ();

    void optimize (int ITERATIONS);

private:
    double distance (int cityi, int cityj);
    bool exists (int cityi, int cityc);
    bool vizited (int antk, int c);
    double PHI (int cityi, int cityj, int antk);

    double length (int antk);

    int city ();
    void route (int antk);
    int valid (int antk, int iteration);

    void updatePHEROMONES ();


    int NUMBEROFANTS, NUMBEROFCITIES, INITIALCITY;
    double ALPHA, BETA, Q, RO, TAUMAX;

    double BESTLENGTH;
    int *BESTROUTE;

    double **GRAPH;
    int **ROUTES;
    double **CITIES, **PHEROMONES, **DELTAPHEROMONES, **PROBS;

    Randoms *randoms;
};


ACO::ACO (int nAnts, int nCities,
          double alpha, double beta, double q, double ro, double taumax,
          int initCity) {
    NUMBEROFANTS 	= nAnts;
    NUMBEROFCITIES 	= nCities;
    ALPHA 			= alpha;
    BETA 			= beta;
    Q 				= q;
    RO 				= ro;
    TAUMAX 			= taumax;
    INITIALCITY		= initCity;

    randoms = new Randoms (21);
}
ACO::~ACO () {
    for(int i=0; i<NUMBEROFCITIES; i++) {
        delete [] GRAPH[i];
        delete [] CITIES[i];
        delete [] PHEROMONES[i];
        delete [] DELTAPHEROMONES[i];
        if(i < NUMBEROFCITIES - 1) {
            delete [] PROBS[i];
        }
    }
    delete [] GRAPH;
    delete [] CITIES;
    delete [] PHEROMONES;
    delete [] DELTAPHEROMONES;
    delete [] PROBS;
}

void ACO::init () {
    GRAPH 			= new double*[NUMBEROFCITIES];
    CITIES 			= new double*[NUMBEROFCITIES];
    PHEROMONES 		= new double*[NUMBEROFCITIES];
    DELTAPHEROMONES = new double*[NUMBEROFCITIES];
    PROBS 			= new double*[NUMBEROFCITIES-1];
    for(int i=0; i<NUMBEROFCITIES; i++) {
        GRAPH[i] 			= new double[NUMBEROFCITIES];
        CITIES[i] 			= new double[2];
        PHEROMONES[i] 		= new double[NUMBEROFCITIES];
        DELTAPHEROMONES[i] 	= new double[NUMBEROFCITIES];
        PROBS[i] 			= new double[2];
        for (int j=0; j<2; j++) {
            CITIES[i][j] = -1.0;
            PROBS[i][j]  = -1.0;
        }
        for (int j=0; j<NUMBEROFCITIES; j++) {
            GRAPH[i][j] 			= 0;
            PHEROMONES[i][j] 		= 0.0;
            DELTAPHEROMONES[i][j] 	= 0.0;
        }
    }

    ROUTES = new int*[NUMBEROFANTS];
    for (int i=0; i<NUMBEROFANTS; i++) {
        ROUTES[i] = new int[NUMBEROFCITIES];
        for (int j=0; j<NUMBEROFCITIES; j++) {
            ROUTES[i][j] = -1;
        }
    }

    BESTLENGTH = (double) INT_MAX;
    BESTROUTE  = new int[NUMBEROFCITIES];
    for (int i=0; i<NUMBEROFCITIES; i++) {
        BESTROUTE[i] = -1;
    }
}


void ACO::connectCITIES (int cityi, int cityj, double distance) {
    GRAPH[cityi][cityj] = distance;
    PHEROMONES[cityi][cityj] = randoms -> Uniforme() * TAUMAX;
    GRAPH[cityj][cityi] = distance;
    PHEROMONES[cityj][cityi] = PHEROMONES[cityi][cityj];
}

void ACO::printPHEROMONES () {
    cout << " FEROMONIOS: " << endl;
    cout << "  | ";
    for (int i=0; i<NUMBEROFCITIES; i++) {
        printf("%5d   ", i);
    }
    cout << endl << "- | ";
    for (int i=0; i<NUMBEROFCITIES; i++) {
        cout << "--------";
    }
    cout << endl;
    for (int i=0; i<NUMBEROFCITIES; i++) {
        cout << i << " | ";
        for (int j=0; j<NUMBEROFCITIES; j++) {
            if (i == j) {
                printf ("%5s   ", "x");
                continue;
            }
            if (exists(i, j)) {
                printf ("%7.3f ", PHEROMONES[i][j]);
            }
            else {
                if(PHEROMONES[i][j] == 0.0) {
                    printf ("%5.0f   ", PHEROMONES[i][j]);
                }
                else {
                    printf ("%7.3f ", PHEROMONES[i][j]);
                }
            }
        }
        cout << endl;
    }
    cout << endl;
}


double ACO::distance (int cityi, int cityj) {
    return (GRAPH[cityi][cityj]);
    /*(double)
        sqrt (pow (CITIES[cityi][0] - CITIES[cityj][0], 2) +
               pow (CITIES[cityi][1] - CITIES[cityj][1], 2)); */
}
bool ACO::exists (int cityi, int cityc) {
    return (GRAPH[cityi][cityc] != 0);
}
bool ACO::vizited (int antk, int c) {
    for (int l=0; l<NUMBEROFCITIES; l++) {
        if (ROUTES[antk][l] == -1) {
            break;
        }
        if (ROUTES[antk][l] == c) {
            return true;
        }
    }
    return false;
}
double ACO::PHI (int cityi, int cityj, int antk) {
    double ETAij = (double) pow (1 / distance (cityi, cityj), BETA);
    double TAUij = (double) pow (PHEROMONES[cityi][cityj],   ALPHA);

    double sum = 0.0;
    for (int c=0; c<NUMBEROFCITIES; c++) {
        if (exists(cityi, c)) {
            if (!vizited(antk, c)) {
                double ETA = (double) pow (1 / distance (cityi, c), BETA);
                double TAU = (double) pow (PHEROMONES[cityi][c],   ALPHA);
                sum += ETA * TAU;
            }
        }
    }
    return (ETAij * TAUij) / sum;
}

double ACO::length (int antk) {
    double sum = 0.0;
    for (int j=0; j<NUMBEROFCITIES-1; j++) {
        sum += distance (ROUTES[antk][j], ROUTES[antk][j+1]);
    }
    return sum;
}

int ACO::city () {
    double xi = randoms -> Uniforme();
    int i = 0;
    double sum = PROBS[i][0];
    while (sum < xi) {
        i++;
        sum += PROBS[i][0];
    }
    return (int) PROBS[i][1];
}

void ACO::route (int antk) {
    ROUTES[antk][0] = INITIALCITY;
    for (int i=0; i<NUMBEROFCITIES-1; i++) {
        int cityi = ROUTES[antk][i];
        int count = 0;
        for (int c=0; c<NUMBEROFCITIES; c++) {
            if (cityi == c) {
                continue;
            }
            if (exists (cityi, c)) {
                if (!vizited (antk, c)) {
                    PROBS[count][0] = PHI (cityi, c, antk);
                    PROBS[count][1] = (double) c;
                    count++;
                }

            }
        }

        // deadlock
        if (0 == count) {
            return;
        }

        ROUTES[antk][i+1] = city();
    }
}
int ACO::valid (int antk, int iteration) {
    for(int i=0; i<NUMBEROFCITIES-1; i++) {
        int cityi = ROUTES[antk][i];
        int cityj = ROUTES[antk][i+1];
        if (cityi < 0 || cityj < 0) {
            return -1;
        }
        if (!exists(cityi, cityj)) {
            return -2;
        }
        for (int j=0; j<i-1; j++) {
            if (ROUTES[antk][i] == ROUTES[antk][j]) {
                return -3;
            }
        }
    }

    if (!exists (INITIALCITY, ROUTES[antk][NUMBEROFCITIES-1])) {
        return -4;
    }

    return 0;
}

void ACO::printGRAPH () {
    cout << " GRAFO: " << endl;
    cout << "  |     ";
    for( int i=0; i<NUMBEROFCITIES; i++) {
        cout << i << "     ";
    }
    cout << endl << "- | ";
    for (int i=0; i<NUMBEROFCITIES; i++) {
        cout << "------";
    }
    cout << endl;
    int count = 0;
    for (int i=0; i<NUMBEROFCITIES; i++) {
        cout << i << " | ";
        for (int j=0; j<NUMBEROFCITIES; j++) {
            if(i == j) {
                cout << "    x ";
            }
            else {
                printf( " %4.0f ", GRAPH[i][j]);
                //cout << GRAPH[i][j] << " ";
            }
            if (GRAPH[i][j] == 1) {
                count++;
            }
        }
        cout << endl;
    }
    cout << endl;
    cout << "Numero de conexoes: " << count << endl << endl;
}
void ACO::printRESULTS () {
    BESTLENGTH += distance (BESTROUTE[NUMBEROFCITIES-1], INITIALCITY);
    cout << " MELHOR ROTA:" << endl;
    for (int i=0; i<NUMBEROFCITIES; i++) {
        cout << BESTROUTE[i] << " ";
    }
    cout << endl << "Comprimento: " << BESTLENGTH << endl;


}

void ACO::updatePHEROMONES () {
    for (int k=0; k<NUMBEROFANTS; k++) {
        double rlength = length(k);
        for (int r=0; r<NUMBEROFCITIES-1; r++) {
            int cityi = ROUTES[k][r];
            int cityj = ROUTES[k][r+1];
            DELTAPHEROMONES[cityi][cityj] += Q / rlength;
            DELTAPHEROMONES[cityj][cityi] += Q / rlength;
        }
    }
    for (int i=0; i<NUMBEROFCITIES; i++) {
        for (int j=0; j<NUMBEROFCITIES; j++) {
            PHEROMONES[i][j] = (1 - RO) * PHEROMONES[i][j] + DELTAPHEROMONES[i][j];
            DELTAPHEROMONES[i][j] = 0.0;
        }
    }
}


void ACO::optimize (int ITERATIONS) {
    for (int iterations=1; iterations<=ITERATIONS; iterations++) {
        cout << flush;
        cout << "ITERACAO " << iterations << " INICIO!" << endl << endl;

        for (int k=0; k<NUMBEROFANTS; k++) {
            cout << " : formiga " << k << " foi liberada!" << endl;
            while (0 != valid(k, iterations)) {
                cout << "  :: liberando formiga " << k << " novamente!" << endl;
                for (int i=0; i<NUMBEROFCITIES; i++) {
                    ROUTES[k][i] = -1;
                }
                route(k);
            }

            for (int i=0; i<NUMBEROFCITIES; i++) {
                cout << ROUTES[k][i] << " ";
            }
            cout << endl;


            double rlength = length(k);
            cout << "  :: rota terminada! "<< endl;
            //: Comprimento: " << rlength << endl;

            if (rlength < BESTLENGTH) {
                BESTLENGTH = rlength;
                for (int i=0; i<NUMBEROFCITIES; i++) {
                    BESTROUTE[i] = ROUTES[k][i];
                }
            }
            cout << " : formiga " << k << " finalizou!" << endl;
        }

        cout << endl << "atualizando FEROMONIOS . . .";
        updatePHEROMONES ();
        cout << " feito!" << endl << endl;
        printPHEROMONES ();

        for (int i=0; i<NUMBEROFANTS; i++) {
            for (int j=0; j<NUMBEROFCITIES; j++) {
                ROUTES[i][j] = -1;
            }
        }

        cout << endl << "ITERACAO " << iterations << " COMPLETA!" << endl << endl;
    }
}


#define ITERATIONS		(int) 5

#define NUMBEROFANTS	(int) 4
#define NUMBEROFCITIES	(int) 8


// if (ALPHA == 0) { stochastic search & sub-optimal route }
#define ALPHA			(double) 0.5
// if (BETA  == 0) { sub-optimal route }
#define BETA			(double) 0.8
// Estimation of the suspected best route.
#define Q				(double) 80
// Pheromones evaporation.
#define RO				(double) 0.2
// Maximum pheromone random number.
#define TAUMAX			(int) 2

#define INITIALCITY		(int) 0

int main() {

    ACO *ANTS = new ACO (NUMBEROFANTS, NUMBEROFCITIES,
                         ALPHA, BETA, Q, RO, TAUMAX,
                         INITIALCITY);

    ANTS -> init();

    ANTS -> connectCITIES (0, 1, 42);
    ANTS -> connectCITIES (0, 2, 61);
    ANTS -> connectCITIES (0, 3, 30);
    ANTS -> connectCITIES (0, 4, 17);
    ANTS -> connectCITIES (0, 5, 82);
    ANTS -> connectCITIES (0, 6, 31);
    ANTS -> connectCITIES (0, 7, 11);

    ANTS -> connectCITIES (1, 2, 14);
    ANTS -> connectCITIES (1, 3, 87);
    ANTS -> connectCITIES (1, 4, 28);
    ANTS -> connectCITIES (1, 5, 70);
    ANTS -> connectCITIES (1, 6, 19);
    ANTS -> connectCITIES (1, 7, 33);

    ANTS -> connectCITIES (2, 3, 20);
    ANTS -> connectCITIES (2, 4, 81);
    ANTS -> connectCITIES (2, 5, 21);
    ANTS -> connectCITIES (2, 6,  8);
    ANTS -> connectCITIES (2, 7, 29);

    ANTS -> connectCITIES (3, 4, 34);
    ANTS -> connectCITIES (3, 5, 33);
    ANTS -> connectCITIES (3, 6, 91);
    ANTS -> connectCITIES (3, 7, 10);

    ANTS -> connectCITIES (4, 5, 41);
    ANTS -> connectCITIES (4, 6, 34);
    ANTS -> connectCITIES (4, 7, 82);

    ANTS -> connectCITIES (5, 6, 19);
    ANTS -> connectCITIES (5, 7, 32);

    ANTS -> connectCITIES (6, 7, 59);


    ANTS -> printGRAPH ();

    ANTS -> printPHEROMONES ();

    ANTS -> optimize (ITERATIONS);

    ANTS -> printRESULTS ();

    return 0;
}

