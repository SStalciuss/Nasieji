#include <iostream>
#include <time.h>
#include <sys/time.h>
#include <cmath>
#include <mpi.h>

using namespace std;

//===== Globalus kintamieji ===================================================

int num_points = 5000;      		// Tasku skaicius (max 50000). Didinant ilgeja matricos skaiciavimo ir sprendinio paieskos laikas
int num_variables = 3;      		// Tasku, kuriuos reikia rasti, skaicius
int num_possible_locations = 57;	// Galimu vietu naujiems objektams skaicius (pirmi taskai is masyvo)

double **points;            		// Masyvas taskams saugoti
double **distance_matrix;   		// Masyvas atstumu matricai saugoti

//===== Funkciju antrastes ====================================================

double get_time();                                          // Funkcija laiko matavimui
void load_data();                                           // Duomenu ikelimo is failo funkcija
double Haversine_distance(double, double, double, double);  // Atstumo skaiciavimo pagal Haversino formule funkcija
double distance_from_matrix(int, int);                      // Atstumo paemimo is atstumu matricos funkcija
bool next_solution(int*); 									// Sekancio sprendinio generavimo funkcija
double evaluate_solution(int*);                             // Funkcija sprendinio tikslo funkcijos reiksmei ivertinti

//=============================================================================

int n_choose_k(int n, int k){
    if (k == 0) return 1;
    return (n * n_choose_k(n - 1, k - 1)) / k;
}

int find_row_for_element_count(int elements) {
    return static_cast<int>(std::floor((-1 + std::sqrt(1 + 8.0 * elements)) / 2.0));
}

int main(int argc, char *argv[]) {
	
	MPI_Init(&argc, &argv);

    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    srand(time(NULL));

	double t_0 = get_time();    // Programos vykdymo pradzios laiko fiksavimas
		
    load_data();                // Duomenu ikelimas is failo
	
	double t_1 = get_time();    // Duomenu ikelimo pabaigos laiko fiksavimas
	
    //-------------------------------------------------------------------------
	// Skaiciuojam atstumu matrica
    // Matrica yra "trikampe", nes atstumai nuo A iki B ir nuo B iki A yra lygus
    //-------------------------------------------------------------------------

    int all_matrix_elements = num_points * (num_points + 1) / 2;
    int element_per_proc = all_matrix_elements / world_size;

    int start = world_rank * element_per_proc;
    int end = start + element_per_proc;

    int start_row = find_row_for_element_count(start);
    int end_row = find_row_for_element_count(end);

    // std::cout << "Processor " << world_rank << " handles rows from " << start_row << " to " << end_row << std::endl;

    int num_elements = 0;
    for (int i = start_row; i < end_row; i++) {
        num_elements += i + 1;
    }

    double* local_distances = new double[num_elements];
    int local_index = 0;
    for (int i = start_row; i < end_row; i++) {
        for (int j=0; j<=i; j++) {
			local_distances[local_index++] = Haversine_distance(points[i][0], points[i][1], points[j][0], points[j][1]);
		}
    }

    int* recvcounts = new int[world_size];
    int* displacements = new int[world_size];

    MPI_Allgather(&num_elements, 1, MPI_INT, recvcounts, 1, MPI_INT, MPI_COMM_WORLD);

    displacements[0] = 0;
    for (int i = 1; i < world_size; i++) {
        displacements[i] = displacements[i - 1] + recvcounts[i - 1];
    }

    double* gathered_matrix = new double[(num_points * (num_points + 1)) / 2];
    MPI_Allgatherv(local_distances, num_elements, MPI_DOUBLE, gathered_matrix, recvcounts, displacements, MPI_DOUBLE, MPI_COMM_WORLD);

    int flat_index = 0;
    distance_matrix = new double*[num_points];
    for (int i=0; i<num_points; i++) {
		distance_matrix[i] = new double[i+1];
		for (int j=0; j<=i; j++) {
			distance_matrix[i][j] = gathered_matrix[flat_index++];
		}
	}

	double t_2 = get_time();    // Matricos skaiciavimo pabaigos laiko fiksavimas

    //-------------------------------------------------------------------------
	// Geriausio sprendinio pasieska pilno perrinkimo algoritmu
    // (angl. Complete Enumeration)
    //-------------------------------------------------------------------------
	
    int *solution = new int[num_variables];       // Masyvas atsitiktinai sugeneruotam sprendiniui saugoti
    int *best_solution = new int[num_variables];  // Masyvas geriausiam rastam sprendiniui saugoti
	double f_solution, f_best_solution = 1e10;     // Atsitiktinio ir geriausio rasto sprendiniu tikslo funkciju reiksmes

	for (int i=0; i<num_variables; i++) solution[i] = i;	// Pradinis sprendinys [0, 1, 2]

    int total_solution_count = n_choose_k(num_possible_locations, num_variables);
    int partition = total_solution_count / world_size;
    int startIndex = world_rank * partition;
    int items_to_proccess = (world_rank == world_size - 1) ? total_solution_count - startIndex : partition;

    for (int i = 0; i < startIndex; i++){
        next_solution(solution);
    }

	int processed = 0;
    while (items_to_proccess > 0 && next_solution(solution)) {
        f_solution = evaluate_solution(solution);
        if (f_solution < f_best_solution) {
            f_best_solution = f_solution;
            for (int j = 0; j < num_variables; j++) {
                best_solution[j] = solution[j];
            }
        }
        items_to_proccess--;
    }

    double f_global_best_solution;

    MPI_Reduce(&f_best_solution, &f_global_best_solution, 1, MPI_DOUBLE_INT, MPI_MINLOC, 0, MPI_COMM_WORLD);
    MPI_Bcast(&f_global_best_solution, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // if (f_best_solution == f_global_best_solution) {
        cout << "Geriausias rastas sprendinys (tasku indeksai duomenu masyve): ";
        for (int i=0; i<num_variables; i++) cout << best_solution[i] << "\t";
        cout << endl;
        cout << "Geriausia reiksme: " << f_best_solution << endl;
    // }

	double t_3 = get_time();    // Sprendinio paieskos pabaigos laiko fiksavimas

    //-------------------------------------------------------------------------
	// Rezultatų spausdinimas
    //-------------------------------------------------------------------------
        
    if(world_rank == 0){
        cout << "Duomenu ikelimo laikas: " << t_1 - t_0 << " s." << endl;
        cout << "Atstumu matricos skaiciavimo laikas: " << t_2 - t_1 << " s." << endl;
        cout << "Sprendinio paieskos laikas: " << t_3 - t_2 << " s." << endl;	
    }
    MPI_Finalize();
}

//=============================================================================
// Papildomos funkcijos. SIU FUNKCIJU LYGIAGRETINTI NEREIKIA.
//=============================================================================

double get_time() {
   struct timeval laikas;
   gettimeofday(&laikas, NULL);
   double rez = (double)laikas.tv_sec+(double)laikas.tv_usec/1000000;
   return rez;
}

//=============================================================================

void load_data() {
	
	//----- Load demand points ------------------------------------------------
	FILE *f;
	f = fopen("lab_data.dat", "r");
	points = new double*[num_points];
	for (int i=0; i<num_points; i++) {
		points[i] = new double[2];
		fscanf(f, "%lf%lf", &points[i][0], &points[i][1]);
	}
	fclose(f);
}

//=============================================================================

double Haversine_distance(double lat1, double lon1, double lat2, double lon2) {
	double dlat = fabs(lat1 - lat2);
	double dlon = fabs(lon1 - lon2);
	double aa = pow((sin((double)dlat/(double)2*0.01745)),2) + cos(lat1*0.01745) *
               cos(lat2*0.01745) * pow((sin((double)dlon/(double)2*0.01745)),2);
	double c = 2 * atan2(sqrt(aa), sqrt(1-aa));
	double d = 6371 * c; 
	return d;
}

//=============================================================================

double distance_from_matrix(int i, int j) {
	if (i >= j)	return distance_matrix[i][j];
	else return distance_matrix[j][i];
}

//=============================================================================

bool next_solution(int *solution) {
    int k = num_variables;
	int n = num_possible_locations;

	int i = k - 1;
    while (i >= 0 && solution[i] == n - k + i) i--;
    if (i < 0) return false;
    solution[i]++;
    for (int j = i + 1; j < k; j++) {
        solution[j] = solution[j - 1] + 1;
    }
    return true;
}

//=============================================================================

double evaluate_solution(int *solution) {
	double distance, min_distance, total_distance = 0;	
	for (int i=0; i<num_points; i++) {
        min_distance = 1e10;
        for (int j=0; j<num_variables; j++) {
		    distance = distance_from_matrix(i, solution[j]);
            if (distance < min_distance) {
                min_distance = distance;
            }
        }
        total_distance += min_distance;		    
    }
	return total_distance;
}

//=============================================================================