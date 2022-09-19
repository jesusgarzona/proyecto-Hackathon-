#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <stdint.h>

/*
Haciendo pruebas, el programa me estaba dando un error de runtime cuando compilaba
con U y W al tiempo (los arrays originales). Era un error absolutamente absurdo que mencionaba
que no encontraba un archivo de Cygwin (y no estaba compilando con Cygwin).
Sin embargo, al compilar únicamente con U el error no ocurría. Reduje el tamaño de M y N a 100 y el
error dejó de ocurrir.

Conclusión: 500x500x8x2 bytes (el peso de ambos arrays) es aproximadamente 3.8 MB, por lo que el programa
estaba haciendo stack overflow (el stack es típicamente de 1-2 MB).

Solución: generar los arrays en el heap.
*/

#define M 500
#define N 500

//resolución de los heat_plate(s) en su formato gráfico
#define GRAPHIC_N 200

#define MAX_STORED_FRAMES 400

#define HEAT_MAX (100.0)
#define HEAT_MIN (0.0)

//aka each 50 simulation steps, a graphic frame is stored
#define GRAPHIC_PRECISION 50

typedef struct {
	double cell[M][N];
} heat_plate;

//pixelbuffer
//color is interpreted as 8 bit intensity. Intensity is relative to heat range
typedef struct {
	// access with [x][y]
	uint8_t cell[GRAPHIC_N][GRAPHIC_N];
} graphic_heat_plate;

/*
Graphic file format:
int16_t <n>: defines how many frames are in the file
then <n> consecutive graphic_heat_plate(s), defined
as their cell values in for(y) { for(x) } order.
*/

double cpu_time () {
	return ( double ) clock () / ( double ) CLOCKS_PER_SEC;
}

static inline int iclamp(int v, int vmin, int vmax) {
	if (v < vmin) return vmin;
	if (v > vmax) return vmax;
	return v;
}

static inline int imin(int a, int b) {
	if (a < b) return a;
	return b;
}

int main ( int argc, char *argv[] ) {	

	double ctime;
	double ctime1;
	double ctime2;
	double diff;
	double epsilon;
	int i;
	int iterations;
	int iterations_print;
	int j;
	double mean;
	char output_file[80];
	int success;
	heat_plate* u = malloc(sizeof(heat_plate));
	heat_plate* w = malloc(sizeof(heat_plate));
	graphic_heat_plate** graphic_history = malloc(sizeof(graphic_heat_plate*) * MAX_STORED_FRAMES );

	printf ( "\n" );
	printf ( "HEATED_PLATE\n" );
	printf ( "  C version\n" );
	printf ( "  A program to solve for the steady state temperature distribution\n" );
	printf ( "  over a rectangular plate.\n" );
	printf ( "\n" );
	printf ( "  Spatial grid of %d by %d points.\n", M, N );
	
	//Read EPSILON from the command line or the user.
	if ( argc < 2 ) 
	{
	printf ( "\n" );
	printf ( "  Enter EPSILON, the error tolerance:\n" );
	success = scanf ( "%lf", &epsilon );
	}
	else
	{
	success = sscanf ( argv[1], "%lf", &epsilon );
	}

	if ( success != 1 )
	{
	printf ( "\n" );
	printf ( "HEATED_PLATE\n" );
	printf ( "  Error reading in the value of EPSILON.\n");
	return 1;
	}

	printf ( "\n" );
	printf ( "  The iteration will be repeated until the change is <= %f\n", epsilon );
	diff = epsilon;
	
	//Read OUTPUT_FILE from the command line or the user.
	if ( argc < 3 ) 
	{
	printf ( "\n" );
	printf ( "  Enter OUTPUT_FILE, the name of the output file:\n" );
	success = scanf ( "%s", output_file );
	}
	else
	{
	success = sscanf ( argv[2], "%s", output_file );
	}

	if ( success != 1 )
	{
	printf ( "\n" );
	printf ( "HEATED_PLATE\n" );
	printf ( "  Error reading in the value of OUTPUT_FILE.\n");
	return 1;
	}

	printf ( "\n" );
	printf ( "  The steady state solution will be written to '%s'.\n", output_file );
	
	//Set the boundary values, which don't change. 
	for ( i = 1; i < M - 1; i++ ) { w->cell[i][0]   = HEAT_MAX; }
	for ( i = 1; i < M - 1; i++ ) { w->cell[i][N-1] = HEAT_MAX; }
	for ( j = 0; j < N; j++ ) 	  { w->cell[M-1][j] = HEAT_MAX; }
	for ( j = 0; j < N; j++ ) 	  { w->cell[0][j]   = HEAT_MIN; }
	
	/*
	//Average the boundary values, to come up with a reasonable
	//initial value for the interior.
	mean = 0.0;
	for ( i = 1; i < M - 1; i++ ) { mean = mean + w->cell[i][0]; }
	for ( i = 1; i < M - 1; i++ ) { mean = mean + w->cell[i][N-1]; }
	for ( j = 0; j < N; j++ ) { mean = mean + w->cell[M-1][j]; }
	for ( j = 0; j < N; j++ ) { mean = mean + w->cell[0][j]; }
	mean = mean / ( double ) ( 2 * M + 2 * N - 4 );
	*/
	 
	//Initialize the interior solution to -----the mean value---- ZERO.
	for ( i = 1; i < M - 1; i++ )
	for ( j = 1; j < N - 1; j++ ) {
		w->cell[i][j] = HEAT_MIN;
	}
	
	int graphic_delay = 0;
	int16_t nframes = 0;
	
	//iterate until the  new solution W differs from the old solution U
	//by no more than EPSILON.
	iterations = 0;
	iterations_print = 1;
	printf ( "\n" );
	printf ( " Iteration  Change\n" );
	printf ( "\n" );
	ctime1 = cpu_time ( );

	while ( epsilon <= diff ) {

		//Save the old solution in U.
		for ( i = 0; i < M; i++ ) 
		for ( j = 0; j < N; j++ ) {
			u->cell[i][j] = w->cell[i][j];
		}
		
		//copy iteration into graphic_heat_plate array
		if (graphic_delay <= 0)
		if (nframes < MAX_STORED_FRAMES) {
			graphic_history[nframes] = malloc(sizeof(graphic_heat_plate));
			for (int py = 0; py < GRAPHIC_N; py++)
			for (int px = 0; px < GRAPHIC_N; px++) {
				float ptg_x = (float)px / (float)(GRAPHIC_N - 1);
				float ptg_y = (float)py / (float)(GRAPHIC_N - 1);
				int sample_x = iclamp(ptg_x * (N - 1), 0, N - 1);
				int sample_y = iclamp(ptg_y * (M - 1), 0, M - 1);
				double heat_sample = u->cell[sample_y][sample_x];
				graphic_history[nframes]->cell[px][py] =
					(uint8_t)iclamp(255 * ((heat_sample - HEAT_MIN) / (HEAT_MAX - HEAT_MIN)) , 0, 255);
			}
			
			graphic_delay = GRAPHIC_PRECISION + 1;
			nframes++;
		}
		
		if (graphic_delay > 0) graphic_delay--;
		
		//Determine the new estimate of the solution at the interior points.
		//The new solution W is the average of north, south, east and west neighbors.
		diff = 0.0;
		for ( i = 1; i < M - 1; i++ )
		for ( j = 1; j < N - 1; j++ ) {
			w->cell[i][j] = ( u->cell[i-1][j] + u->cell[i+1][j] + u->cell[i][j-1] + u->cell[i][j+1] ) / 4.0;

			if ( diff < fabs (w->cell[i][j] - u->cell[i][j])  ) {
				diff = fabs ( w->cell[i][j] - u->cell[i][j] );
			}
		}
		
		iterations++;
		if ( iterations == iterations_print )
		{
		  printf ( "  %8d  %f\n", iterations, diff );
		  iterations_print = 2 * iterations_print;
		}
	}
	
	ctime2 = cpu_time ( );
	ctime = ctime2 - ctime1;

	printf ( "\n" );
	printf ( "  %8d  %f\n", iterations, diff );
	printf ( "\n" );
	printf ( "  Error tolerance achieved.\n" );
	printf ( "  CPU time = %f\n", ctime );
	
	//Write the solution to the output file.
	FILE* f = fopen(output_file, "wb");
	if (f) {
		fwrite(&nframes, sizeof(int16_t), 1, f);
		for (int k = 0; k < nframes; k++) {
			for (int py = 0; py < GRAPHIC_N; py++)
			for (int px = 0; px < GRAPHIC_N; px++) {
				fwrite(&(graphic_history[k]->cell[px][py]), sizeof(uint8_t), 1, f);
			}
		}
		
		fclose(f);
		printf ( "\n" );
		printf ("  Solution written to the output file '%s'\n", output_file );	
	}



	printf ( "\n" );
	printf ( "HEATED_PLATE:\n" );
	printf ( "  Normal end of execution.\n" );
	printf ( "\n" );

	return 0;
}