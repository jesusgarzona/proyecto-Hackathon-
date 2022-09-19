#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <stdint.h>
#include "mpi.h"

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
#define GRAPHIC_N 100

#define MAX_STORED_FRAMES 100

#define HEAT_MAX (100.0)
#define HEAT_MIN (0.0)

//aka each 50 simulation steps, a graphic frame is stored
#define GRAPHIC_PRECISION 200

#define OUTPUT_FILE "aber.data"

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

static inline int imax(int a, int b) {
	if (a > b) return a;
	return b;
}

int main ( int argc, char *argv[] ) {	

  int size, rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); 

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
	int success;
	heat_plate* u = malloc(sizeof(heat_plate));
	heat_plate* w = malloc(sizeof(heat_plate));
	graphic_heat_plate** graphic_history = malloc(sizeof(graphic_heat_plate*) * MAX_STORED_FRAMES );
 
  if (argc < 2) {
    if (rank == 0) {
      printf("Usage: heated_plate <epsilon>\n");
    }
    return 0;
  }
  
  epsilon = atof(argv[1]);
	diff = epsilon;
	
	//Set the boundary values, which don't change. 
	for ( i = 1; i < M - 1; i++ ) { w->cell[i][0]   = HEAT_MAX; }
	for ( i = 1; i < M - 1; i++ ) { w->cell[i][N-1] = HEAT_MAX; }
	for ( j = 0; j < N; j++ ) 	  { w->cell[M-1][j] = HEAT_MAX; }
	for ( j = 0; j < N; j++ ) 	  { w->cell[0][j]   = HEAT_MIN; }
	 
	//Initialize the interior solution to -----the mean value---- ZERO.
	for ( i = 1; i < M - 1; i++ )
	for ( j = 1; j < N - 1; j++ ) {
		w->cell[i][j] = HEAT_MIN;
	}
	
	int graphic_delay = 0;
	int16_t nframes = 0;
 
  //each rank will do the calculations for a particular horizontal slice of the surface.
  //Each slice has the whole X range, but only a partial section of the Y range.
  int ystart = imax(1, (rank * M) / size);
  int yend   = imin(M - 1, ((rank + 1) * M) / size);
  
  //same thing but used in graphic stuff 
  int ystart_graph =  (rank * GRAPHIC_N) / size;
  int yend_graph   = ((rank + 1) * GRAPHIC_N) / size;  
	
	//iterate until the  new solution W differs from the old solution U
	//by no more than EPSILON.
	iterations = 0;
	iterations_print = 1;
	ctime1 = cpu_time ( );

	while ( epsilon <= diff ) {

		//Save the old solution in U.
		for ( i = 0; i < M; i++ ) 
		for ( j = 0; j < N; j++ ) {
			u->cell[i][j] = w->cell[i][j];
		}
		
		//copy iteration into graphic_heat_plate array.
    // any rank will only do its particular slice
		if (graphic_delay <= 0)
		if (nframes < MAX_STORED_FRAMES) {
			graphic_history[nframes] = malloc(sizeof(graphic_heat_plate));
			for (int py = ystart_graph; py < yend_graph; py++)
			for (int px = 0; px < GRAPHIC_N; px++) {
				float ptg_x = (float)px / (float)(GRAPHIC_N - 1);
				float ptg_y = (float)py / (float)(GRAPHIC_N - 1);
				int sample_x = iclamp(ptg_x * (N - 1), 0, N - 1);
				int sample_y = iclamp(ptg_y * (M - 1), 0, M - 1);
				double heat_sample = u->cell[sample_y][sample_x];
				graphic_history[nframes]->cell[px][py] =
					(uint8_t)iclamp(255 * ((heat_sample - HEAT_MIN) / (HEAT_MAX - HEAT_MIN)) , 0, 255);
			}
      
      //once the graphic frame slice was calculated, the data broadcasts begin.
      //ranks will broadcast in order (0, 1, 2, etc) 
      
      for (int broadcast_rank = 0; broadcast_rank < size; broadcast_rank++) {
        if (rank == broadcast_rank) {
          for (int target_rank = 0; target_rank < size; target_rank++) {
            if (target_rank != rank) {
              MPI_Send(&ystart_graph, 1, MPI_INT, target_rank, 0, MPI_COMM_WORLD);
              MPI_Send(&yend_graph, 1, MPI_INT, target_rank, 0, MPI_COMM_WORLD);
              //data will be sent as for (y) { for (x) }
              for (int py = ystart_graph; py < yend_graph; py++)
              for (int px = 0; px < GRAPHIC_N; px++) {
                MPI_Send(&( graphic_history[nframes]->cell[px][py] ), 1, MPI_CHAR, target_rank, 0, MPI_COMM_WORLD);
              }
            }
          }
        }
        else {
          int other_ystart_graph;
          int other_yend_graph;
          MPI_Recv(&other_ystart_graph, 1, MPI_INT, broadcast_rank, 0, MPI_COMM_WORLD, NULL);
          MPI_Recv(&other_yend_graph, 1, MPI_INT, broadcast_rank, 0, MPI_COMM_WORLD, NULL);
          for (int py = other_ystart_graph; py < other_yend_graph; py++)
          for (int px = 0; px < GRAPHIC_N; px++) {
            MPI_Recv(&( graphic_history[nframes]->cell[px][py] ), 1, MPI_CHAR, broadcast_rank, 0, MPI_COMM_WORLD, NULL);
          }
        }
      }
			
			graphic_delay = GRAPHIC_PRECISION + 1;
			nframes++;
		}
		
		if (graphic_delay > 0) graphic_delay--;
		
		//Determine the new estimate of the solution at the interior points.
		//The new solution W is the average of north, south, east and west neighbors.
   
   
		diff = 0.0;
		for ( i = ystart; i < yend; i++ )
		for ( j = 1; j < N - 1; j++ ) {
			w->cell[i][j] = ( u->cell[i-1][j] + u->cell[i+1][j] + u->cell[i][j-1] + u->cell[i][j+1] ) / 4.0;

			if ( diff < fabs (w->cell[i][j] - u->cell[i][j])  ) {
				diff = fabs ( w->cell[i][j] - u->cell[i][j] );
			}
		}
   
    //Same broadcasts as before. For now just send grid data. Rank 0 will compare <diff> values after this
    for (int broadcast_rank = 0; broadcast_rank < size; broadcast_rank++) {
      if (rank == broadcast_rank) {
        for (int target_rank = 0; target_rank < size; target_rank++) {
          if (target_rank != rank) {
            MPI_Send(&ystart, 1, MPI_INT, target_rank, 0, MPI_COMM_WORLD);
            MPI_Send(&yend, 1, MPI_INT, target_rank, 0, MPI_COMM_WORLD);
            //data will be sent as for (y) { for (x) }
            for (int py = ystart; py < yend; py++)
            for (int px = 1; px < N - 1; px++) {
              MPI_Send(&( w->cell[py][px] ), 1, MPI_DOUBLE, target_rank, 0, MPI_COMM_WORLD);
            }
          }
        }
      }
      else {
        int other_ystart;
        int other_yend;
        MPI_Recv(&other_ystart, 1, MPI_INT, broadcast_rank, 0, MPI_COMM_WORLD, NULL);
        MPI_Recv(&other_yend, 1, MPI_INT, broadcast_rank, 0, MPI_COMM_WORLD, NULL);
        for (int py = other_ystart; py < other_yend; py++)
        for (int px = 1; px < N - 1; px++) {
          MPI_Recv(&( w->cell[py][px] ), 1, MPI_DOUBLE, broadcast_rank, 0, MPI_COMM_WORLD, NULL);
        }
      }
    }
    
    //now everyone sends their <diff> to rank 0 and it compares
    if (rank > 0) {
      MPI_Send(&( diff ), 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    else {
      for (int other_rank = 1; other_rank < size; other_rank++) {
        double other_diff;
        MPI_Recv(&( other_diff ), 1, MPI_DOUBLE, other_rank, 0, MPI_COMM_WORLD, NULL);
        if (other_diff > diff) {
          diff = other_diff;
        }
      }
    }
    
    //now rank 0 sends actual diff value to everyone
    if (rank > 0) {
      MPI_Recv(&( diff ), 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, NULL);
    }
    else {
      for (int other_rank = 1; other_rank < size; other_rank++) {
        MPI_Send(&( diff ), 1, MPI_DOUBLE, other_rank, 0, MPI_COMM_WORLD);
      }
    }
		
		iterations++;
    if (rank == 0)
    if ( iterations == iterations_print )
    {
      printf ( "  %8d  %f\n", iterations, diff );
      iterations_print = 2 * iterations_print;
    }
	}

	if (rank == 0) {
    ctime2 = cpu_time ( );
    ctime = ctime2 - ctime1;
    
    printf ( "\n" );
    printf ( "  %8d  %f\n", iterations, diff );
    printf ( "\n" );
    printf ( "  Error tolerance achieved.\n" );
    printf ( "  CPU time = %f\n", ctime );
	
  	//Write the solution to the output file.
  	FILE* f = fopen(OUTPUT_FILE, "wb");
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
  		printf ("  Solution written to the output file '%s'\n", OUTPUT_FILE );	
  	}
  
  	printf ( "\n" );
  	printf ( "HEATED_PLATE:\n" );
  	printf ( "  Normal end of execution.\n" );
  	printf ( "\n" );
  }
  
  MPI_Finalize();
  return 0;
}