/**
 *  \file mandelbrot_ms.cc
 *
 *  \brief Implement your parallel mandelbrot set in this file.
 */

#include <iostream>
#include <cstdlib>
#include "render.hh"
#include <mpi.h>

#define WIDTH 1000
#define HEIGHT 1000

using namespace std;

int mandelbrot(double x, double y) {
  int maxit = 511;
  double cx = x;
  double cy = y;
  double newx, newy;

  int it = 0;
  for (it = 0; it < maxit && (x*x + y*y) < 4; ++it) {
    newx = x*x - y*y + cx;
    newy = 2*x*y + cy;
    x = newx;
    y = newy;
  }
  return it;
}


int main(int argc, char* argv[]) {
  int rank = 0;
  int np = 0;
  char hostname[MPI_MAX_PROCESSOR_NAME+1];
  int namelen = 0;
  const int MSG_TAG = 1;

  FILE *fp = NULL; /* output file, only valid on rank 0 */

  int len = 0;

  MPI_Init (&argc, &argv);	/* starts MPI */
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);	/* Get process id */
  MPI_Comm_size (MPI_COMM_WORLD, &np);	/* Get number of processes */
  MPI_Get_processor_name (hostname, &namelen); /* Get hostname of node */
  printf ("Hello, world! [Host:%s -- Rank %d out of %d]\n", hostname, rank, np);
  
  
  /* Initialization */
  double minX = -2.1;
  double maxX = 0.7;
  double minY = -1.25;
  double maxY = 1.25;
  /* Master takes the input from the user */
  int height, width;
  if(rank == 0){
	  if (argc == 3) {
		height = atoi (argv[1]);
		width = atoi (argv[2]);
		assert (height > 0 && width > 0);
	  } else {
		fprintf (stderr, "usage: %s <height> <width>\n", argv[0]);
		fprintf (stderr, "where <height> and <width> are the dimensions of the image.\n");
		return -1;
	  }
  }
    
  height = 128;
  width = 128;
  
  double it = (maxY - minY)/height;
  double jt = (maxX - minX)/width;

  /* Master divides the work */
  

  /* Master sends work to all the slaves */   
  if (rank == 0) {
  	for(int i = 1; i< np; i++){
	    int value = i-1;
	    MPI_Send (&value, 1, MPI_INT, i, MSG_TAG, MPI_COMM_WORLD);	  	
  	}
  } 
  
  /* Each slave recieves its work */
  int* msgbuf = NULL;
  if(rank != 0){
		msgbuf = new int[1];
		MPI_Status stat;
		MPI_Recv (msgbuf, 1, MPI_INT, 0, MSG_TAG, MPI_COMM_WORLD, &stat);  
  }
  int lrow = msgbuf[0];
    
  /* Slaves perform the work */  
  int BUFFER_LENGTH = 1 * width;
  if(rank != 0 ){	  
	  int* ldata = (int *) malloc(BUFFER_LENGTH * sizeof(int));
	  double x, y;
	  y = minY + (rank-1) * it;
	  for (int i = lrow; i < 1; ++i) {
		x = minX;
		for (int j = 0; j < width; ++j) {
		  ldata[i*width+j] = mandelbrot(x, y)/512.0;
		  x += jt;
		}
		y += it;
	  }
	  /* Each slave sends the results to master */
	  MPI_Send (ldata, width, MPI_INT, 0, MSG_TAG, MPI_COMM_WORLD);	  	
  }

  int** data = NULL;
  /* Master recieves slaves work*/   
  if(rank == 0) {
    data = new int*[width];
    // data = (int *) malloc( BUFFER_LENGTH * sizeof(int));
  	for(int i = 1; i< np; i++){
  		MPI_Status stat;
		MPI_Recv (data[i], width, MPI_INT, i, MSG_TAG, MPI_COMM_WORLD, &stat);  
  	}
  }
  
  /* Master gathers the results and sends more work for the slave */
  MPI_Barrier (MPI_COMM_WORLD); /* Synchronize the nodes */
 
  /* Once the work is completed, master renders the image and outputs the result */  
  if(rank == 0) {
	  gil::rgb8_image_t img(height, width);
	  auto img_view = gil::view(img);
	  for (int i = 0; i < height; ++i) {
		for (int j = 0; j < width; ++j) {
		  img_view(j, i) = render(data[j][i]);
		}
	  }
	  gil::png_write_view("mandelbrot.png", const_view(img));  
  }
  
  MPI_Finalize (); 
}
