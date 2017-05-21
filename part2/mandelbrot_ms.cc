/**
 *  \file mandelbrot_ms.cc
 *
 *  \brief Implement your parallel mandelbrot set in this file.
 */

#include <iostream>
#include <cstdlib>
#include "render.hh"
#include <mpi.h>

using namespace std;
void printdata(double* data, int len);
int mandelbrot(double x, double y);

void printdata(double* data, int len){
	for(int i=0; i< len; i++){
		cout<< data[i] << " ";
	}
	cout<< endl;
}

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
  const int MSG_TAG2 = 2;
  FILE *fp = NULL; /* output file, only valid on rank 0 */

  int len = 0;  
  
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
  
//  blocksize = 8;
//  height = (np-1)*blocksize;
//  width = (np-1)*blocksize;
  
  double it = (maxY - minY)/height;
  double jt = (maxX - minX)/width;

  MPI_Init (&argc, &argv);	/* starts MPI */
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);	/* Get process id */
  MPI_Comm_size (MPI_COMM_WORLD, &np);	/* Get number of processes */
  MPI_Get_processor_name (hostname, &namelen); /* Get hostname of node */
  printf ("Hello, world! [Host:%s -- Rank %d out of %d]\n", hostname, rank, np);

  /* code requires the data and processors to match */
  assert (height % (np-1) == 0);

  int blocksize = height/(np-1);

  /* Master divides the work */
  

  /* Master sends work to all the slaves */   
  MPI_Request requests[np-1];
  MPI_Status statuses[np-1];
  if (rank == 0) {
    cout << "Master is trying to send." << endl;
  	for(int i = 1; i< np; i++){
  	    cout << "Master is sending to " << i << endl;
	    int sendbuf = i-1;
	    MPI_Isend (&sendbuf, 1, MPI_INT, i, MSG_TAG, MPI_COMM_WORLD, &requests[i-1]);
  	}
//  MPI_Waitall(np-1, requests, statuses);
  } 
      
  /* Each slave recieves its work */
  int lrow = 0;
  if(rank != 0){
        cout << "Slave "<< rank << " is trying to recieve." << endl;
		MPI_Status stat;
		MPI_Recv (&lrow, 1, MPI_INT, 0, MSG_TAG, MPI_COMM_WORLD, &stat);
       cout << "Slave recieved the value " << lrow << endl;
  }  
 
  cout << "Initial send is complete for process " << rank << endl; 
  /* Slaves perform the work */  
  int BUFFER_LENGTH = blocksize * width;
  if(rank != 0 ){
          cout << "Slave " << rank << " is starting the work" << endl;	  
	  double* ldata = new double[BUFFER_LENGTH];
	  double x, y;
	  y = minY + (rank-1) * blocksize * it;
	  for (int i = 0; i < blocksize; i++) {
		x = minX;
		for (int j = 0; j < width; ++j) {
		  ldata[i*width+j] = mandelbrot(x, y)/512.0;
		  x += jt;
		}
		y += it;
	  }
	  /* Each slave sends the results to master */
	  cout << "Slave "<< rank << " is sending data to master" << endl;
//	  printdata(ldata, BUFFER_LENGTH);
      MPI_Send (ldata, BUFFER_LENGTH, MPI_DOUBLE, 0, MSG_TAG2, MPI_COMM_WORLD);
  }

  cout << "Height " << height << " width " << width << " rank " << rank << endl;
  double* data = NULL;
  /* Master recieves slaves work*/   

  if(rank == 0) {
    cout << 126 << endl; 
    data = new double[height*width];
    cout << 128 << endl; 

  	for(int i = 1; i< np; i++){
        cout << "Master is recieving data from slave " << i << endl;
		MPI_Irecv (data + (BUFFER_LENGTH * (i-1)), BUFFER_LENGTH, MPI_DOUBLE, i, MSG_TAG2, MPI_COMM_WORLD, &requests[i-1]);  
		cout << "Master recieved data from slave " << i << endl;
  	}
  	MPI_Waitall(np-1, requests, statuses);
  }
  

//  cout << "Rank " << rank << " before barrier." << endl;
  /* Master gathers the results and sends more work for the slave */
//  MPI_Barrier (MPI_COMM_WORLD); /* Synchronize the nodes */
//  cout << "Rank " << rank << " after barrier." << endl;

  /* Once the work is completed, master renders the image and outputs the result */  
  if(rank == 0) {
//      cout << 149 << endl;
	  gil::rgb8_image_t img(height, width);
//	  cout << 151 << endl;

	  auto img_view = gil::view(img);
//	  cout << 154 << endl;

	  for (int i = 0; i < height; ++i) {
//	  	cout << 157 << endl;
		for (int j = 0; j < width; ++j) {
//	      cout << 159 << endl;
		  img_view(j, i) = render(data[j*width+i]);
//		  cout << 161 << endl;
		}
	  }
//	  cout << 164 << endl;
	  gil::png_write_view("mandelbrot.png", const_view(img));  
  }
  cout << "Rank " << rank << " before finalize." << endl;  
  MPI_Finalize (); 
  cout << "Rank " << rank << " after finalize." << endl;  
}
