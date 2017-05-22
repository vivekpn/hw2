/**
 *  \file mandelbrot_serial.cc
 *  \brief Lab 2: Mandelbrot set serial code
 */


#include <iostream>
#include <cstdlib>
#include <mpi.h>

#include "render.hh"

using namespace std;

#define WIDTH 1000
#define HEIGHT 1000

int
mandelbrot(double x, double y) {
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



int
main(int argc, char* argv[]) {
  
  int height, width;
  if (argc == 3) {
    height = atoi (argv[1]);
    width = atoi (argv[2]);
    assert (height > 0 && width > 0);
  } else {
    fprintf (stderr, "usage: %s <height> <width>\n", argv[0]);
    fprintf (stderr, "where <height> and <width> are the dimensions of the image.\n");
    return -1;
  }

   double t_start, t_elapsed;
  gil::rgb8_image_t img(height, width);
  auto img_view = gil::view(img);

  /*Starting the parallel section  */
  MPI_Init (&argc, &argv);	/* starts MPI */
  int rank = 0;
  int np = 0;
  char hostname[MPI_MAX_PROCESSOR_NAME+1];
  int namelen = 0;

  FILE *fp = NULL; /* output file, only valid on rank 0 */

  int len = 0;

  MPI_Comm_rank (MPI_COMM_WORLD, &rank);	/* Get process id */
  MPI_Comm_size (MPI_COMM_WORLD, &np);	/* Get number of processes */
  MPI_Get_processor_name (hostname, &namelen); /* Get hostname of node */
  float* msgbuf = NULL;
  float* buff =NULL;
  
  const int MSG_TAG = 1;
  //printf ("Hello, world! [Host:%s -- Rank %d out of %d]\n", hostname, rank, np);
  double minX = -2.1;
  double maxX = 0.7;
  double minY = -1.25;
  double maxY = 1.25;
  MPI_Barrier (MPI_COMM_WORLD); /* Synchronize the nodes */
    if(rank ==0)
    t_start = MPI_Wtime (); /* Start timer */
  

  //printf ("Height = %d  Width= %d  rank= %d\n",height,width,rank);
  double it = (maxY - minY)/height;
  double jt = (maxX - minX)/width;
  double x, y;
    
 // Handle for uneven data or procesors 
  int lheight = height/np; 

  if(rank ==0)
  {
   //printf ("msgbuff allocated by rank %d \n", rank);
   msgbuf = (float *)malloc (height*width * sizeof (float));
   //printf ("msgbuff allocated by rank %d \n", rank);
  }

  //printf ("lheight = %d rank= %d \n",lheight, rank);
  int  MAX_BUFLEN = lheight *width;
  buff = (float *)malloc (MAX_BUFLEN * sizeof (int));
   
   
   y = minY+ (lheight)*rank*it;
   for (int i = 0; i < lheight; ++i) {
   x = minX;
    for (int j = 0; j < width; ++j) {
    buff[i*width+j]= (mandelbrot(x, y)/512.0);
      x += jt;
    }
    y += it;
   }
  
   //printf ("mandralbot caluculations over  rank =%d\n",rank);

/*   if (rank == 0) {
    MPI_Status stat;
    memcpy(msgbuf,buff,sizeof(float)*MAX_BUFLEN);
   
    for(int i=1;i<np;i++)
    MPI_Recv (&msgbuf[i*lheight*width], MAX_BUFLEN, MPI_FLOAT, i, MSG_TAG, MPI_COMM_WORLD,&stat);
    
  } else {
    MPI_Status stat;
    MPI_Send (buff, MAX_BUFLEN, MPI_FLOAT, 0, MSG_TAG, MPI_COMM_WORLD);
    printf ("sent  data  rank =%d\n",rank);
  }
*/
  //Gather all the data
  MPI_Gather(buff,MAX_BUFLEN,MPI_FLOAT,msgbuf,MAX_BUFLEN,MPI_FLOAT,0,MPI_COMM_WORLD);    
  //printf ("Gather completed by rank =%d\n",rank);

  if(rank ==0)
  { 
  //printf ("exited parallel section rank =%d\n",rank);
  
 // for (int i = 0; i < height; ++i) {
 for (int i = 0; i < lheight*np; ++i) {   
    for (int j = 0; j < width; ++j) {
      //Image j and i conflict
      img_view(j, i) = render(msgbuf[i*width+j]);
      //printf("%f ",msgbuf[i*width+j]);
    } 
      //printf("\n");
  }
 
 y = minY+lheight*np*it;   
 for (int i = lheight*np; i < height; ++i) {
    x = minX;   
    for (int j = 0; j < width; ++j) {
      //Image j and i conflict
      img_view(j, i) = render(mandelbrot(x, y)/512.0);
      x += jt;
    } 
      y += it;
  }
 
  //printf ("Image computations complete! rank =%d\n",rank);

  //Do rendering and image work 
  gil::png_write_view("mandelbrot.png", const_view(img));
  //printf ("Writing the image file  rank =%d\n",rank);
  }

  MPI_Barrier (MPI_COMM_WORLD); /* Synchronize the nodes */
  if(rank ==0)
  {
  t_elapsed = MPI_Wtime () - t_start; /* Stop timer */
   printf ("Mandelbrot total execution time: %1.2f seconds", t_elapsed);
  }
	
  MPI_Finalize ();
}

/* eof */


