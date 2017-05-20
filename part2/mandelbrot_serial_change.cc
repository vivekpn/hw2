/**
 *  \file mandelbrot_serial.cc
 *  \brief Lab 2: Mandelbrot set serial code
 */


#include <iostream>
#include <cstdlib>

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
  double minX = -2.1;
  double maxX = 0.7;
  double minY = -1.25;
  double maxY = 1.25;
  
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

  double it = (maxY - minY)/height;
  double jt = (maxX - minX)/width;
  double x, y;
  
  gil::rgb8_image_t img(height, width);
  auto img_view = gil::view(img);

  /*Starting the parallel section  */

  int rank = 0;
  int np = 0;
  char hostname[MPI_MAX_PROCESSOR_NAME+1];
  int namelen = 0;

  FILE *fp = NULL; /* output file, only valid on rank 0 */

  int* msgbuf = NULL;
  int len = 0;
  const int MSG_TAG = 1;
  
  
  MPI_Init (&argc, &argv);	/* starts MPI */
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);	/* Get process id */
  MPI_Comm_size (MPI_COMM_WORLD, &np);	/* Get number of processes */
  MPI_Get_processor_name (hostname, &namelen); /* Get hostname of node */
  printf ("Hello, world! [Host:%s -- Rank %d out of %d]\n", hostname, rank, np);
  
  if(rank ==0)
  {
   int* msgbuf = malloc (height*width * sizeof (int));
  }
 /* Handle for uneven data or procesors */
  int lheight = height/np; 
  
   int  MAX_BUFLEN = lheight *width;
 // /*Slaves should recive all the rows from master */ 
   int* buff = malloc (MAX_BUFLEN * sizeof (int));

   y = minY+ rank*it;
   for (int i = 0; i < lheight; ++i) {
   x = minX;
    for (int j = 0; j < width; ++j) {
      buff[i][j]= mandelbrot(x, y)/512.0);
      x += jt;
    }
    y += it;
   }
   
   if (rank == 0) {
    MPI_Status stat;
    memcpy(msgbuff,buff,sizeof(int)*MAX_BUFLEN);
 
    for(int i=1;i<np;i++)
    MPI_Recv (msgbuf[i*lheight], MAX_BUFLEN, MPI_INT, i, MSG_TAG, MPI_COMM_WORLD);
  } else {
    MPI_Status stat;
    MPI_Send (buff, MAX_BUFLEN, MPI_INT, 0, MSG_TAG, MPI_COMM_WORLD, &stat);
  }

  
  MPI_Finalize ();
  for (int i = 0; i < height; ++i) {
    for (int j = 0; j < width; ++j) {
      /* Image j and i conflict*/
      img_view(j, i) = render(msgbuf[j][i]);
    }
  } 

  /*Do rendering and image work */
  gil::png_write_view("mandelbrot.png", const_view(img));
}

/* eof */
