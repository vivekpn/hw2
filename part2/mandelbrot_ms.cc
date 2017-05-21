/**
 *  \file mandelbrot_ms.cc
 *
 *  \brief Implement your parallel mandelbrot set in this file.
 */

#include <iostream>
#include <cstdlib>
#include "render.hh"
#include <mpi.h>
#include "timer.c"

using namespace std;
void printdata(double* data, int len);
int mandelbrot(double x, double y);

void printdata(double* data, int len)
{
    for (int i = 0; i < len; i++) {
        cout << data[i] << " ";
    }
    cout << endl;
}

int mandelbrot(double x, double y)
{
    int maxit = 511;
    double cx = x;
    double cy = y;
    double newx, newy;

    int it = 0;
    for (it = 0; it < maxit && (x * x + y * y) < 4; ++it) {
        newx = x * x - y * y + cx;
        newy = 2 * x * y + cy;
        x = newx;
        y = newy;
    }
    return it;
}

int main(int argc, char* argv[])
{
    int rank = 0;
    int np = 0;
    char hostname[MPI_MAX_PROCESSOR_NAME + 1];
    int namelen = 0;
    const int MSG_TAG_DATA = 1;
    const int MSG_TAG_COMPLETED = 0;

    /* Initialization */
    double minX = -2.1;
    double maxX = 0.7;
    double minY = -1.25;
    double maxY = 1.25;
    /* Master takes the input from the user */
    int height, width;

    
    if (rank == 0) {
        if (argc == 3) {
            height = atoi(argv[1]);
            width = atoi(argv[2]);
            assert(height > 0 && width > 0);
        }
        else {
            fprintf(stderr, "usage: %s <height> <width>\n", argv[0]);
            fprintf(stderr, "where <height> and <width> are the dimensions of the image.\n");
            return -1;
        }
    }

    struct stopwatch_t* timer;
    if(rank == 0){
		  stopwatch_init ();
			timer = stopwatch_create ();
			stopwatch_start (timer);
    }

    double it = (maxY - minY) / height;
    double jt = (maxX - minX) / width;

    MPI_Init(&argc, &argv); /* starts MPI */
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); /* Get process id */
    MPI_Comm_size(MPI_COMM_WORLD, &np); /* Get number of processes */
    MPI_Get_processor_name(hostname, &namelen); /* Get hostname of node */
    printf("Hello, world! [Host:%s -- Rank %d out of %d]\n", hostname, rank, np);

    /* code requires the data and processors to match */

    int blocksize = height / (np - 1);

    /* Master divides the work */

    MPI_Request requests[np - 1];
    MPI_Status statuses[np - 1];
    int row = 0;
    map<int, int> prows;
    int INITIAL_NUMBERS = 1;
    double* data = NULL;

    if (rank == 0) {
        for (int i = 1; i < np; i++) {
            prows[i] = INITIAL_NUMBERS;
        }
        /* Master sends work to all the slaves */
        cout << "Master is alloting work to slaves." << endl;
        for (int i = 1; i < np; i++) {
            cout << "Master is sending to " << i << endl;
            int sendbuf[2];
            sendbuf[0] = row;
            sendbuf[1] = prows[i];
            row = row + prows[i];
            // prows[i] = prows[i] * 2;
            MPI_Isend(&sendbuf, 2, MPI_INT, i, MSG_TAG_DATA, MPI_COMM_WORLD, &requests[i - 1]);
        }
        data = new double[height * width];
        double recvbuf[width + 1];
        while (row < height) {
            /* Master recieves slaves work*/
            cout << "Master is waiting to recieve data from a slave " << endl;
            MPI_Status status;            
            MPI_Recv(recvbuf, width + 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            int slave = status.MPI_SOURCE;
            cout << "Master recieved data from slave " << slave << endl;
            memcpy(data + lround(recvbuf[width] * width), recvbuf, width * prows[slave] * sizeof(double));

            /* Master sends new work to all the slaves */
            int sendbuf[2];
            sendbuf[0] = row;
            sendbuf[1] = prows[slave];
            MPI_Send(&sendbuf, 2, MPI_INT, slave, MSG_TAG_DATA, MPI_COMM_WORLD);
            row = row + prows[slave];
        }

        /* Collect the final results. */
				cout << "Master is starting to collect the final results." << endl;
        for (int i = 1; i < np; i++) {
        		MPI_Status status;
            MPI_Recv(recvbuf, width + 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            cout << "Collected the final results from " << status.MPI_SOURCE << endl;            
            // store the received result:
            memcpy(data + lround(recvbuf[width] * width), recvbuf, width * sizeof(double));
            cout << "Master is sending done to the slave " << status.MPI_SOURCE << endl;
            int sendbuf[2];
            sendbuf[0] = 0;
            sendbuf[1] = 0;
            MPI_Send(&sendbuf, 2, MPI_INT, status.MPI_SOURCE, MSG_TAG_COMPLETED, MPI_COMM_WORLD);
        }
    }

    int lrow;
    int lsize;
    if (rank != 0) {
        while (true) {
            /* Each slave recieves its work */
            cout << "Slave " << rank << " is recieving." << endl;
            MPI_Status stat;
            int rec_buf[2];
            MPI_Recv(&rec_buf, 2, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
            lrow = rec_buf[0];
            lsize = rec_buf[1];
            cout << "Slave "<< rank <<" recieved the rows starting from " << lrow << " of size " << lsize << endl;
            if (lsize == 0) {
                break;
            }
            /* Slaves perform the work */
            int BUFFER_LENGTH = width * lsize + 1;
            cout << "Slave " << rank << " is starting the work" << endl;
            double* ldata = new double[BUFFER_LENGTH];
            double x, y;
            y = minY + lrow * it;
            for (int i = 0; i < lsize; i++) {
                x = minX;
                for (int j = 0; j < width; ++j) {
                    ldata[i * width + j] = mandelbrot(x, y) / 512.0;
                    x += jt;
                }
                y += it;
            }
            /* Each slave sends the results to the master */
            ldata[BUFFER_LENGTH - 1] = lrow;
            cout << "Slave " << rank << " is sending data to master" << endl;
            MPI_Send(ldata, BUFFER_LENGTH, MPI_DOUBLE, 0, MSG_TAG_DATA, MPI_COMM_WORLD);
        }
    }

    cout << "Initial send is complete for process " << rank << endl;
    cout << "Height " << height << " width " << width << " rank " << rank << endl;

    /* Once the work is completed, master renders the image and outputs the result */
    if (rank == 0) {
        gil::rgb8_image_t img(height, width);
        auto img_view = gil::view(img);
        for (int i = 0; i < height; ++i) {
            for (int j = 0; j < width; ++j) {
                img_view(j, i) = render(data[j * width + i]);
            }
        }
        gil::png_write_view("mandelbrot.png", const_view(img));
    }
    cout << "Rank " << rank << " before finalize." << endl;
    MPI_Finalize();
    cout << "Rank " << rank << " after finalize." << endl;
    
   if(rank == 0){
        long double t_qs = stopwatch_stop (timer);
		    printf ("Mandelbrot total execution time: %Lg seconds", t_qs);
	      stopwatch_destroy (timer);
    }

}
