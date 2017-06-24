#include "mpi.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define NUM_POINTS 524288

unsigned int X_axis[NUM_POINTS];
unsigned int Y_axis[NUM_POINTS];

struct retval compute_dimension(int num_quadrants, int numprocs);
int power(int x, int y);
void seperate_data_x_axis(int x_num, int y_num);
void seperate_data_y_axis(int x_num, int y_num);
void partition(int offset, int x_num, int y_num, int num_quadrants);
void quadrants_coord(int x_num, int y_num, int offset, int num_quadrants);
unsigned int find_kth(unsigned int *x, int n, int k, unsigned int *y);
void swap(unsigned int array[], int i, int j);

int i = 0;
int j = 0;
int k = 0;
int z = 0;

double global_cost = 0;

int numprocs;  /* Number of processors to use */
int myid;

int num_quadrants;

struct retval
{
    int x;
    int y;
};

struct retval compute_dimenssion(int num_quadrants, int numprocs){
    /**************************
    # This function compute the dimension r.x, r.y
    # Input:
    #   num_quadrants: the number of quadrants we need to implement
    #   numprocs: the number of processors we have
    # Output:
    #   r: type retval, r.x is the x dimension, r.y is the y dimension.
    **************************/
    struct retval r;
    int i = 1;
    int num = 0;
    if(num_quadrants < numprocs){
        while(i < num_quadrants){
            i = i * 2;
            num ++;
        }
    }else{
        while(i < numprocs){
            i = i * 2;
            num ++;
        }
    }
    r.x = num /2;
    r.y = num - num/2;
    return r;
}

int power(int x, int y){
    int num = 1;
    int i;
    for(i = 0; i< y; i++){
        num = num * x;
    }
    return num;
}

void seperate_data_x_axis(int x_num, int y_num){
    /**************************
    # This function separate data in x-coordinates
    # Input:
    #   x_num: the number of processor firstly partition the data.
    #   y_num: after first partition, each processor then partition data into y_num parts, each part 
    #          is sent to a processor.
    # Output:
    **************************/
    int points = NUM_POINTS / x_num;
    if (myid == 0){
        int i;
        int x_arr[x_num];
        for(i = 0; i<x_num - 1; i++){
            int x_pivot = find_kth(X_axis + i * points, points * (x_num-i), points, Y_axis + i * points);
            int k = i * points;
            int j = NUM_POINTS -1;

            while (k < j && k < i * points + points/2) {
                if (X_axis[k] > x_pivot) {
                    while (X_axis[j] > x_pivot) {
                        j--;
                    }
                    if (k < j) {
                        swap(X_axis, k, j);
                        swap(Y_axis, k, j);
                    } 
                }
                k++;
            }
            x_arr[i] =  x_pivot;
            if (i > 0){
                MPI_Send(X_axis + i * points, points, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(Y_axis + i * points, points, MPI_INT, i, 0, MPI_COMM_WORLD);
            }  
        }
        if(x_num > 1){
            MPI_Send(X_axis + (x_num-1) * points, points, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(Y_axis + (x_num-1) * points, points, MPI_INT, i, 0, MPI_COMM_WORLD);       
        }
    }else{
         MPI_Recv(X_axis + myid * points, points, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
         MPI_Recv(Y_axis + myid * points, points, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

void seperate_data_y_axis(int x_num, int y_num){
    /**************************
    # This function separate data in y-coordinates
    # Input:
    #   x_num: the number of processor firstly partition the data.
    #   y_num: after first partition, each processor then partition data into y_num parts, each part 
    #          is sent to a processor.
    # Output:
    **************************/
    int points = NUM_POINTS / (x_num * y_num);
    int interval = NUM_POINTS / x_num ;

    if (myid < x_num){
        int i;
        // int y_arr[x_num * y_num];
        for(i = 0; i<y_num - 1; i++){ 
            int y_pivot = find_kth(Y_axis + myid * interval + i * points, points * (y_num -i), points, X_axis + myid * interval + i * points);
            int k = i * points;
            int j = interval -1; 

            while (k < j && k < i * points + points/2) {
                if (Y_axis[myid * interval + k] > y_pivot) {
                    while (Y_axis[myid * interval + j] > y_pivot) {
                        j--;
                    }
                    if (k < j) {
                        swap(X_axis, myid * interval + k, myid * interval + j);
                        swap(Y_axis, myid * interval + k, myid * interval + j);
                    } 
                }
                k++;
            }
            // y_arr[myid * y_num + i] = y_pivot; 
            if (i > 0){
                MPI_Send(X_axis + myid * interval + i * points, points, MPI_INT, myid + x_num * i, 0, MPI_COMM_WORLD);
                MPI_Send(Y_axis + myid * interval + i * points, points, MPI_INT, myid + x_num * i, 0, MPI_COMM_WORLD);
            }         
        }
        if (y_num > 1){
            MPI_Send(X_axis + myid * interval + (y_num-1) * points, points, MPI_INT, myid + x_num * (y_num-1), 0, MPI_COMM_WORLD);
            MPI_Send(Y_axis + myid * interval + (y_num-1) * points, points, MPI_INT, myid + x_num * (y_num-1), 0, MPI_COMM_WORLD);
        }
    }else{
        MPI_Recv(X_axis + (myid % x_num) * interval + (myid / x_num) * points, points, MPI_INT, myid % x_num , 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MPI_Recv(Y_axis + (myid % x_num) * interval + (myid / x_num) * points, points, MPI_INT, myid % x_num , 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }

    // if (myid > 0 ){
    //     MPI_Send(y_arr + myid * y_num , y_num, MPI_INT, 0, 0, MPI_COMM_WORLD);
    // }else{
    //     int i;
    //     for(i = 1; i< x_num; i++){
    //         MPI_Recv(y_arr + myid * y_num , y_num, MPI_INT, 0, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    //     }
    // }
}


void partition(int offset, int x_num, int y_num, int num_quadrants){
    /**************************
    # This function does partition in a certain processor
    # Input:
    #   offset: the address of data this processor received minus base address(X_axis, Y_axis).
    #   x_num: the number of processor firstly partition the data.
    #   y_num: after first partition, each processor then partition data into y_num parts, each part 
    #          is sent to a processor.
    # Output:
    #   pivot_arr: record bisection coordinates in this processor.
    **************************/

    int x_cut = 0;
    int quadrants = 1;
    // the points in each processor
    int proc_points = NUM_POINTS / (x_num * y_num);
    // the quadrants in each processor
    int proc_quad = num_quadrants / (x_num * y_num);
    // printf("%d\n",num_quadrants );
    //record the bisection coordinates
    int pivot_arr[proc_quad];
    int block_points;

    if (myid == 0){
        printf("\nThe FINAL Patition Begins!\nmyid : 0\nproc_quad : %d \nnum_quadrants : %d \n",proc_quad, num_quadrants );
    }

    // the quadrants in each processor is less than proc_quad
    while (proc_quad > quadrants) { 
        block_points = NUM_POINTS / (quadrants * x_num * y_num);

        //bisection on x-coordinate
        if (!x_cut) {
            for (i = 0; i < quadrants; i++) {
                //find the median
                int x_pivot = find_kth(X_axis + offset + i * block_points, block_points, block_points/2 - 1, Y_axis + offset + i * block_points);
                // if (myid == 0){
                //     printf("my id  =%d, value is %d \n", myid, x_pivot);
                // }
                k = i * block_points;
                j = i * block_points + block_points -1;
                //values less than median are put on the left half of the array
                while (k < j && k < i * block_points + block_points/2) {
                    if (X_axis[offset + k] > x_pivot) {
                        while (X_axis[offset + j] > x_pivot) {
                            j--;
                        }
                        if (k < j) {
                            swap(X_axis, k, j);
                            swap(Y_axis, k, j);
                        } 
                    }
                    k++;
                }
            }
            // next partition on y-coordinate
            x_cut = 1;  
        } else {
          //bisection on Y coordinate
            for (i = 0; i < quadrants; i++) {
                int y_pivot = find_kth(Y_axis + offset + i * block_points, block_points, block_points/2 - 1, X_axis + offset + i * block_points);

                k = offset + i * block_points;
                j = offset + i * block_points + block_points - 1;
                
                while (k < j && k < offset + i * block_points + block_points/2) {
                    if (Y_axis[k] > y_pivot) {
                        while (Y_axis[j] > y_pivot) {
                            j--;
                        }
                        if (k < j) {
                            swap(X_axis, k, j);
                            swap(Y_axis, k, j);
                        } 
                    }
                    k++;
                }
            }
            x_cut = 0;
        }
        quadrants *= 2;
    }
}

void quadrants_coord(int x_num, int y_num, int offset, int num_quadrants){
    /**************************
    # This function calculates coordinates of each quadrants in a certain processor.
    # Input:
    #   myid: the id of processor.
    #   x_num: the number of processor firstly partition the data.
    #   y_num: after first partition, each processor then partition data into y_num parts, each part 
    #          is sent to a processor.
    # Output:
    **************************/


    // # of points in each quadrants
    int block_points = NUM_POINTS / num_quadrants;
    // # of quadrants in each processors
    int proc_quad = num_quadrants / (x_num * y_num);
    // # of points in each processor
    int proc_points = NUM_POINTS / (x_num * y_num);
    // local variables
    // each quadrants have a min_x, min_y, max_x, max_y
    int min_x[proc_quad];
    int min_y[proc_quad];
    int max_x[proc_quad];
    int max_y[proc_quad];


    int i = 0;

    // for each quadrant
    for (i=0; i<proc_quad; i++){
        // initialize min_x, max_x, min_y, max_y
        min_x[i] = X_axis[offset + i * block_points];
        max_x[i] = X_axis[offset + i * block_points];
        min_y[i] = Y_axis[offset + i * block_points];
        max_y[i] = Y_axis[offset + i * block_points];

        for (j = offset; j < offset + block_points; j++){
            if(min_x[i] > X_axis[j]){
                min_x[i] = X_axis[j];
            }
            if(max_x[i] < X_axis[j]){
                max_x[i] = X_axis[j];
            }
            if(min_y[i] > Y_axis[j]){
                min_y[i] = Y_axis[j];
            }
            if(max_y[i] < Y_axis[j]){
                max_y[i] = Y_axis[j];
            }            
        }
    }

    for(i = 0; i < proc_quad; i++){
        printf("Quadrant: top_left: (%d, %d) top_right: (%d, %d) bottom_left: (%d, %d) bottom_right: (%d, %d) \n", min_x[i],max_y[i],max_x[i], max_y[i], min_x[i], min_y[i], max_x[i], min_y[i]);
    }    

}




double compute_cost(int offset, int x_num, int y_num, int num_quadrants){
    /**************************
    # This function compute the cost
    # Input:
    #   offset: the address of data this processor received minus base address(X_axis, Y_axis).
    #   x_num: the number of processor firstly partition the data.
    #   y_num: after first partition, each processor then partition data into y_num parts, each part 
    #          is sent to a processor.
    #   num_quadrants: the number of quadrants we should implement.
    # Output:
    #   
    **************************/
    int num_each_procs = NUM_POINTS/ (x_num * y_num);
    int qudrants_each_procs = num_quadrants / (x_num * y_num);
    int points =  num_each_procs / qudrants_each_procs;

    double local_cost = 0;

    int j;
    for (j = 0; j < points - 1; j++) {
        int k;
        for (k = j+1; k < points; k++) {
            int x1 = j;
            int x2 = k;
            int y1 = j;
            int y2 = k;
            int i;
            for (i = 0; i < qudrants_each_procs; i++){
                double diff_x = abs(X_axis[offset + i * points + x1] - X_axis[offset + i * points +  x2]);
                double diff_y = abs(Y_axis[offset + i * points + y1] - Y_axis[offset + i * points +  y2]);
                local_cost += sqrt((double)diff_x * diff_x + diff_y * diff_y);
            }
        }
        return local_cost;
   }  
}

unsigned int find_kth(unsigned int *x, int n, int k, unsigned int *y) {
    /**************************
    # This function finds the kth smallest value in x
    # Input:
    #   x: the base address of array X_axis
    #   n: the size of array
    #   k: the kth smallest value 
    #   y: the base address of array Y_axis
    # Output:
    #   the kth smallest value in x
    **************************/
    if (n == 1 && k == 0) return x[0];
    //divide the array into n/5 subarrays of 5 elements and find the medians for each on the subarrays
    int m = (n + 4)/5;
    unsigned int *medians =  (unsigned int *)malloc(m * sizeof(int));

    int i1 = 0;
    int j0 = 0;
    int j1 = 0;
    for (i1 = 0; i1<m; i1++) {
        if (5*i1 + 4 < n) { 
            unsigned int *w = x + 5*i1; 
            unsigned int *w1 = y + 5*i1; 
           
            for (j0=0; j0<3; j0++) {
                int jmin = j0; 
                for (j1=j0+1; j1<5; j1++) {
                    if (w[j1] < w[jmin]) 
                        jmin = j1;
                }
                swap(w, j0, jmin);
                swap(w1, j0, jmin);
            }
           
            medians[i1] = w[2];
           
        } else {
          
            medians[i1] = x[5*i1]; 
        }
    }

    //find the median of medians
    int pivot = find_kth(medians, m, m/2, medians);
    free(medians);
  

   
    for (i1=0; i1<n; i1++) {
        if (x[i1] == pivot) {
            swap(x, i1, n-1);
            swap(y, i1, n-1);
            break;
        }
    }

 
    int store = 0;
    for (i1=0; i1<n-1; i1++) {
        if (x[i1] < pivot) {
            swap(x, i1, store);
            swap(y, i1, store);
            store++;
        }
    }

    swap(x, store, n-1);
    swap(y, store, n-1);
    
    if (store == k) {
      
        return pivot;
    } else if (store > k) {
      
        return find_kth(x, store, k, y);
    } else {
        return find_kth(x+store+1, n-store-1, k-store-1, y+store+1);
    }
}


void swap(unsigned int array[], int i, int j) {
    /**************************
    # This function swap the value of array[i] and array[j]
    # Input:
    #   array[]: the array
    #   i: array[i]
    #   j: array[j] 
    # Output:
    **************************/
    int temp = array[i];
    array[i] = array[j];
    array[j] = temp;
}

int main(argc,argv)
int argc;
char *argv[];
{
    int num_quadrants;
    int  namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];

    /*Time Variables*/
    double startwtime = 0.0, endwtime;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPI_Get_processor_name(processor_name,&namelen);
    
    if (argc != 2)
    {
        fprintf (stderr, "Usage: recursive_bisection <#of quadrants>\n");
        MPI_Finalize();
        exit (0);
    }
    
    fprintf (stderr,"Process %d on %s\n", myid, processor_name);
    
    num_quadrants = atoi (argv[1]);
    
    if (myid == 0)
        fprintf (stdout, "Extracting %d quadrants with %d processors \n", num_quadrants, numprocs);
        
        if (myid == 0)
        {
            int i;
            
            srand (10000);
            
            for (i = 0; i < NUM_POINTS; i++)
                X_axis[i] = (unsigned int)rand();

            
            for (i = 0; i < NUM_POINTS; i++)
                Y_axis[i] = (unsigned int)rand();

            //start timer at process 0
            printf("\nComputing Parallely Using MPI.\n\n");
            startwtime = MPI_Wtime();
        }
    
    // compute x_d, y_d
    struct retval r = compute_dimenssion(num_quadrants, numprocs);
    int x_d = power(2, r.x);
    int y_d = power(2, r.y);
    if (myid == 0){
        printf("\n processor mesh  =  %d * %d\n", x_d,y_d);
    }
    if (myid < x_d){
    seperate_data_x_axis(x_d, y_d);
    }
    

    MPI_Barrier(MPI_COMM_WORLD);
    if (!(num_quadrants < numprocs && myid >= num_quadrants)){
        seperate_data_y_axis(x_d, y_d);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // partition in each processor
    int interval  = NUM_POINTS / x_d ;
    int points = NUM_POINTS / (x_d * y_d);
    int offset = (myid % x_d) * interval + (myid / x_d) * points;
    if (!(num_quadrants < numprocs && myid >= num_quadrants)){
        partition(offset, x_d, y_d, num_quadrants);
    }

    if (myid == 0){
    printf("\n Programm has finished seperating data\n", x_d,y_d);
    }

    MPI_Barrier(MPI_COMM_WORLD);
 
    if (!(num_quadrants < numprocs && myid >= num_quadrants)){
        quadrants_coord(x_d, y_d, offset, num_quadrants);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    double local_cost = 0.0;
    if (!(num_quadrants <= numprocs && myid >= num_quadrants)){
        local_cost = compute_cost(offset, x_d, y_d, num_quadrants);
    }
    
    MPI_Reduce(&local_cost, &global_cost, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


        
    // MPI_Barrier(MPI_COMM_WORLD);

    if (myid == 0) {
        // printf("\n processor mesh  =  %d * %d\n", x_d,y_d);
        //end timer at process 0
        endwtime = MPI_Wtime();
        printf("\nelapsed time = %f\n", endwtime - startwtime);
        printf("\nTotal cost:  %lf \n", global_cost);

    }
    
    MPI_Finalize();

    return 0;
}

