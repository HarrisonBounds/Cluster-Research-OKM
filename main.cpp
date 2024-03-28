#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <cstdlib>
#include <iostream>
#include <random>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <limits.h>


typedef struct
{
 double **data;   /* Attribute data [N x D]                 */
 int num_points;  /* Number of points {> 0}                 */
 int num_dims;    /* Number of attributes {> 0}             */
 int num_categories; /* Number of categories {>= 0}               */
 int *category;      /* Class membership of each point [N x 1] */
} Data_Set;

typedef struct
{
  int num_runs;
  double *init_sse;
  double *fin_sse;
  int *num_iters;
} Run;

typedef struct
{
 double min, max, mean, std;
} Dbl_Stats;

typedef struct
{
 int min, max;
 double mean, std;
} Int_Stats;


#define MAX_CHAR 32
#define MAX_DIST 195075

int POW2[] = { 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536 };

//To allocate memory for the 'double **data' variable, you need the following function:
double **
alloc_double_2d ( const int num_points, const int num_dims )
{
 int ip;
 double **data;

 if ( num_points <= 0 || num_dims <= 0 )
  {
   return NULL;
  }

 data = ( double ** ) malloc ( num_points * sizeof ( double * ) );
 data[0] = ( double * ) malloc ( num_points * num_dims * sizeof ( double ) );

 for ( ip = 0 + 1; ip < num_points; ip++ )
  {
   data[ip] = data[ip - 1] + num_dims;
  }

 return data;
}

Data_Set *
alloc_data_set ( const int num_points, const int num_dims, const int num_categories )
{
 Data_Set *data_set;

 if ( num_points <= 0 || num_dims <= 0 || num_categories < 0 )
  {
   fprintf ( stderr, "Number of points (%d) must be > 0\n", num_points );
   fprintf ( stderr, "Number of attributes (%d) must be > 0\n", num_dims );
   fprintf ( stderr, "Number of categories (%d) must be >= 0\n", num_categories );
   abort ( );
  }

 data_set = ( Data_Set * ) malloc ( sizeof ( Data_Set ) );

 data_set->num_points = num_points;
 data_set->num_dims = num_dims;
 data_set->num_categories = num_categories;

 data_set->data = alloc_double_2d ( num_points, num_dims );

 data_set->category = NULL;

 if ( 0 < num_categories )
  {
   data_set->category = ( int * ) malloc ( num_points * sizeof ( int ) );
  }

 return data_set;
}

Run *
alloc_run ( const int num_runs )
{
  Run *run;

 if ( num_runs <= 0)
  {
   fprintf ( stderr, "Number of runs (%d) must be > 0\n", num_runs );
   abort ( );
  }

 run = ( Run * ) malloc ( sizeof ( Run ) );

 run->num_runs = num_runs;

 run->init_sse = ( double * ) malloc ( num_runs * sizeof ( double ));
 if ( run->init_sse == NULL ) 
 {
    free( run );
    fprintf ( stderr, "Error allocating init_sse\n" ); 
    return NULL;
 }

 run->fin_sse = ( double * ) malloc ( num_runs * sizeof ( double ));
 if ( run->fin_sse == NULL ) 
 {
    free( run->init_sse );
    free( run );
    fprintf ( stderr, "Error allocating fin_sse\n" );  
    return NULL;
 }

 run->num_iters = ( int * ) malloc ( num_runs * sizeof ( int ));
 if ( run->num_iters == NULL )
 {
    free ( run->init_sse );
    free ( run->fin_sse );
    free ( run );
    fprintf ( stderr, "Error allocating num_iters\n" );  
    return NULL;
 }

 return run;
}

void
free_double_2d ( double **data )
{
 if ( data != NULL )
  {
   free ( data[0] );
   free ( data );
  }
}

void
free_data_set ( Data_Set *data_set )
{
 if ( data_set != NULL )
  {
   free_double_2d ( data_set->data );

   if ( data_set->category != NULL )
    {
     free ( data_set->category );
    }

   free ( data_set );
  }
}

void
free_run ( Run *run )
{
  if ( run->init_sse != NULL)
  {
    free ( run->init_sse );
  }

  if ( run->fin_sse != NULL)
  {
    free ( run->fin_sse );
  }

  if ( run->num_iters != NULL)
  {
    free ( run->num_iters );
  }

  free ( run );
}

void
print_vector ( const double *vec, const int num_dims )
{
 int id;

 for ( id = 0; id < num_dims - 1; id++ )
  {
   printf ( "%f ", vec[id] );
  }
 
 printf ( "%f\n", vec[num_dims - 1] );
}

void
print_data_set ( const Data_Set *data_set )
{
 int ip;

 for ( ip = 0; ip < data_set->num_points; ip++ )
  {
   print_vector ( data_set->data[ip], data_set->num_dims );
  }
}

Data_Set *
load_data_set ( const char *filename )
{
 bool have_cat_labels = false;
 char *buf;
 char *str;
 int ip, id;
 int num_points = 0;
 int num_dims = 0;
 int num_categories = 0;
 int buf_size;
 int label;
 int num_values;
 double *point;
 FILE *fp;
 Data_Set *data_set;

 fp = fopen ( filename, "r" );
 if ( fp == NULL )
  {
   fprintf ( stderr, "Error opening file '%s' !\n", filename );
   abort ( );
  }

 num_values = fscanf ( fp, "%d %d %d\n", &num_points, &num_dims, &num_categories );
 if ( num_values != 3 )
  {
   fprintf ( stderr, "Error reading file '%s' line 1:\n", filename );
   fprintf ( stderr, "Expected 3 integers (# points, #attributes, #categories), found %d !\n", num_values );
   abort ( );
  }

 if ( num_points <= 0 || num_dims <= 0 || num_categories < 0 )
  {
   fprintf ( stderr, "Error reading file '%s' line 1:\n", filename );
   fprintf ( stderr, "Number of points (%d) must be > 0\n", num_points );
   fprintf ( stderr, "Number of attributes (%d) must be > 0\n", num_dims );
   fprintf ( stderr, "Number of categories (%d) must be >= 0\n", num_categories );
   abort ( );
  }

 if ( 0 < num_categories )
  {
   /* Last column contains class labels */
   have_cat_labels = true;
   num_dims--;
  }

 /* Allocate */
 data_set = alloc_data_set ( num_points, num_dims, num_categories );

 buf_size = num_dims * MAX_CHAR;
 buf = ( char * ) malloc ( buf_size * sizeof ( char ) );

 for ( ip = 0; ip < num_points; ip++ )
  {
   point = data_set->data[ip];

   /* Read an entire line from the file */
   if ( fgets ( buf, buf_size, fp ) == NULL )
    {
     fprintf ( stderr, "Error reading file '%s' line %d !\n", filename, ip + 1 );
     abort ( );
    }

   /* Read the components of this point */
   for ( str = strtok ( buf, " ,;" ), id = 0;
         str != NULL && id < num_dims;
         str = strtok ( NULL, " ,;" ), id++ )
    {
     point[id] = atof ( str );
    }

   /* Line contains fewer than NUM_DIMS attributes or class label is missing */
   if ( id != num_dims || ( have_cat_labels && str == NULL ) )
    {
     fprintf ( stderr, "Error reading file '%s' line %d (expected %d values, found %d)!\n", filename, ip, num_dims, id );
     abort ( );
    }
   
   if ( have_cat_labels )
    {
     /* Last attribute is the class label */
     label = atoi ( str );

     if ( num_categories <= label || label < 0 )
      {
       fprintf ( stderr, "Class label (%d) must be in [%d,%d]!\n", label, 0, num_categories - 1 );
       abort ( );
      }

     data_set->category[ip] = label;
    }
  }

 free ( buf );

 fclose ( fp );

 return data_set;
}


void
reset_double_2d ( double **data, const int num_points, const int num_dims )
{
 memset ( data[0], 0, num_points * num_dims * sizeof ( double ) );
}

/* Reset clusters */
void
reset_data_set ( const Data_Set *data_set )
{
 /*
 memset ( data_set->data[0], 0, 
          data_set->num_points * data_set->num_dims * sizeof ( double ) );
  */
 reset_double_2d ( data_set->data, data_set->num_points, data_set->num_dims );

 if ( data_set->category != NULL )
  {
   memset ( data_set->category, 0, data_set->num_points * sizeof ( int ) );
  }
}


double 
sqr_euc_dist( const double *vec_a, const double *vec_b, const int num_dims)
{
    double delta, dist = 0.0;

    for (int i = 0; i < num_dims; i++)
    {
        delta = vec_b[i] - vec_a[i];
        dist += delta * delta;
    }
    return dist;
}

double
calc_sse(const Data_Set* data_set, const Data_Set* clusters, const int numClusters)
{
  double minDist, dist, sse = 0.0;
  int minIndex;
  for (int i = 0; i < data_set->num_points; i++)
		{
      
			minDist = DBL_MAX;
			minIndex = 0;

			for (int j = 0; j < numClusters; j++)
			{
				dist = sqr_euc_dist(clusters->data[j], data_set->data[i], data_set->num_dims); 
				if (dist < minDist)
				{
					minDist = dist;
					minIndex = j;
				}
			}
      sse += minDist;
    }
    return sse;
}

void 
data_set_compare(const Data_Set* data_set_a, Data_Set* data_set_b)
{
  for (int i = 0; i < data_set_a->num_points; i++)
  {
    for (int k = 0; k < data_set_a->num_dims; k++)
    {
      if ( fabs( data_set_a->data[i][k] - data_set_b->data[i][k] ) > 1e-12 )
      {
        printf("Point%d\n", i);
        print_vector(data_set_a->data[i], data_set_a->num_dims);
        print_vector(data_set_b->data[i], data_set_b->num_dims);
      }
    }
  }
}


/* Batch k-means Algorithm */
void
bkm(const Data_Set* data_set, Data_Set* clusters, const int numClusters, int* numIters, double* sse)
{
  
	int numChanges, minIndex, mySize;
  double minDist, dist;

  Data_Set *temp = alloc_data_set(numClusters, data_set->num_dims, 0);
  int *size = (int *) malloc ( numClusters * sizeof ( int ) );
  int *member = ( int * ) calloc ( data_set->num_points, sizeof ( int ) );


  *numIters = 0;

	do
	{
		numChanges = 0;
		*sse = 0.0;
    memset( size, 0, numClusters * sizeof( int )); //reset the sizes
    reset_data_set(temp); // Reset data set

		for (int i = 0; i < data_set->num_points; i++)
		{
      
			minDist = DBL_MAX;
			minIndex = 0;

			for (int j = 0; j < numClusters; j++)
			{
				dist = sqr_euc_dist(clusters->data[j], data_set->data[i], data_set->num_dims); 
				if (dist < minDist)
				{
					minDist = dist;
					minIndex = j;
				}
			}

      
      //Assignment
      for (int k = 0; k < data_set->num_dims; k++)
      {
        temp->data[minIndex][k] += data_set->data[i][k];
      } 
  
      size[minIndex] += 1;

      
      if (minIndex != member[i])
			{
				numChanges += 1;
				member[i] = minIndex;

			}

			*sse += minDist;
		}

		(*numIters)++;
    //printf("Iteration %d: SSE = %f: Num Changes = %d\n", *numIters, *sse, numChanges);
	  if ( numChanges == 0 ) break;
    
    
    //Update via batch k-means
    for (int j = 0; j < numClusters; j++) 
    {
      
      mySize = size[j];

      if (mySize != 0) 
      {
        for (int k = 0; k < data_set->num_dims; k++)
        {
            clusters->data[j][k] = temp->data[j][k] / mySize; 
        }       
        
                
      }
    }
	} while ( true );

	free_data_set ( temp );
  free( size );
  free ( member );

}

Data_Set*
compute_centroid(const Data_Set* data_set, const int numClusters) {
    double *sum = new double[data_set->num_dims](); // Initialize sum for each dimension
    Data_Set* centroid = alloc_data_set(numClusters, data_set->num_dims, 0);

    // Calculate sum of each dimension across all data points
    for (int i = 0; i < data_set->num_points; i++)
    {
        for (int j = 0; j < numClusters; j++)
        {
            sum[j] += data_set->data[i][j];
        }
    }

    // Compute centroid for each dimension
    for (int j = 0; j < numClusters; j++) 
    {
      for (int k = 0; k < data_set->num_dims; k++)
      {
        centroid->data[j][k] = sum[j] / data_set->num_points;
      }
    }

    delete[] sum;
    return centroid;
}

/*Incremental Batch K-Means*/

void
ibkm(const Data_Set* data_set, const int numClusters, int* numIters, double* sse)
{
  int num_splits, index;
  int numChanges, minIndex, mySize;
  double minDist, dist;
  num_splits = ( int ) ( log ( numClusters ) / log ( 2 ) + 0.5 ); /* round */

  Data_Set *temp = alloc_data_set(numClusters, data_set->num_dims, 0);
  int *size = (int *) malloc ( numClusters * sizeof ( int ) );
  int *member = ( int * ) calloc ( data_set->num_points, sizeof ( int ) );

  *numIters = 0;

  memcpy(temp->data[0], data_set->data[5], data_set->num_dims * sizeof(double)); //replace this with centroid

  /*Start iteration counter (t) = 0*/
  for ( int t = 0; t < num_splits; t++ )
  {
    for ( int n = POW2[t] - 1; n < POW2[t + 1] - 1; n++ )
    {
      for (int j = 0; j < numClusters; j++)
      {
        index = 2 * n + 1;
        temp->data[index][j] = data_set->data[j];
      }
    }
  }
  
}


/*Online K-Means Algorithm*/
void
okm(const Data_Set* data_set, Data_Set* clusters, const int numClusters, const double lr_exp)
{
  int minIndex, mySize, rand_data_point;
  double minDist, dist, learningRate, sse;

  int *size = ( int * ) calloc ( numClusters, sizeof ( int ) );

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> rand_index(0, data_set->num_points - 1);

  for (int i = 0 ; i < data_set->num_points; i++)
  {
    rand_data_point = rand_index(gen);
    
    minDist = DBL_MAX;
    minIndex = 0;
    
    for (int j = 0; j < numClusters; j++)
    {

      dist = sqr_euc_dist(clusters->data[j], data_set->data[rand_data_point], data_set->num_dims);
      if (dist < minDist)
      {
        minDist = dist;
        minIndex = j;
      }
    }

    learningRate = pow(++size[minIndex], -lr_exp);
    
    for (int k = 0; k < data_set->num_dims; k++)
    {
      clusters->data[minIndex][k] += learningRate * ( data_set->data[rand_data_point][k] - clusters->data[minIndex][k]);
    } 

  } 

  free( size );
}

void
comp_stats ( const Run *run, Dbl_Stats *init_sse_stats, Dbl_Stats *fin_sse_stats, Int_Stats *num_iters_stats)
{
  init_sse_stats->std = fin_sse_stats->std = num_iters_stats->std = 0.0;

if ( run->num_runs == 1 )
{
 init_sse_stats->min = init_sse_stats->max = init_sse_stats->mean = run->init_sse[0];
 fin_sse_stats->min = fin_sse_stats->max = fin_sse_stats->mean = run->fin_sse[0];
 num_iters_stats->min = num_iters_stats->max = num_iters_stats->mean = run->num_iters[0];
 
 return;
}

init_sse_stats->min = fin_sse_stats->min = DBL_MAX;
num_iters_stats->min = INT_MAX;

init_sse_stats->max = fin_sse_stats->max = 0.0;
num_iters_stats->max = 0;

init_sse_stats->mean = fin_sse_stats->mean = num_iters_stats->mean =  0.0;

  for(int i = 0; i < run->num_runs; i++)
  {
    //Calculate min
    if (run->init_sse[i] < init_sse_stats->min)
    {
      init_sse_stats->min = run->init_sse[i];
    }

    if (run->fin_sse[i] < fin_sse_stats->min)
    {
      fin_sse_stats->min = run->fin_sse[i];
    }

    if (run->num_iters[i] < num_iters_stats->min)
    {
      num_iters_stats->min = run->num_iters[i];
    }

    //Calculate max
    if (run->init_sse[i] > init_sse_stats->max)
    {
      init_sse_stats->max = run->init_sse[i];
    }

    if (run->fin_sse[i] > fin_sse_stats->max)
    {
      fin_sse_stats->max = run->fin_sse[i];
    }

    if (run->num_iters[i] > num_iters_stats->max)
    {
      num_iters_stats->max = run->num_iters[i];
    }

    //Calculate mean
    init_sse_stats->mean += run->init_sse[i];
    fin_sse_stats->mean += run->fin_sse[i];
    num_iters_stats->mean += run->num_iters[i];

  }

  init_sse_stats->mean /= run->num_runs;
  fin_sse_stats->mean /= run->num_runs;
  num_iters_stats->mean /= run->num_runs;

  //Standard Deviation
  for ( int i = 0; i < run->num_runs; i++)
  {
    init_sse_stats->std += (run->init_sse[i] - init_sse_stats->mean) * (run->init_sse[i] - init_sse_stats->mean);
    fin_sse_stats->std += (run->fin_sse[i] -  fin_sse_stats->mean) * (run->fin_sse[i] -  fin_sse_stats->mean);
    num_iters_stats->std += (run->num_iters[i] - num_iters_stats->mean) * (run->num_iters[i] - num_iters_stats->mean);
  }

  init_sse_stats->std = sqrt(init_sse_stats->std / (run->num_runs - 1));
  fin_sse_stats->std = sqrt(fin_sse_stats->std / (run->num_runs - 1));
  num_iters_stats->std = sqrt(num_iters_stats->std / (run->num_runs - 1));

  
}

Data_Set* kmeanspp(const Data_Set *data_set, const int numClusters, std::mt19937& gen)
{
    double dist;
    double *distances = (double *) malloc ( data_set->num_points * sizeof ( double ) );
    double dist_sum = 0.0;
    int i, j;
    Data_Set *centers = alloc_data_set(numClusters, data_set->num_dims, 0);

    //Select random data to use as initial center for K-Means++ clustering
    std::uniform_int_distribution<> rand_index(0, data_set->num_points - 1);

    //definte a stochastic threshold r_thresh
    std::uniform_real_distribution<> rand_interval(0, 1);
    double r_thresh = dist_sum * rand_interval(gen);

    //Copy the data point from data_set as the first centerS
    memcpy(centers->data[0], data_set->data[rand_index(gen)], data_set->num_dims * sizeof(double));

    /*Assign L2 distance from first center to each datum*/
    for (i = 0; i < data_set->num_points; i++)
    {
        dist = sqr_euc_dist(centers->data[0], data_set->data[i], data_set->num_dims);

        dist_sum += ( distances[i] = dist );

    }

    //Assign remaining centers
    for (j = 1; j < numClusters; j++)
    {
        //find next center
        for (i = 0; i < data_set->num_points; i++)
        {
            if (r_thresh < distances[i])
            { 
              break; 
            }
            else 
            { 
              r_thresh -= distances[i]; 
            }
        }

        memcpy(centers->data[j], data_set->data[i], data_set->num_dims * sizeof(double));

        if (j == numClusters - 1) //finished populating centers
        {
            break;
        }

        //Update minimum distances
        for (i = 0; i < data_set->num_points; i++)
        {
            dist = sqr_euc_dist(centers->data[j], data_set->data[i], data_set->num_dims);
            
            if (dist < distances[i])
            {
                dist_sum -= (distances[i] - dist);
                distances[i] = dist;
            }
        }
    }

    free ( distances );
    return centers;

}
 
Data_Set*
rand_sel(const Data_Set *data_set, const int numClusters, std::mt19937& gen)
{
    int rand_int;
    Data_Set *center = alloc_data_set(numClusters, data_set->num_dims, 0);
    // std::random_device rd;
    // std::mt19937 gen(rd());
    std::uniform_int_distribution<int> distribution(0, data_set->num_points - 1);
    

    for (int j = 0; j < numClusters; j++)
    {
     rand_int = distribution(gen);
     for(int k = 0; k < data_set->num_dims; k++) 
     {
      center->data[j][k] = data_set->data[rand_int][k];
     }
    }

    return center;
}

Data_Set*
maximin(const Data_Set *data_set, const int numClusters, std::mt19937& gen)
	{	
    
    double *d = ( double * ) malloc ( data_set->num_points * sizeof ( double ) );
    Data_Set *centers = alloc_data_set(numClusters, data_set->num_dims, 0);
    int next_center;
    double dist, maxDist;
    
    std::uniform_int_distribution<> rand_index(0, data_set->num_points - 1);

    memcpy(centers->data[0], data_set->data[rand_index(gen)], data_set->num_dims * sizeof(double));
    
    
    /*Calculate the remaining centers*/
    for (int j = 0 + 1; j < numClusters; j++)
    {
      maxDist = -DBL_MAX;
      next_center = 0;

      for (int i = 0; i < data_set->num_points; i++)
      {
        
        dist = sqr_euc_dist(centers->data[j-1], data_set->data[i], data_set->num_dims);
       

        if (dist < d[i])
        {
          d[i] = dist;
        }

        if (maxDist < d[i])
        {
          maxDist = d[i];
          next_center = i;
        }
      }

      memcpy(centers->data[j], data_set->data[next_center], data_set->num_dims * sizeof(double));

    }

    free( d );
    return centers;

}

Data_Set*
kmeanspp_okm(const Data_Set *data_set, const int numClusters, const double lr_exp, std::mt19937& gen)
{
  Data_Set *centers;

  centers = kmeanspp(data_set, numClusters, gen);
  okm(data_set, centers, numClusters, lr_exp);

  return centers;
}

Data_Set*
rand_sel_okm(const Data_Set *data_set, const int numClusters, const double lr_exp, std::mt19937& gen)
{
  Data_Set *centers;

  centers = rand_sel(data_set, numClusters, gen);
  okm(data_set, centers, numClusters, lr_exp);

  return centers;
}

Data_Set*
maximin_okm(const Data_Set *data_set, const int numClusters, const double lr_exp, std::mt19937& gen)
{
  Data_Set *centers;

  centers = maximin(data_set, numClusters, gen);
  okm(data_set, centers, numClusters, lr_exp);

  return centers;
}

Data_Set*
rand_sel_bkm(const Data_Set *data_set, const int numClusters, int* numIters, double* sse, std::mt19937& gen)
{
  Data_Set *centers;

  centers = rand_sel(data_set, numClusters, gen);
  bkm(data_set, centers, numClusters, numIters, sse);

  return centers;
}

Data_Set*
maximin_bkm(const Data_Set *data_set, const int numClusters, int* numIters, double* sse, std::mt19937& gen)
{
  Data_Set *centers;

  centers = maximin(data_set, numClusters, gen);
  bkm(data_set, centers, numClusters, numIters, sse);

  return centers;
}

Data_Set*
kmeanspp_bkm(const Data_Set *data_set, const int numClusters, int* numIters, double* sse, std::mt19937& gen)
{
  Data_Set *centers;

  centers = kmeanspp(data_set, numClusters, gen);
  bkm(data_set, centers, numClusters, numIters, sse);

  return centers;
}

Data_Set*
rand_sel_okm_bkm(const Data_Set *data_set, const int numClusters, int* numIters, double* sse, const double lr_exp, std::mt19937& gen)
{
  Data_Set *centers;

  centers = rand_sel(data_set, numClusters, gen);
  okm(data_set, centers, numClusters, lr_exp);
  bkm(data_set, centers, numClusters, numIters, sse);

  return centers;
}

Data_Set*
maximin_okm_bkm(const Data_Set *data_set, const int numClusters, int* numIters, double* sse, const double lr_exp, std::mt19937& gen)
{
  Data_Set *centers;

  centers = maximin(data_set, numClusters, gen);
  okm(data_set, centers, numClusters, lr_exp);
  bkm(data_set, centers, numClusters, numIters, sse);

  return centers;
}

Data_Set*
kmeanspp_okm_bkm(const Data_Set *data_set, const int numClusters, int* numIters, double* sse, const double lr_exp, std::mt19937& gen)
{
  Data_Set *centers;

  centers = kmeanspp(data_set, numClusters, gen);
  okm(data_set, centers, numClusters, lr_exp);
  bkm(data_set, centers, numClusters, numIters, sse);

  return centers;
}

Run*
rep_rand_sel_okm(const Data_Set *data_set, const int numClusters, const int numReps, const double lr_exp, std::mt19937& gen)
{
  Data_Set *temp_centers;
  Run *run = alloc_run(numReps);

  for (int i = 0; i < numReps; i++)
  {
    temp_centers = rand_sel(data_set, numClusters, gen);
    run->init_sse[i] = calc_sse(data_set, temp_centers, numClusters);

    okm(data_set, temp_centers, numClusters, lr_exp);
    run->fin_sse[i] = calc_sse(data_set, temp_centers, numClusters);

    //OKM always iterates once
    run->num_iters[i] = 1;

    free_data_set ( temp_centers );
   }
  
  return run;
}

Run*
rep_maximin_okm(const Data_Set *data_set, const int numClusters, const int numReps, const double lr_exp, std::mt19937& gen)
{
  Data_Set *temp_centers;
  Run *run = alloc_run(numReps);

  for (int i = 0; i < numReps; i++)
  {
    temp_centers = maximin(data_set, numClusters, gen);
    run->init_sse[i] = calc_sse(data_set, temp_centers, numClusters);

    okm(data_set, temp_centers, numClusters, lr_exp);
    run->fin_sse[i] = calc_sse(data_set, temp_centers, numClusters);

    //OKM always iterates once
    run->num_iters[i] = 1;

    free_data_set ( temp_centers );
  }

  return run;
}

Run*
rep_kmeanspp_okm(const Data_Set *data_set, const int numClusters, const int numReps, const double lr_exp, std::mt19937& gen)
{
  Data_Set *temp_centers;
  Run *run = alloc_run(numReps);

  for (int i = 0; i < numReps; i++)
  {
    temp_centers = kmeanspp(data_set, numClusters, gen);
    run->init_sse[i] = calc_sse(data_set, temp_centers, numClusters);

    okm(data_set, temp_centers, numClusters, lr_exp);
    run->fin_sse[i] = calc_sse(data_set, temp_centers, numClusters);

    //OKM always iterates once
    run->num_iters[i] = 1;

    free_data_set ( temp_centers );
  }

  return run;
}

Run*
rep_rand_sel_bkm(const Data_Set *data_set, const int numClusters, const int numReps, std::mt19937& gen)
{
  
  Data_Set *temp_centers;
  Run *run = alloc_run(numReps);
  double sse = 0.0;
  int numIters = 0;

  for (int i = 0; i < numReps; i++)
  {
    temp_centers = rand_sel(data_set, numClusters, gen);
    run->init_sse[i] = calc_sse(data_set, temp_centers, numClusters);

    bkm(data_set, temp_centers, numClusters, &numIters, &sse);
    run->fin_sse[i] = sse;

    run->num_iters[i] = numIters;

    free_data_set ( temp_centers );
  }

  return run;
}

Run*
rep_maximin_bkm(const Data_Set *data_set, const int numClusters, const int numReps, std::mt19937& gen)
{
  
  Data_Set *temp_centers;
  Run *run = alloc_run(numReps);
  double sse = 0.0;
  int numIters = 0;

  for (int i = 0; i < numReps; i++)
  {
    temp_centers = maximin(data_set, numClusters, gen);
    run->init_sse[i] = calc_sse(data_set, temp_centers, numClusters);

    bkm(data_set, temp_centers, numClusters, &numIters, &sse);
    run->fin_sse[i] = sse;

    run->num_iters[i] = numIters;

    free_data_set ( temp_centers );
  }

  return run;
}

Run*
rep_kmeanspp_bkm(const Data_Set *data_set, const int numClusters, const int numReps, std::mt19937& gen)
{
  Data_Set *temp_centers;
  Run *run = alloc_run(numReps);
  double sse = 0.0;
  int numIters = 0;

  for (int i = 0; i < numReps; i++)
  {
    temp_centers = kmeanspp(data_set, numClusters, gen);
    run->init_sse[i] = calc_sse(data_set, temp_centers, numClusters);

    bkm(data_set, temp_centers, numClusters, &numIters, &sse);
    run->fin_sse[i] = sse;

    run->num_iters[i] = numIters;

    free_data_set ( temp_centers );
  }

  return run;
}

Run*
rep_rand_sel_okm_bkm(const Data_Set *data_set, const int numClusters, const int numReps, const double lr_exp, std::mt19937& gen)
{
  Data_Set *temp_centers;
  Run *run = alloc_run(numReps);
  double sse = 0.0;
  int numIters = 0;

  for (int i = 0; i < numReps; i++)
  {
    temp_centers = rand_sel(data_set, numClusters, gen);
    okm(data_set, temp_centers, numClusters, lr_exp);

    run->init_sse[i] = calc_sse(data_set, temp_centers, numClusters);
    bkm(data_set, temp_centers, numClusters, &numIters, &sse);

    run->fin_sse[i] = sse;
    run->num_iters[i] = numIters;

    free_data_set ( temp_centers );
  }

  return run;
}

Run*
rep_maximin_okm_bkm(const Data_Set *data_set, const int numClusters, const int numReps, const double lr_exp, std::mt19937& gen)
{
  
  Data_Set *temp_centers;
  Run *run = alloc_run(numReps);
  double sse = 0.0;
  int numIters = 0;

  for (int i = 0; i < numReps; i++)
  {
    temp_centers = maximin(data_set, numClusters, gen);
    okm(data_set, temp_centers, numClusters, lr_exp);

    run->init_sse[i] = calc_sse(data_set, temp_centers, numClusters);
    bkm(data_set, temp_centers, numClusters, &numIters, &sse);

    run->fin_sse[i] = sse;
    run->num_iters[i] = numIters;

    free_data_set ( temp_centers );
  }

  return run;
}

Run*
rep_kmeanspp_okm_bkm(const Data_Set *data_set, const int numClusters, const int numReps, const double lr_exp, std::mt19937& gen)
{
  
  Data_Set *temp_centers;
  Run *run = alloc_run(numReps);
  double sse = 0.0;
  int numIters = 0;

  for (int i = 0; i < numReps; i++)
  {
    temp_centers = kmeanspp(data_set, numClusters, gen);
    okm(data_set, temp_centers, numClusters, lr_exp);
    
    run->init_sse[i] = calc_sse(data_set, temp_centers, numClusters);
    bkm(data_set, temp_centers, numClusters, &numIters, &sse);

    run->fin_sse[i] = sse;
    run->num_iters[i] = numIters;

    free_data_set ( temp_centers );
  }

  return run;
}

void print_stats(const Dbl_Stats *init_sse_stats, const Dbl_Stats *fin_sse_stats, const Int_Stats *num_iters_stats)
{
  printf("init_sse: [%.2g, %.2g], %.2g +/- %.2g ;\nfina_sse: [%.2g, %.2g], %.2g +/- %.2g ;\nnum_iter: [%d, %d], %.2g +/- %.2g\n", init_sse_stats->min, init_sse_stats->max, init_sse_stats->mean, init_sse_stats->std, fin_sse_stats->min, fin_sse_stats->max, fin_sse_stats->mean, fin_sse_stats->std, num_iters_stats->min, num_iters_stats->max, num_iters_stats->mean, num_iters_stats->std);
}



int 
main(int argc, char *argv[])
{
  const char* filenames[] = { "data/s1.txt", "data/s2.txt", "data/s3.txt", "data/s4.txt", "data/a1.txt", "data/a2.txt", "data/a3.txt", "data/dim032.txt", "data/dim064.txt", "data/dim128.txt", "data/dim256.txt", "data/dim512.txt", "data/dim1024.txt"};
  int filenamesLength = sizeof(filenames) / sizeof(filenames[0]);
  const int numClusters[] = { 15, 15, 15, 15, 20, 35, 50, 16, 16, 16, 16, 16, 16};

  Run* run;
  Dbl_Stats init_sse_stats;
  Dbl_Stats fin_sse_stats;
  Int_Stats num_iters_stats;

  const double lr_exp = 0.5;
  const int numReps = 100;

  // Initialize Mersenne Twister RNG
  std::random_device rd;
  std::mt19937 gen(rd());

  printf("# Reps = %d, lr_exp = %.2g\n", numReps, lr_exp);

  for ( int i = 0; i < filenamesLength; i++)
  {
    printf("\nFilename: %s, # Clusters: %d\n", filenames[i], numClusters[i]);
    Data_Set* data_set = load_data_set(filenames[i]);
    
    run = rep_rand_sel_okm(data_set, numClusters[i], numReps, lr_exp, gen);
    comp_stats(run, &init_sse_stats, &fin_sse_stats, &num_iters_stats);
    printf("rep_rand_sel_okm: \n");
    print_stats(&init_sse_stats, &fin_sse_stats, &num_iters_stats);
    free(run);

    run = rep_maximin_okm(data_set, numClusters[i], numReps, lr_exp, gen);
    comp_stats(run, &init_sse_stats, &fin_sse_stats, &num_iters_stats);
    printf("\nrep_maximin_okm: \n");
    print_stats(&init_sse_stats, &fin_sse_stats, &num_iters_stats);
    free(run);

    run = rep_kmeanspp_okm(data_set, numClusters[i], numReps, lr_exp, gen);
    comp_stats(run, &init_sse_stats, &fin_sse_stats, &num_iters_stats);
    printf("\nrep_kmeanspp_okm: \n");
    print_stats(&init_sse_stats, &fin_sse_stats, &num_iters_stats);
    free(run); 

    run = rep_rand_sel_bkm(data_set, numClusters[i], numReps, gen);
    comp_stats(run, &init_sse_stats, &fin_sse_stats, &num_iters_stats);
    printf("\nrep_rand_sel_bkm: \n");
    print_stats(&init_sse_stats, &fin_sse_stats, &num_iters_stats);
    free(run);
    
    run = rep_maximin_bkm(data_set, numClusters[i], numReps, gen);
    comp_stats(run, &init_sse_stats, &fin_sse_stats, &num_iters_stats);
    printf("\nrep_maximin_bkm: \n");
    print_stats(&init_sse_stats, &fin_sse_stats, &num_iters_stats);
    free(run);

    run = rep_kmeanspp_bkm(data_set, numClusters[i], numReps, gen);
    comp_stats(run, &init_sse_stats, &fin_sse_stats, &num_iters_stats);
    printf("\nrep_kmeanspp_bkm: \n");
    print_stats(&init_sse_stats, &fin_sse_stats, &num_iters_stats);
    free(run);

    run = rep_rand_sel_okm_bkm(data_set, numClusters[i], numReps, lr_exp, gen);
    comp_stats(run, &init_sse_stats, &fin_sse_stats, &num_iters_stats);
    printf("\nrep_rand_sel_okm_bkm: \n");
    print_stats(&init_sse_stats, &fin_sse_stats, &num_iters_stats);
    free(run);

    run = rep_maximin_okm_bkm(data_set, numClusters[i], numReps, lr_exp, gen);
    comp_stats(run, &init_sse_stats, &fin_sse_stats, &num_iters_stats);
    printf("\nrep_maximin_okm_bkm: \n");
    print_stats(&init_sse_stats, &fin_sse_stats, &num_iters_stats);
    free(run);

    run = rep_kmeanspp_okm_bkm(data_set, numClusters[i], numReps, lr_exp, gen);
    comp_stats(run, &init_sse_stats, &fin_sse_stats, &num_iters_stats);
    printf("\nrep_kmeanspp_okm_bkm: \n");
    print_stats(&init_sse_stats, &fin_sse_stats, &num_iters_stats);
    free(run);

    free_data_set(data_set);

  }

  printf("DONE");

  return 0;
}

