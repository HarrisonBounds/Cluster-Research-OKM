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


typedef struct
{
 double **data;   /* Attribute data [N x D]                 */
 int num_points;  /* Number of points {> 0}                 */
 int num_dims;    /* Number of attributes {> 0}             */
 int num_categories; /* Number of categories {>= 0}               */
 int *category;      /* Class membership of each point [N x 1] */
} Data_Set;


#define MAX_CHAR 32
#define MAX_DIST 195075

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


/*Online K-Means Algorithm*/
void
okm(const Data_Set* data_set, Data_Set* clusters, const int numClusters, const double lr_exp)
{
  int minIndex, mySize, rand_data_point;
  double minDist, dist, learningRate, sse;

  int *size = (int *) malloc ( numClusters * sizeof ( int ) );

  memset( size, 0, numClusters * sizeof( int )); //reset the sizes

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

    mySize = ++size[minIndex];
    
    learningRate = pow(mySize, -lr_exp);
    
    for (int k = 0; k < data_set->num_dims; k++)
    {
      clusters->data[minIndex][k] += learningRate * ( data_set->data[rand_data_point][k] - clusters->data[minIndex][k]);
    } 

  } 

  free( size );
}

Data_Set* kmeanspp(const Data_Set *data_set, const int numClusters)
{
    double dist;
    double *distances = (double *) malloc ( data_set->num_points * sizeof ( double ) );
    double dist_sum = 0.0;
    int i, j;
    Data_Set *centers = alloc_data_set(numClusters, data_set->num_dims, 0);

    //Select random data to use as initial center for K-Means++ clustering
    std::random_device rd;
    std::mt19937 gen(rd());
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

        distances[i] = dist;
        dist_sum += dist;

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
            free( distances );
            return centers;
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
rand_sel(const Data_Set *data_set, const int numClusters)
{
    int rand_int;
    Data_Set *center = alloc_data_set(numClusters, data_set->num_dims, 0);
    std::random_device rd;
    std::mt19937 gen(rd());
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
maximin(const Data_Set *data_set, const int numClusters)
	{	
    
    double *d = ( double * ) malloc ( data_set->num_points * sizeof ( double ) );
    double *sums = ( double * ) malloc ( data_set->num_dims * sizeof ( double ) );
    Data_Set *centers = alloc_data_set(numClusters, data_set->num_dims, 0);
    int next_center;
    double dist, maxDist;
    
    for(int k = 0; k < data_set->num_dims; k++) sums[k] = 0.0;

    /*Select the first center arbitrarily*/
		for(int i = 0; i < data_set->num_points; i++)
    {
      /*Set distances to 'infinity'*/
      d[i] = DBL_MAX;
			for(int k = 0; k < data_set->num_dims; k++)
      {
        sums[k] += data_set->data[i][k];
      }
		}

    /* First center is given by the centroid of the data set */
    for(int k = 0; k < data_set->num_dims; k++)
    {
      centers->data[0][k] = sums[k] / data_set->num_points;
    }
    
    
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
    free ( sums );
    return centers;

}

Data_Set*
kmeanspp_okm(const Data_Set *data_set, const int numClusters, const double lr_exp)
{
  Data_Set *centers;

  centers = kmeanspp(data_set, numClusters);
  okm(data_set, centers, numClusters, lr_exp);

  return centers;
}

Data_Set*
rand_sel_okm(const Data_Set *data_set, const int numClusters, const double lr_exp)
{
  Data_Set *centers;

  centers = rand_sel(data_set, numClusters);
  okm(data_set, centers, numClusters, lr_exp);

  return centers;
}

Data_Set*
maximin_okm(const Data_Set *data_set, const int numClusters, const double lr_exp)
{
  Data_Set *centers;

  centers = maximin(data_set, numClusters);
  okm(data_set, centers, numClusters, lr_exp);

  return centers;
}

Data_Set*
rand_sel_bkm(const Data_Set *data_set, const int numClusters, int* numIters, double* sse)
{
  Data_Set *centers;

  centers = rand_sel(data_set, numClusters);
  bkm(data_set, centers, numClusters, numIters, sse);

  return centers;
}

Data_Set*
maximin_bkm(const Data_Set *data_set, const int numClusters, int* numIters, double* sse)
{
  Data_Set *centers;

  centers = maximin(data_set, numClusters);
  bkm(data_set, centers, numClusters, numIters, sse);

  return centers;
}

Data_Set*
kmeanspp_bkm(const Data_Set *data_set, const int numClusters, int* numIters, double* sse)
{
  Data_Set *centers;

  centers = kmeanspp(data_set, numClusters);
  bkm(data_set, centers, numClusters, numIters, sse);

  return centers;
}

Data_Set*
rand_sel_okm_bkm(const Data_Set *data_set, const int numClusters, int* numIters, double* sse, const double lr_exp)
{
  Data_Set *centers;

  centers = rand_sel(data_set, numClusters);
  okm(data_set, centers, numClusters, lr_exp);
  bkm(data_set, centers, numClusters, numIters, sse);

  return centers;
}

Data_Set*
maximin_okm_bkm(const Data_Set *data_set, const int numClusters, int* numIters, double* sse, const double lr_exp)
{
  Data_Set *centers;

  centers = maximin(data_set, numClusters);
  okm(data_set, centers, numClusters, lr_exp);
  bkm(data_set, centers, numClusters, numIters, sse);

  return centers;
}

Data_Set*
kmeanspp_okm_bkm(const Data_Set *data_set, const int numClusters, int* numIters, double* sse, const double lr_exp)
{
  Data_Set *centers;

  centers = kmeanspp(data_set, numClusters);
  okm(data_set, centers, numClusters, lr_exp);
  bkm(data_set, centers, numClusters, numIters, sse);

  return centers;
}


int 
main(int argc, char *argv[])
{
   const char* filenames[] = { "data/ecoli.txt", "data/glass.txt", "data/ionosphere.txt", "data/iris_bezdek.txt",  "data/landsat.txt", "data/letter_recognition.txt", "data/segmentation.txt", "data/vehicle.txt", "data/wine.txt", "data/yeast.txt"};
    int filenamesLength = sizeof(filenames) / sizeof(filenames[0]);
    const int numClusters[] = { 8, 6, 2, 3, 6, 26, 7, 4, 3, 10 };
    Data_Set *initCenters;
    Data_Set* centers;
    int rand_sel_bkm_numIters, maximin_bkm_numIters, kmeanspp_bkm_numIters, rand_sel_okm_bkm_numIters, maximin_okm_bkm_numIters, kmeanspp_okm_bkm_numIters;
    double rand_sel_bkm_sse, maximin_bkm_sse, kmeanspp_bkm_sse, rand_sel_okm_bkm_sse, maximin_okm_bkm_sse, kmeanspp_okm_bkm_sse;
    double rand_sel_okm_sse, maximin_okm_sse, kmeanspp_okm_sse;
    const double lr_exp = 0.5;

    srand ( time ( NULL ) );

    // Data_Set* data_set = load_data_set(filenames[3]);

    // centers = rand_sel_bkm(data_set, numClusters[3], &rand_sel_bkm_numIters, &rand_sel_bkm_sse);
    // free_data_set(centers);

    // printf("\nSSE (rand_sel + bkm) = %g\n", rand_sel_bkm_sse);

    // centers = maximin_bkm(data_set, numClusters[3], &maximin_bkm_numIters, &maximin_bkm_sse);
    // free_data_set(centers);

    // printf("\nSSE (maximin + bkm) = %g\n", maximin_bkm_sse);

    // centers = kmeanspp_bkm(data_set, numClusters[3], &kmeanspp_bkm_numIters, &kmeanspp_bkm_sse);
    // free_data_set(centers);

    // printf("\nSSE (kmeanspp + bkm) = %g\n", kmeanspp_bkm_sse);


    for (int i = 0; i < filenamesLength; i++){
        printf("Filename: %s, Number of Clusters: %d\n", filenames[i], numClusters[i]);
        Data_Set* data_set = load_data_set(filenames[i]);

        centers = rand_sel_okm(data_set, numClusters[i], lr_exp);
        rand_sel_okm_sse = calc_sse(data_set, centers, numClusters[i]);
        free_data_set(centers);

        centers = maximin_okm(data_set, numClusters[i], lr_exp);
        maximin_okm_sse = calc_sse(data_set, centers, numClusters[i]);
        free_data_set(centers);

        centers = kmeanspp_okm(data_set, numClusters[i], lr_exp);
        kmeanspp_okm_sse = calc_sse(data_set, centers, numClusters[i]);
        free_data_set(centers);

        centers = rand_sel_bkm(data_set, numClusters[i], &rand_sel_bkm_numIters, &rand_sel_bkm_sse);
        free_data_set(centers);

        centers = maximin_bkm(data_set, numClusters[i], &maximin_bkm_numIters, &maximin_bkm_sse);
        free_data_set(centers);

        centers = kmeanspp_bkm(data_set, numClusters[i], &kmeanspp_bkm_numIters, &kmeanspp_bkm_sse);
        free_data_set(centers);

        centers = rand_sel_okm_bkm(data_set, numClusters[i], &rand_sel_okm_bkm_numIters, &rand_sel_okm_bkm_sse, lr_exp);
        free_data_set(centers);

        centers = maximin_okm_bkm(data_set, numClusters[i], &maximin_okm_bkm_numIters, &maximin_okm_bkm_sse, lr_exp);
        free_data_set(centers);

        centers = kmeanspp_okm_bkm(data_set, numClusters[i], &kmeanspp_okm_bkm_numIters, &kmeanspp_okm_bkm_sse, lr_exp);
        free_data_set(centers);

        printf("SSE (rand_sel + okm) = %g : SSE (maximin + okm) = %g : SSE (kmeanspp + okm) = %g\n", rand_sel_okm_sse, maximin_okm_sse, kmeanspp_okm_sse);
        printf("SSE (rand_sel + bkm) = %g [%d iters]: SSE (maximin + bkm) = %g [%d iters] : SSE (kmeanspp + bkm) = %g [%d iters]\n", rand_sel_bkm_sse, rand_sel_bkm_numIters, maximin_bkm_sse, maximin_bkm_numIters, kmeanspp_bkm_sse, kmeanspp_bkm_numIters);
        printf("SSE (rand_sel + okm + bkm) = %g [%d iters]: SSE (maximin + okm + bkm) = %g [%d iters]: SSE (kmeanspp + okm + bkm) = %g [%d iters]\n", rand_sel_okm_bkm_sse, rand_sel_okm_bkm_numIters, maximin_okm_bkm_sse, maximin_okm_bkm_numIters, kmeanspp_okm_bkm_sse, kmeanspp_okm_bkm_numIters);
        
        free_data_set(data_set);
        printf("Done processing file: %s\n\n", filenames[i]);
    }

    printf("FINAL");

    return 0;
}

