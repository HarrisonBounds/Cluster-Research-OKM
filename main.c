#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <float.h>


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
    double dist = 0.0;

    for (int i = 0; i < num_dims; i++)
    {
        dist += (vec_b[i] - vec_a[i]) * (vec_b[i] - vec_a[i]);
    }

    return dist;
}

/* Jancey Algorithm */
void
lkm(const Data_Set* data_set, Data_Set* clusters, const int numClusters, int* numIters, double* mse)
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
		*mse = 0.0;
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

			*mse += minDist;
		}
    
    
    
    //Update via batch k-means
    for (int j = 0; j < numClusters; j++) 
    {
      
      mySize = temp[j].num_points;

      if (mySize != 0) 
      {
        for (int k = 0; k < data_set->num_dims; k++)
        {
            clusters->data[j][k] = temp->data[j][k] / mySize; 
        }       
        
                
      }
    }

		(*numIters)++;

		//cout << "Iteration " << *numIters << ": SSE = " << *sse << " [" << "# changes = " << numChanges << "]" << endl;

	} while (numChanges != 0);

	free_data_set ( temp );
  free( size );
  free ( member );

}

Data_Set*
maximin(const Data_Set *data_set, const int numClusters)
	{	
    
    double *d = ( double * ) malloc ( data_set->num_points * sizeof ( double ) );
    double *sums = ( double * ) malloc ( data_set->num_dims * sizeof ( double ) );
    Data_Set *center = alloc_data_set(numClusters, data_set->num_dims, 0);
    int next_center;
    double dist, maxDist;

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
      center->data[0][k] = sums[k] / data_set->num_points;
    }
    
    
    /*Calculate the remaining centers*/
    for (int j = 0 + 1; j < numClusters; j++)
    {
      maxDist = -DBL_MAX;
      next_center = 0;

      for (int i = 0; i < data_set->num_points; i++)
      {
        
        dist = sqr_euc_dist(center->data[j-1], data_set->data[i], data_set->num_dims);
       

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

      center->data[j] = data_set->data[next_center];

    }

    free( d );
    free ( sums );
    return center;

}


int 
main(int argc, char *argv[])
{
    const char* filenames[] = { "data/ecoli.txt", "data/glass.txt", "data/iris_bezdek.txt",  "data/yeast.txt", "data/landsat.txt", "data/letter_recognition.txt", "data/segmentation.txt", "data/vehicle.txt", "data/wine.txt"};
    int filenamesLength = sizeof(filenames) / sizeof(filenames[0]);
    const int numClusters = 64;
    Data_Set *initCenters;
    int iters;
    double mse;

    for (int i = 0; i < filenamesLength; i++){
        Data_Set* data_set = load_data_set(filenames[i]);
        initCenters = maximin(data_set, numClusters);
        lkm(data_set, initCenters, numClusters, &iters, &mse);
        free_data_set(data_set);
        printf("Done processing file: %s\n", filenames[i]);
    }

    printf("FINAL");

    return 0;
}

