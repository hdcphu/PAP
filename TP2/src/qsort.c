#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <omp.h>

#include <x86intrin.h>

#define NBEXPERIMENTS 7

static long long unsigned int experiments[NBEXPERIMENTS];

/* 
  quick sort -- sequential, parallel --
*/

static unsigned int N;

typedef int *array_int;

static array_int X;

void init_array(array_int T)
{
  register int i;

  for (i = 0; i < N; i++)
  {
    T[i] = N - i;
  }
}

void print_array(array_int T)
{
  register int i;

  for (i = 0; i < N; i++)
  {
    printf("%d ", T[i]);
  }
  printf("\n");
}

int is_sorted(array_int T)
{
  register int i;

  for (i = 1; i < N; i++)
  {
    /* test designed specifically for our usecase */
    if (T[i - 1] + 1 != T[i])
      return 0;
  }
  return 1;
}

long long unsigned int average(long long unsigned int *exps)
{
  unsigned int i;
  long long unsigned int s = 0;

  for (i = 2; i < (NBEXPERIMENTS - 2); i++)
  {
    s = s + exps[i];
  }

  return s / (NBEXPERIMENTS - 2);
}

static int compare(const void *x, const void *y)
{
  /* + TODO: comparison function to be used by qsort()*/

  /* cast x and y to int* before comparing */
  int v = (*((int*) x) - *((int*) y));
  return v < -1 ? -1 : v > 1 ? 1 : v ;
}

void sequential_qsort_sort(int *T, const int size)
{

  /* + TODO: sequential sorting based on libc qsort() function */

  qsort(T, size, sizeof(int), compare);

  return;
}

/* 
   Merge two sorted chunks of array T!
   The two chunks are of size size
   First chunck starts at T[0], second chunck starts at T[size]
*/
void merge(int *T, const int size)
{
  int *X = (int *)malloc(2 * size * sizeof(int));

  int i = 0;
  int j = size;
  int k = 0;

  while ((i < size) && (j < 2 * size))
  {
    if (T[i] < T[j])
    {
      X[k] = T[i];
      i = i + 1;
    }
    else
    {
      X[k] = T[j];
      j = j + 1;
    }
    k = k + 1;
  }

  if (i < size)
  {
    for (; i < size; i++, k++)
    {
      X[k] = T[i];
    }
  }
  else
  {
    for (; j < 2 * size; j++, k++)
    {
      X[k] = T[j];
    }
  }

  memcpy(T, X, 2 * size * sizeof(int));
  free(X);

  return;
}

void parallel_qsort_sort(int *T, const int size)
{

  /* + TODO: parallel sorting based on libc qsort() function +
     * sequential merging */

  int NB_CHUNKS = 8;
  register unsigned int* start;
  register unsigned int i;
  register unsigned int subsize = size / NB_CHUNKS;
#pragma omp parallel for schedule(static, NB_CHUNKS) private(i, start)
  for (i = 0; i < NB_CHUNKS; i++)
  {
    start = T + i * subsize;
    qsort(start, subsize, sizeof(int), compare);
  }
  
  unsigned int c = NB_CHUNKS;
  unsigned int s = subsize;
  unsigned int* end = T + size;
  while (c > 1)
  {    
// #pragma omp parallel for schedule(dynamic) private(s, start)
    for (start = T; start < end; start = start + s + s)
    {      
      merge(start, s);      
      // printf("%d %ui", start, s);
    }
    c = c >> 1;
    s = s << 1;
    // for (int a = 0; a < N; a++)
    // {
    //   printf("%d ", T[a]);
    // }
    // printf("\n");
  }
  // for (i = 0; i < chunks; i++)
  // {
    
  // }
  
}

void parallel_qsort_sort1(int *T, const int size)
{

  /* TODO: parallel sorting based on libc qsort() function +
     * PARALLEL merging */
}

int main(int argc, char **argv)
{
  unsigned long long int start, end, residu;
  unsigned long long int av;
  unsigned int exp;

  if (argc != 2)
  {
    fprintf(stderr, "qsort N \n");
    exit(-1);
  }

  N = 1 << (atoi(argv[1]));
  X = (int *)malloc(N * sizeof(int));

  printf("--> Sorting an array of size %u\n", N);

  // int a[3] = {5, 2, 5};
  // printf(" compare %i \n", compare(a, a + 1 ));
  // printf(" compare %i \n", compare(a + 1, a ));
  // printf(" compare %i \n", compare(a + 2, a ));

  start = _rdtsc();
  end = _rdtsc();
  residu = end - start;

  printf("sequential sorting ...\n");

  for (exp = 0; exp < NBEXPERIMENTS; exp++)
  {
    init_array(X);

    start = _rdtsc();

    sequential_qsort_sort(X, N);

    end = _rdtsc();
    experiments[exp] = end - start;

    if (!is_sorted(X))
    {
      fprintf(stderr, "ERROR: the array is not properly sorted\n");
      exit(-1);
    }
  }

  av = average(experiments);
  printf("\n qsort serial\t\t %Ld cycles\n\n", av - residu);

  printf("parallel (seq merge) ...\n");

  for (exp = 0; exp < NBEXPERIMENTS; exp++)
  {
    init_array(X);

    start = _rdtsc();

    parallel_qsort_sort(X, N);

    end = _rdtsc();
    experiments[exp] = end - start;

    if (!is_sorted(X))
    {
      fprintf(stderr, "ERROR: the array is not properly sorted\n");
      exit(-1);
    }
  }

  av = average(experiments);
  printf("\n qsort parallel (seq merge) \t %Ld cycles\n\n", av - residu);

  printf("parallel ...\n");

  for (exp = 0; exp < NBEXPERIMENTS; exp++)
  {
    init_array(X);

    start = _rdtsc();

    parallel_qsort_sort1(X, N);

    end = _rdtsc();
    experiments[exp] = end - start;

    if (!is_sorted(X))
    {
      fprintf(stderr, "ERROR: the array is not properly sorted\n");
      exit(-1);
    }
  }

  av = average(experiments);
  printf("\n qsort parallel \t %Ld cycles\n\n", av - residu);

  //   print_array (X) ;
}
