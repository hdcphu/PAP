#include <stdio.h>
#include <omp.h>

#include <x86intrin.h>

#define NBEXPERIMENTS 7
static long long unsigned int experiments[NBEXPERIMENTS];

#define NBTHREADS 2
// #define NBCHUNKS 4

typedef char BOOLEAN;
#define TRUE 1
#define FALSE 0
/*
  bubble sort -- sequential, parallel --
*/

static unsigned int N;

typedef int *array_int;

/* the X array will be sorted  */

static array_int X;

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

/*
  test if T is properly sorted
 */
int is_sorted(array_int T)
{
  register int i;

  for (i = 1; i < N; i++)
  {
    /* test designed specifically for our usecase */
    if (T[i - 1] + 1 != T[i])
    {
      printf("*** error at %d \n", i);
      return 0;
    }
  }
  return 1;
}

void sequential_bubble_sort(int *T, const int size)
{
  /* TODO: sequential implementation of bubble sort */
  BOOLEAN swapped = FALSE;
  int temp;
  do
  {
    swapped = FALSE;
    for (int i = 0; i < size; i++)
    {
      if (T[i] > T[i + 1])
      {
        temp = T[i];
        T[i] = T[i + 1];
        T[i + 1] = temp;
        swapped = TRUE;
      }
    }
  } while (swapped == TRUE);

  return;
}


void parallel_bubble_sort(int *T, const int size)
{
  /* TODO: parallel implementation of bubble sort */

  
  unsigned int NBCHUNKS = 2;
  unsigned int size_of_chunk = size / NBCHUNKS;

  register int swapped;
  register unsigned int i;
  register unsigned int j;
  register unsigned int start;
  // register
  unsigned int end;
  register int t;
  // #pragma omp parallel shared(T) private(i, j, t, start, end, swapped)
  do
  {
    swapped = FALSE;
#if 1
#pragma omp parallel for schedule(static, NBCHUNKS) private(i, j, t, start, end) reduction(|| : swapped)
    for (i = 0; i < NBCHUNKS; i++)
    {
      start = i * size_of_chunk;
      end = start + size_of_chunk;

      for (j = start; j < end; j++)
      {
        if (T[j] > T[j + 1])
        {
          t = T[j];
          T[j] = T[j + 1];
          T[j + 1] = t;
          swapped = TRUE;
        }
      }
    }
#else
#pragma omp parallel for schedule(static, size_of_chunk) private(i, j, t) shared(swapped)
    for (j = 0; j < size_of_chunk; j++)
    {
      for (i = j; i < size; i += NBCHUNKS)
      {
        // printf(" %i", T[i + 1]);
        if (T[i] > T[i + 1])
        {
          // printf("swapped at T[%i] = %i\n", i + 1, T[i + 1]);
          t = T[i];
          T[i] = T[i + 1];
          T[i + 1] = t;
          swapped = TRUE;
        }
      }
    }
#endif
    // for (int a = 0; a < N; a++)
    // {
    //   printf("%d ", T[a]);
    // }
    // printf("\n");

    // #pragma omp parallel for schedule(static, NBCHUNKS - 1) private(i, t, start, end) shared(swapped)
    end = size_of_chunk;
    for (; end < size; end += size_of_chunk)
    {
      if (T[end - 1] > T[end])
      {
        t = T[end];
        T[end] = T[end - 1];
        T[end - 1] = t;
        swapped = TRUE;
        // printf("swapped at T[%i] = %i\n", end, T[end]);
        // for (int a = 0; a < N; a++)
        // {
        //   printf("%d ", T[a]);
        // }
        // printf("\n");
      }
    }
  } while (swapped == TRUE);
  // printf("--------- XONG ----\n");
  return;
}

int main(int argc, char **argv)
{
  unsigned long long int start, end, residu;
  unsigned long long int av;
  unsigned int exp;

  /* the program takes one parameter N which is the size of the array to
     be sorted. The array will have size 2^N */
  if (argc != 2)
  {
    fprintf(stderr, "bubble N \n");
    exit(-1);
  }

  N = 1 << (atoi(argv[1]));
  X = (int *)malloc(N * sizeof(int));

  printf("--> Sorting an array of size %u\n", N);

  start = _rdtsc();
  end = _rdtsc();
  residu = end - start;

  for (exp = 0; exp < NBEXPERIMENTS; exp++)
  {
    init_array(X);

    start = _rdtsc();

    sequential_bubble_sort(X, N);

    end = _rdtsc();
    experiments[exp] = end - start;

    /* verifying that X is properly sorted */
    if (!is_sorted(X))
    {
      fprintf(stderr, "ERROR: the array is not properly sorted\n");
      // print_array(X);
      exit(-1);
    }
  }

  av = average(experiments);

  printf("\n bubble serial \t\t\t %Ld cycles\n\n", av - residu);
  // print_array(X);

  // #pragma omp parallel num_threads(NBTHREADS)
  for (exp = 0; exp < NBEXPERIMENTS; exp++)
  {
    init_array(X);
    start = _rdtsc();

    parallel_bubble_sort(X, N);

    end = _rdtsc();
    experiments[exp] = end - start;

    /* verifying that X is properly sorted */
    if (!is_sorted(X))
    {
      fprintf(stderr, "ERROR: the array is not properly sorted\n");
      print_array(X);
      exit(-1);
    }
  }

  av = average(experiments);
  printf("\n bubble parallel \t %Ld cycles\n\n", av - residu);

  // print_array (X) ;
}
