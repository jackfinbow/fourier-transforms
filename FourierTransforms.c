#include <stdio.h>
#include <math.h>


/* define struct of type complex
 * contains real and imaginary parts, of type double
 */
typedef struct complex_t
{
  double Re;
  double Im;
} complex;

// function h1(t) of type complex
complex h1(double t_k)
{
  complex h;
  h.Re = cos(t_k) + cos(5*(t_k));  //real part of h1
  h.Im = sin(t_k) + sin(5*(t_k));  //imaginary part of h1
  return h;
}

// function h2(t) of type complex
complex h2(double t_k)
{
  complex h;
  h.Re = exp(-(t_k - M_PI) * (t_k - M_PI) / 2);  //real part of h2
  h.Im = 0.0;  //imaginary part of h2
  return h;
}

// exponential used in discrete fourier transform
complex e(int n, int k, int N)
{
  complex e;
  e.Re = cos(2 * M_PI * n * k / N);  //real part of exponential
  e.Im = -sin(2 * M_PI * n * k / N);  //imaginary part of exponential
  return e;
}

// sampling function
void sample(complex (*p_func)(double), complex *p_array, int funcNumber, int N)
{
  int k;  // summation integer
  double t_k; // t_k = (k * T) / N
  complex h;
  FILE *fp; // file pointer

  // create text file h.txt and open for writing
  fp = fopen("h.txt", "w");

  /* sample h(t) function N times from k = 0 to k = N - 1
   * N is number of sample values
   */
  for(k = 0; k < N; k++)
  {
    t_k = k * 2 * M_PI / N;  // define t_k
    h = p_func(t_k);  // set complex h to function pointer

    // print real and imaginary parts of h, and t_k, to file
    fprintf(fp, "%.2f, ", h.Re);
    fprintf(fp, "%.2f, ", h.Im);
    fprintf(fp, "%.2f\n", t_k);

    // assign same real and imaginary parts of h to arrays for use in DFT
    (*(p_array + k)).Re = h.Re;
    (*(p_array + k)).Im = h.Im;
  }

  fclose(fp);  // close h.txt

  // if h(t) is h1(t), rename file to h1.txt
  if(funcNumber == 1)
  {
    rename("h.txt", "h1.txt");
  }

  // if h(t) is h2(t), rename file to h2.txt
  if(funcNumber == 2)
  {
    rename("h.txt", "h2.txt");
  }
}

// function to load in h3(t) data
void load_h3(complex *p_array)
{
  int N = 200;  // N, number of samples
  FILE *fp;  // file pointer
  int k_[N];  // define 2 local arrays for k and t_k to read in data from file
  double t_k_[N];

  fp = fopen("h3.txt", "r");  //open h3.txt for reading
  // loop through all lines of the file
  for(int j = 0; j < N; j++)
  {
    // scan each line of the text file into the ith element of the arrays, delimiter ", "
    fscanf(fp, "%d, %lf, %lf, %lf", &k_[j], &t_k_[j], &(*(p_array + j)).Re, &(*(p_array + j)).Im);
  }
  // pointers are used to save Re & Im to arrays in main
}

// discrete fourier transform (DFT) function
void dft(complex *p_array1, complex (*p_exp)(int, int, int), complex *p_array2, int funcNumber, int N)
{
  double t_k, w_n; // t_k & w_n for use in functions, for time and frequency
  complex e, H;  // local variables of type complex. H is DFT of h
 
  for(int n = 0; n < N; n++)  // loop n from n = 0 to n = N - 1
  {
    for(int k = 0; k < N; k++)  // for each n value, sum up values over all k from k = 0 to k = N - 1
    {
      t_k = k * 2 * M_PI / N;
      // set complex e equal to the pointer to the exponential function
      e = p_exp(n, k, N);

      H.Re = H.Re + ((*(p_array1 + k)).Re * e.Re) - ((*(p_array1 + k)).Im * e.Im);  // Im * Im = Re (i^2 = -1)
      H.Im = H.Im + ((*(p_array1 + k)).Re * e.Im) + ((*(p_array1 + k)).Im * e.Re);
    }  // use pointers to arrays from sample function to apply DFT and sum up over all k values

    //once k = N - 1 reached, save Re and Im part of H to output arrays in main, using pointers
    (*(p_array2 + n)).Re = H.Re;
    (*(p_array2 + n)).Im = H.Im;
    H.Re = 0; H.Im = 0;  // reset Re & Im parts of H to zero for next n value

    // if function called is h1(t) or h2(t)
    if(funcNumber == 1 || funcNumber == 2)
    {
      w_n = n * 2 * M_PI / N;
      printf("For w%d = %.2f, H%d(w) = %.2f + i%.2f\n", n, w_n, n, (*(p_array2 + n)).Re, (*(p_array2 + n)).Im);
    }  // print values of w and H(w) to the screen for all n values
  }
}

// inverse discrete fourier transform (IDFT) function
void idft(complex(*p_exp)(int, int, int), complex *p_array, int funcNumber, int skipValue, int N)
{
  double t_k;
  complex e, h_;  // h_ (h') is the IDFT of H
  complex max_H[200] = {0};  // local array to save 4 values of H_3 with largest amplitude for use in IDFT
  FILE *fp;

  fp = fopen("h_.txt", "w");  // create file h_.txt and open for writing
  
  if(funcNumber == 3)  // if function called is h3(t)
  {
    int maxn; // maximum n
    double amp;
    double maxAmp = 0;  // initially set max amplitude value to zero

    // loop m to find the 4 values of largest amplitude
    for(int m = 0; m < 4; m++)
    {
      // for each m value, cycle through every element of array produced in DFT function
      for(int n = 0; n < N; n++)
      {
        //calculate amplitude for each element of array by square rooting sum of Re part squared and Im part squared
        amp = sqrt(((*(p_array + n)).Re * (*(p_array + n)).Re) + ((*(p_array + n)).Im * (*(p_array + n)).Im));
        
        // updates max n and max amplitude if larger than current stored value
        if(amp > maxAmp)
        {
          maxAmp = amp;
          maxn = n;
        }
      }

      // save Re & Im part of nth element of array to local array
      max_H[maxn].Re = (*(p_array + maxn)).Re;
      max_H[maxn].Im = (*(p_array + maxn)).Im;

      // set nth element to zero after saving to new array so that this value is not found on next cycle through array
      (*(p_array + maxn)).Re = 0;
      (*(p_array + maxn)).Im = 0;
      maxAmp = 0;
      maxn = 0;
    }
  }

  for(int k = 0; k < N; k++)
  {
    if(funcNumber == 3) // if function called is h3(t)
    {
      /* use new local max H array with 4 values of greatest amplitude, and all other values zero,
       * to apply IDFT & sum up over all n values
       */
      for(int n = 0; n < N; n++)
      {
        // set complex e equal to pointer to exponential function
        e = p_exp(n, k, N);

        h_.Re = h_.Re + (max_H[n].Re * e.Re) + (max_H[n].Im * e.Im);
        h_.Im = h_.Im - (max_H[n].Re * e.Im) + (max_H[n].Im * e.Re);
      }
         
    }

    else  // if function called is h1(t) or h2(t)
    {
      for(int n = 0; n < N; n++)
      {
        e = p_exp(n, k, N);
      
        if(n == skipValue)  // if n value equals the skip value
        {
          continue;  // increment n value to n value after skip value, & return to condition part of for loop
        }

        h_.Re = h_.Re + ((*(p_array + n)).Re * e.Re) + ((*(p_array + n)).Im * e.Im);
        h_.Im = h_.Im - ((*(p_array + n)).Re * e.Im) + ((*(p_array + n)).Im * e.Re);
      }
      // use pointers to arrays from DFT function to apply IDFT & sum up over all n values
    }

    // once all n or m values have been summed up, divide Re & Im parts by N
    h_.Re = h_.Re / N;
    h_.Im = h_.Im / N;
    t_k = k * 2 * M_PI / N;
    // print Re & Im parts of h'(t), & t_k to file
    fprintf(fp, "%.2f, ", h_.Re);
    fprintf(fp, "%.2f, ", h_.Im);
    fprintf(fp, "%.2f\n", t_k);
    h_.Re = 0; h_.Im = 0;
  }

  fclose(fp);

  if(funcNumber == 1)
  {
    rename("h_.txt", "h1_.txt");  // if h'(t) is h1'(t), rename file to h1_.txt
  }

  if(funcNumber == 2)
  {
    rename("h_.txt", "h2_.txt");  // if h'(t) is h2'(t), rename file to h2_.txt
  }

  if(funcNumber == 3)
  {
    rename("h_.txt", "h3_.txt");  // if h'(t) is h3'(t), rename file to h3_.txt
  }
}

int main()
{ 
  int N = 100;  // number of samples

  // complex function pointers to h1(t), h2(t) and exponential function
  complex (*p_h1)(double);
  complex (*p_h2)(double);
  complex (*p_e)(int, int, int);

  // set address the pointers point to, as the address of the functions
  p_h1 = &h1;
  p_h2 = &h2;
  p_e = &e;

  // define complex arrays for use in sampling & DFT functions
  complex h1_n[N], h2_n[N];
  // define complex arrays for use in DFT & IDFT function
  complex H1_n[N], H2_n[N];

  //define complex pointers for the 4 arrays
  complex *p_h1;
  complex *p_h2;
  complex *p_H1;
  complex *p_H2;
  /* set address the pointers point to as the address of the arrays
   * the pointer points to the first element in the array
   */
  p_h1 = h1_n;
  p_h2 = h2_n;
  p_H1 = H1_n;
  p_H2 = H2_n;

  // sampling h1
  sample(p_h1, p_h1, 1, N);

  // sampling h2
  sample(p_h2, p_h2, 2, N);

  // discrete fourier transform of h1
  printf("For the first function:\n");

  dft(p_h1, p_e, p_H1, 1, N);

  printf("\n");

  // discrete fourier transform of h2
  printf("For the second function:\n");

  dft(p_h2, p_e, p_H2, 2, N);

  // inverse discrete fourier transform of H1
  idft(p_e, p_H1, 1, 1, N);

  // inverse discrete fourier transform of H2
  idft(p_e, p_H2, 2, 0, N);

  N = 200;  // set number of samples to 200

  complex h3_n[N];  // complex array for use in sampling & DFT functions
  complex H3_n[N];  // complex arrays for use in DFT & IDFT function
  complex *p_h3;  // complex pointers for the 2 arrays
  complex *p_H3;
  p_h3 = h3_n;  // set address the pointers point to as the address of the arrays
  p_H3 = H3_n;

  // load h3 data
  load_h3(p_h3);

  // discrete fourier transform of h3
  dft(p_h3, p_e, p_H3, 3, N);

  // inverse discrete fourier transform of H3
  idft(p_e, p_H3, 3, 0, N);

  return 0;
}
