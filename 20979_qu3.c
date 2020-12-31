#include <stdio.h>
#include <math.h>


/*define struct of type complex*/
typedef struct complex_t
{
  double Re;  //complex contains real and imaginary parts, of type double
  double Im;
} complex;

/*define h1(t) as function of type complex*/
complex h_1(double t_k)
{
  complex h;
  h.Re = cos(t_k) + cos(5*(t_k));  //real part of h1
  h.Im = sin(t_k) + sin(5*(t_k));  //imaginary part of h1
  return h;
}

/*define h2(t) as function of type complex*/
complex h_2(double t_k)
{
  complex h;
  h.Re = exp(-(t_k - M_PI) * (t_k - M_PI) / 2);  //real part of h2
  h.Im = 0.0;  //imaginary part of h2
  return h;
}

/*define exponential used in discrete fourier transform*/
complex e(int n, int k, int N)
{
  complex e;
  e.Re = cos(2 * M_PI * n * k / N);  //real part of exponential
  e.Im = -sin(2 * M_PI * n * k / N);  //imaginary part of exponential
  return e;
}

/*define sampling function*/
void sample(complex (*p_func)(double), complex *p_array, int i, int N)
{
  int k;
  double t_k;
  complex h;
  FILE *fp;

  fp = fopen("h.txt", "w");  //create text file h.txt and open for writing

  for(k=0; k<N; k++)  //sample h(t) function N times from k=0 to k=N-1
  {
    t_k = k * 2 * M_PI / N;  //define t_k
    h = p_func(t_k);  //set complex h to function pointer

    fprintf(fp, "%.2f, ", h.Re);  //print real part of h to file
    fprintf(fp, "%.2f, ", h.Im);  //print imaginary part of h to file
    fprintf(fp, "%.2f\n", t_k);  //print t_k values to file

    (*(p_array + k)).Re = h.Re;  //assign same real and imaginary parts
    (*(p_array + k)).Im = h.Im;  //of h to arrays for use in DFT
  }

  fclose(fp);  //close h.txt

  if(i==1)
  {
    rename("h.txt", "h1.txt");  //if h(t) is h1(t), rename file to h1.txt
  }

  if(i==2)
  {
    rename("h.txt", "h2.txt");  //if h(t) is h2(t), rename file to h2.txt
  }
}

/*define function to load in h3(t) data*/
void load_h_3(complex *p_array)
{
  int j; int N = 200;  //define variable for the for loop and N, number of samples
  FILE *fp;  //define file pointer
  int k_[N];  //define 2 local arrays for k and t_k to read in data from file
  double tk_[N];

  fp = fopen("h3.txt", "r");  //open h3.txt for reading

    for(j=0; j<N; j++)  //loop through all lines of the file
    {
      fscanf(fp, "%d, %lf, %lf, %lf", &k_[j], &tk_[j], &(*(p_array + j)).Re, &(*(p_array + j)).Im);
    }  //scan each line of the text file into the ith element of the arrays, delimiter ", "
}  //pointers are used to save Re & Im to arrays in main

/*define discrete fourier transform function*/
void dft(complex *p_array1, complex (*p_exp)(int, int, int), complex *p_array2, int i, int N)
{
  int n, k;  //define integers for the for loops
  double t_k, w_n; //define t_k & w_n for use in functions
  complex e, H;  //define local variables of type complex
 
  for(n=0; n<N; n++)  //loop n from n=0 to n=N-1
  {
    for(k=0; k<N; k++)  //for each n value, sum up values over all k from k=0 to k=N-1
    {
      t_k = k * 2 * M_PI / N;
      e = p_exp(n, k, N);  //set complex e equal to the pointer to the exponential function

      H.Re = H.Re + ((*(p_array1 + k)).Re * e.Re) - ((*(p_array1 + k)).Im * e.Im);  //Im * Im = Re (i^2 = -1)
      H.Im = H.Im + ((*(p_array1 + k)).Re * e.Im) + ((*(p_array1 + k)).Im * e.Re);
    }  //use pointers to arrays from sample function to apply DFT and sum up over all k values

    (*(p_array2 + n)).Re = H.Re;  //once k=N-1 reached, save Re and Im part of H
    (*(p_array2 + n)).Im = H.Im;  //to output arrays in main, using pointers
    H.Re = 0; H.Im = 0;  //reset Re & Im parts of H to zero for next n value

    if(i==1 || i==2)  //if function called is h1(t) or h2(t)
    {
    w_n = n * 2 * M_PI / N;
    printf("For w%d = %.2f, H%d(w) = %.2f + i%.2f\n", n, w_n, n, (*(p_array2 + n)).Re, (*(p_array2 + n)).Im);
    }  //print values of w and H(w) to the screen for all n values
  }
}

/*define inverse discrete fourier transform function*/
void idft(complex(*p_exp)(int, int, int), complex *p_array, int i, int j, int N)
{
  int n, k, m;  //define variable for the for loops
  double t_k;  //define t_k
  complex e, h_;  //define local variable of type complex
  complex max_H[200] = {0};  //define local array to save 4 values of H3 with largest amplitude for use in IDFT
  FILE *fp;  //define general file pointer

  fp = fopen("h_.txt", "w");  //create file h_.txt and open for writing
  
  if(i==3)  //if function called is h3(t)
  {
    int maxn;
    double amp;
    double maxamp = 0;  //initially set max amplitude value to zero

    for(m=0; m<4; m++)  //loop m from m=0 to m=3 to find the 4 values of largest amplitude
    {
      for(n=0; n<N; n++)  //for each m value, cycle through every element of array produced in DFT function
      {
        amp = sqrt(((*(p_array + n)).Re * (*(p_array + n)).Re) + ((*(p_array + n)).Im * (*(p_array + n)).Im));
        //calculate amplitude for each element of array by square rooting sum of Re part squared and Im part squared

        if(amp > maxamp)
        {
          maxamp = amp;  //if calculated value of amplitude, redefine maxamp variable to this value
          maxn = n;  //save n value, corresponding to element number, for max amplitude
        }
      }

      max_H[maxn].Re = (*(p_array + maxn)).Re;  //save Re & Im part of nth element of array to local array
      max_H[maxn].Im = (*(p_array + maxn)).Im;

      (*(p_array + maxn)).Re = 0;  //set nth element to zero after saving to new array
      (*(p_array + maxn)).Im = 0;  //so that this value is not found on next cycle through array
      maxamp = 0;  //reset maxamp and maxn to zero after each successive max amplitude found
      maxn = 0;
    }
  }

  for(k=0; k<N; k++)  //loop k from k=0 to k=N-1
  {
    if(i==3)  //if function called is h3(t)
    {
      for(n=0; n<N; n++)  //loop through all n values, from n=0 to n=N-1
      {
        e = p_exp(n, k, N);  //set complex e equal to pointer to exponential function

        h_.Re = h_.Re + (max_H[n].Re * e.Re) + (max_H[n].Im * e.Im);
        h_.Im = h_.Im - (max_H[n].Re * e.Im) + (max_H[n].Im * e.Re);
      }  //use new local max H array with 4 values of greatest amplitude, and all other values zero,
         //to apply IDFT & sum up over all n values
    }

    else  //if function called is h1(t) or h2(t)
    {
      for(n=0; n<N; n++)  //loop through all n values, from n=0 to n=N-1
      {
        e = p_exp(n, k, N);  //set complex e equal to pointer to exponential function
      
        if(n==j)  //if n value equals the skip value
        {
          continue;  //increment n value to n value after skip value, & return to condition part of for loop
        }

        h_.Re = h_.Re + ((*(p_array + n)).Re * e.Re) + ((*(p_array + n)).Im * e.Im);
        h_.Im = h_.Im - ((*(p_array + n)).Re * e.Im) + ((*(p_array + n)).Im * e.Re);
      }  //use pointers to arrays from DFT function to apply IDFT & sum up over all n values
    }

    h_.Re = h_.Re/N;  //once all n or m values have been summed up, divide Re & Im parts by N
    h_.Im = h_.Im/N;
    t_k = k * 2 * M_PI / N;
    fprintf(fp, "%.2f, ", h_.Re);  //print Re & Im parts of h'(t), & t_k to file
    fprintf(fp, "%.2f, ", h_.Im);
    fprintf(fp, "%.2f\n", t_k);
    h_.Re = 0; h_.Im = 0;  //reset Re & Im parts of h' to zero for next k value
  }

  fclose(fp);  //close fp

  if(i==1)
  {
    rename("h_.txt", "h1_.txt");  //if h'(t) is h1'(t), rename file to h1_.txt
  }

  if(i==2)
  {
    rename("h_.txt", "h2_.txt");  //if h'(t) is h2'(t), rename file to h2_.txt
  }

  if(i==3)
  {
    rename("h_.txt", "h3_.txt");  //if h'(t) is h3'(t), rename file to h3_.txt
  }
}

int main()
{ 
  int N = 100;  //number of samples

  complex (*p_h_1)(double);  //define complex function pointers to h1(t) & h2(t)
  complex (*p_h_2)(double);
  complex (*p_e)(int, int, int);  //define complex function pointers to exponential function
  p_h_1 = &h_1;  //set address the pointers point to, as the address of the functions
  p_h_2 = &h_2;
  p_e = &e;

  complex h1_n[N], h2_n[N];  //define complex arrays for use in sampling & DFT functions
  complex H1_n[N], H2_n[N];  //define complex arrays for use in DFT & IDFT function

  complex *p_h1;  //define complex pointers for the 4 arrays
  complex *p_h2;
  complex *p_H1;
  complex *p_H2;
  p_h1 = h1_n;  //set address the pointers point to as the address of the arrays
  p_h2 = h2_n;  //the pointer points to the first element in the array
  p_H1 = H1_n;
  p_H2 = H2_n;

  /*sampling h1*/
  sample(p_h_1, p_h1, 1, N);

  /*sampling h2*/
  sample(p_h_2, p_h2, 2, N);

  /*discrete fourier transform of h1*/
  printf("For the first function:\n");

  dft(p_h1, p_e, p_H1, 1, N);

  printf("\n");

  /*discrete fourier transform of h2*/
  printf("For the second function:\n");

  dft(p_h2, p_e, p_H2, 2, N);

  /*inverse discrete fourier transform of H1*/
  idft(p_e, p_H1, 1, 1, N);

  /*inverse discrete fourier transform of H2*/
  idft(p_e, p_H2, 2, 0, N);

  N = 200;  //set number of samples to 200

  complex h3_n[N];  //define complex array for use in sampling & DFT functions
  complex H3_n[N];  //define complex arrays for use in DFT & IDFT function
  complex *p_h3;  //define complex pointers for the 2 arrays
  complex *p_H3;
  p_h3 = h3_n;  //set address the pointers point to as the address of the arrays
  p_H3 = H3_n;

  /*load h3 data*/
  load_h_3(p_h3);

  /*discrete fourier transform of h3*/
  dft(p_h3, p_e, p_H3, 3, N);

  /*inverse discrete fourier transform of H3*/
  idft(p_e, p_H3, 3, 0, N);

  return 0;  //return 0 to the shell upon successful run
}

