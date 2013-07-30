#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Helper functions to set command line options. */
#include <getopt.h>
#include <string.h>

#define DEFAULT_FILE "trial.csv"


/* Check to see if the order is right between min and max. */
void checkOrder(double *theMin,double *theMax)
{
  double tmp;
  if(*theMax < *theMin)
    {
      tmp     = *theMax;
      *theMax = *theMin;
      *theMin = tmp;
    }

}

/* Routine to calculate the step size for a given range and number of iterations. */
double calcDelta(double theMin,double theMax,int number)
{
  return((theMax-theMin)/((double)number));
}




/* Routine to generate two normally distributed random numbers */
inline void randNormal(double nu[]) {
  double tmp = sqrt(-2.0*log(drand48()));
  double trig = 2.0*M_PI*drand48();
  nu[0] = tmp*sin(trig);
  nu[1] = tmp*cos(trig);
}




/* The function to evaluate for the differential equation. */
inline void functionEvaluation(double t, 
															 double *x, 
															 double *result,
															 double dt,
															 double dw,
															 double P,
															 double alpha,
															 double gamma)
{
  result[0] = x[0]*(1.0-x[1]);
  result[1] = x[1]-x[0];
}

int main(int argc,char **argv)
{

  /* Define the basic run time variables. */
   double dt  = 1.0e-4;
	 double sdt = sqrt(dt);
   double t;
   double initialTime =  0.0;
   double finalTime   = 10.0;
   double x[2];
	 double dW[2];
   int lupe;
   int numberIters;
   int skipPrint = 25;    /* Set this to 4 to print out every fourth value. */

	 /* define some book keeping variables. */
	 double stochasticIntegral;  // The integral used for the sol. to the pop eqn.
	 double z;                   // The transformed (linear) sol. to the pop eqn.
	 double W;                   // The random walk.

	 /* Define the estimated parameters for the problem. */
	 double r1 = 1.7;     // Deer max reproduction rate
	 double h = .16;      // Harvest rate of the deer
	 double F = 28000;    // Carrying capacity of the deer.
	 double rho = .04;    // Bond fund rate of growth: log(1+rate); 
	 double beta = 9;     // cost due to deer collisions .003*3000 */
	 double g = .05;      // Net target rate of growth of the fund.

	 /* Scaled parameter values. */
	 double rtilde;      // scaled growth rate
	 double ftilde;      // scaled carrying capacity
	 double a;           // exp  exponent for sol. to deep eqn.
	 double b;           // exp. exponent for integration term
	 double g0;          // int. constant for deer pop. solution.



   /* define the parameters ranges*/
   double Pmin     = 2000000.0;
   double alphaMin = 0.0;
   double gammaMin = 0.0;

   double Pmax     = 250000.0;
   double alphaMax = 0.1;
   double gammaMax = 0.1;

   double deltaP;
   double deltaAlpha;
   double deltaGamma;

   int numP     = 100;
   int numAlpha = 100;
   int numGamma = 100;
   int lupeP,lupeAlpha,lupeGamma;

   /* define the parameters */
   double P;
   double alpha;
   double gamma;


   /* Define the output parameters. */
   char outFile[1024];
   strcpy(outFile,DEFAULT_FILE);

  /* File stuff */
  FILE *fp;

  /* Set the step values for the parameters. */
  deltaP     = calcDelta(Pmin,Pmax,numP);
  deltaAlpha = calcDelta(alphaMin,alphaMax,numAlpha);
  deltaGamma = calcDelta(gammaMin,gammaMax,numGamma);

  /* Set the number of iterations used in the main loop. */
  numberIters = (int)((finalTime-initialTime)/dt);

  printf("Starting iteration. %d iterations.\n",numberIters);
  fp = fopen(outFile,"w");
  fprintf(fp,"time,x,y,P,alpha,gamma\n");


  for(lupeP=0;lupeP<=numP;++lupeP)
    {
      P = Pmin + deltaP*((double)lupeP);

      for(lupeAlpha=0;lupeAlpha<=numAlpha;++lupeAlpha)
        {
          alpha = alphaMin + deltaAlpha*((double)lupeAlpha);

          for(lupeGamma=0;lupeGamma<=numGamma;++lupeGamma)
            {
              gamma = gammaMin + deltaGamma*((double)lupeGamma);

							/* set the scaled parameters */
							rtilde = r1-h;            // scaled growth rate
							ftilde = (rtilde/r1)*F;   // scaled carrying capacity
							a = (alpha*alpha)-rtilde; // exp  exponent for sol. to deep eqn.
							b = (0.5*(alpha*alpha))-a;// exp. exponent for integration term
							g0 = 1.0-(rtilde/b);      // int. constant for deer pop. solution.


              /* set the initial conditions. */
              x[0] = ftilde;
              //x[1] = y0; TODO
							stochasticIntegral = 0.0;

              /* Start the loop. */
              for(lupe=0;lupe<numberIters;++lupe)
                {
                  /* Set the time step. */
                  t = ((double)lupe)*dt;

									/* Calc. two normally distributed random numbers */
									if(lupe%2==0)
										randNormal(dW); // calc. a new set of random numbers.
									else
										dW[0] = dW[1];  // shift the 2nd number into the first slot.
									
									//stochasticIntegral += ; TODO
									//z = (rtilde/b) + g0*exp((-b*s)-(alpha*W)) + exp((-b*s)-(alpha*W)).*((-alpha*rtilde)/b).*stochasticIntegral; 
									//x[0] = ;
									
								}

						}

				}

		}

	fclose(fp);
	return(0);
}

