/* ***************************************************************

	 This code is subject to the MIT License (MIT)

	 Copyright (c) 2013 Kelly Black, Candace Liu, Elizabeth Sweeney, Lindong Zhou

	 Permission is hereby granted, free of charge, to any person
	 obtaining a copy of this software and associated documentation
	 files (the "Software"), to deal in the Software without
	 restriction, including without limitation the rights to use, copy,
	 modify, merge, publish, distribute, sublicense, and/or sell copies
	 of the Software, and to permit persons to whom the Software is
	 furnished to do so, subject to the following conditions:

	 The above copyright notice and this permission notice shall be
	 included in all copies or substantial portions of the Software.

	 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
	 EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
	 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
	 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
	 BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
	 ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
	 CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
	 SOFTWARE.

	 *****************************************************************

	 Code to run a Monte Carlo simulation of the systems of equations
	 representing the population of deer in one area and the associated
	 funds in an account to offset insurance liabilities.

	 ***************************************************************** */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* Helper functions to set command line options. */
#include <getopt.h>
#include <string.h>

#define DEFAULT_FILE "trial.csv"
#define DEBUG

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





int main(int argc,char **argv)
{

  /* Define the basic run time variables. */
	 double dt;
	 double sdt;
	 double t;
   double initialTime =  0.0;
   double finalTime   = 10.0;
   double m[2];
	 double dW[2];
   int lupe;
	 int timeLupe;
	 int numberIters     = 100000;
   int numberTimeSteps = 1000000;

	 /* define some book keeping variables. */
	 double stochasticIntegral;  // The integral used for the sol. to the pop eqn.
	 double z;                   // The transformed (linear) sol. to the pop eqn.
	 double W;                   // The random walk.

	 /* Define the estimated parameters for the problem. */
	 double r1   = 0.053;    // Deer max reproduction rate
	 double h    = 0.16;     // Harvest rate of the deer
	 double F    = 28000.0;  // Carrying capacity of the deer.
	 double rho  = 0.04;     // Bond fund rate of growth: log(1+rate); 
	 double beta = 9.0;      // cost due to deer collisions .003*3000 */
	 double g    = 0.05;     // Net target rate of growth of the fund.

	 /* Scaled parameter values. */
	 double rtilde;      // scaled growth rate
	 double ftilde;      // scaled carrying capacity
	 double a;           // exp  exponent for sol. to deep eqn.
	 double g0;          // int. constant for deer pop. solution.



   /* define the parameters ranges*/
   double Pmin     = 430000.0;
   double alphaMin = 0.0;

   double Pmax     = 530000.0;
   double alphaMax = 0.15;

   double deltaP;
   double deltaAlpha;

   int numP     = 1000;
   int numAlpha = 1000;
   int lupeP,lupeAlpha;

   /* define the parameters */
   double P;
   double alpha;

	 /* Statistical values */
	 float sumX  = 0.0;
	 float sumX2 = 0.0;
	 float sumM  = 0.0;
	 float sumM2 = 0.0;


   /* Define the output parameters. */
   char outFile[1024];
   strcpy(outFile,DEFAULT_FILE);

  /* File stuff */
  FILE *fp;

  /* Set the step values for the parameters. */
  deltaP     = calcDelta(Pmin,Pmax,numP);
  deltaAlpha = calcDelta(alphaMin,alphaMax,numAlpha);

  /* Set the number of iterations used in the main loop. */
  dt  = ((finalTime-initialTime)/((float)numberTimeSteps));
	sdt = sqrt(dt);

#ifdef DEBUG
  printf("Starting iteration. %d iterations.\n",numberTimeSteps);
#endif

	/* Open the output file and print out the header. */
  fp = fopen(outFile,"w");
  fprintf(fp,"time,P,alpha,x,m,sumx,sumx2,summ,summ2,N\n");
	//fprintf(fp,"time,x,m\n");

	/* Set the seed for the random number generator. */
	srand48(time(NULL));


	/* Go through and run the simulations for all possible values of the parameters. */
  for(lupeP=0;lupeP<=numP;++lupeP)
    {
      P = Pmin + deltaP*((double)lupeP);

      for(lupeAlpha=0;lupeAlpha<=numAlpha;++lupeAlpha) 
        {
          alpha = alphaMin + deltaAlpha*((double)lupeAlpha);


#ifdef DEBUG
					/* print a notice */
					printf("%f,%f,%f\n",
								 dt*((float)numberTimeSteps),P,alpha);
#endif

					/* Determine and set the scaled parameters */
					rtilde = r1-h;                     // scaled growth rate
					ftilde = (rtilde/r1)*F;            // scaled carrying capacity
					a      = rtilde-0.5*(alpha*alpha); // exp  exponent for sol. to deep eqn.
					g0     = 0.5*alpha*alpha/a;        // int. constant for deer pop. solution.
					// todo - keep track of 1/a. 
					// keep track of exp(*) or factor it out appropriately?

					/* Start the loop for the multiple simulations. */
					sumX  = 0.0;
					sumX2 = 0.0;
					sumM  = 0.0;
					sumM2 = 0.0;
					for(lupe=0;lupe<numberIters;++lupe)
						{
							/* set the initial conditions. */
							W    = 0.0;
							m[0] = ftilde;
							m[1] = (P-beta*ftilde)/(g-rho);
							stochasticIntegral = 0.0;

							for(timeLupe=0;timeLupe<numberTimeSteps;++timeLupe)
								{
									/* Set the time step. */
									t = ((double)timeLupe)*dt;

									/* Calc. two normally distributed random numbers */
									if(timeLupe%2==0)
										randNormal(dW); // calc. a new set of random numbers.
									else
										dW[0] = dW[1];    // shift the 2nd number into the first slot.
									dW[0] *= sdt;       // scale the change in W to have the proper variance.

									// Update the integral and then update the population and fund balance.
									stochasticIntegral += 
										exp(a*t+alpha*W)*(dW[0] + 0.5*alpha*(dW[0]*dW[0]-dt));
									z = rtilde/a - 
										exp(-a*t-alpha*W)*(g0 + ((alpha*rtilde)/a)*stochasticIntegral); 
									m[0] = ftilde/z;
									m[1] += (rho*m[1]+P-beta*m[0])*dt - beta*m[0]*dW[0] 
										- 0.5*alpha*beta*m[0]*(dW[0]*dW[0]-dt);

									//fprintf(fp,"%f,%f,%f\n",t,m[0],m[1]);

									W += dW[0];
								}
							//fprintf(fp,"%f,%f,%f\n",t,m[0],m[1]);

							// Update the tally used for the statistical ensemble
							sumX  += m[0];
							sumX2 += m[0]*m[0];
							sumM  += m[1]*1.0E-1;
							sumM2 += m[1]*m[1]*1.0E-2;
						}

					//  fprintf(fp,"time,P,alpha,x,m,sumx,sumx2,summ,summ2,N\n");
					fprintf(fp,"%f,%f,%f,%f,%f,%f,%f,%f,%f,%d\n",
									dt*((float)numberTimeSteps),
									P,alpha,m[0],m[1],sumX,sumX2,sumM,sumM2,numberIters);


				}

		}

	fclose(fp);
	return(0);
}

