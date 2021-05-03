/*! \file gene_bruit_rayleigh_scalaire.c
    \brief sous-programme generant un scalaire suivant une loi de rayleigh
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>


void gene_bruit_rayleigh_scalaire(float param_loi, float *scalaire_rayleigh_real,float *scalaire_rayleigh_imag)

{
 unsigned int g_seed;
 static int compteur_appel=-1;
 time_t temps;
 float var_uni1,var_uni2;

 compteur_appel=compteur_appel+1;
 if (compteur_appel==0)
   {
     temps=time(0);
     temps=temps-(3666*24*365*(2005-1970));
     temps=(time_t)fmod((double)temps,100000.); /* 100000 valeurs possibles pour le germe */
     g_seed=(unsigned int)temps;

     srand(g_seed);
/*      fprintf(stderr,"initialisation du germe (g_seed=%d)\n",g_seed); */
   }
 if (compteur_appel>1000000) compteur_appel=-1; /* tous les 1000000 d appels a la fonction on change le germe */

	 do{
	   var_uni1=((float)rand())/((float)RAND_MAX);
	 }while(var_uni1 == 0.);

/* variable uniforme entre 0 et 1 */
	 var_uni2=((float)rand())/((float)RAND_MAX);
/* variable gaussienne */
	 *(scalaire_rayleigh_real)=param_loi* (float)(sqrt( (-2. * log((double) (var_uni1 + 1.e-20)) )) * cos((double)(2.*M_PI*var_uni2)) ); 
	 *(scalaire_rayleigh_imag)=param_loi* (float)(sqrt( (-2. * log((double) (var_uni1 + 1.e-20)) )) * sin((double)(2.*M_PI*var_uni2)) ); 



return;
}
