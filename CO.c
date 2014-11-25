#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*       CO vibrational and rotational transition frequency program */
/*       initial writing:  1994 by James M Anderson                 */
/*       uses data from Autihier et al., 1993, J Mol. Spec.         */
/*       at present, each run of the program results in one         */
/*       frequency output                                           */
/*                                                                  */
/*       to compile, use an ANSI C compiler and add the math library*/
/*       for instance,                                              */
/*       cc -o CO.program CO.c -lm                                  */



extern int	main(void);
extern double	Energy(int v_i, int J_i, double massC, double massO);




int main(void)
{

	int number;
	double m_C, m_O, E, c, wl;
	int v, J;
	
	printf("Program to calculate the line frequencies of\n");
	printf("CO transitions\n\n");
	printf("Enter the isotope of Carbon\n");
	scanf("%d",&number);
	
	c = 2.99792458E+10;

	if(number == 9) m_C = 9.031039;
	else if(number == 10) m_C = 10.016856;
	else if(number == 11) m_C = 11.011433;
	else if(number == 12) m_C = 12.000000;
	else if(number == 13) m_C = 13.003355;
	else if(number == 14) m_C = 14.003242;
	else if(number == 15) m_C = 15.010599;
	else m_C = (double)(abs(number)) + 1.0E-30;
	
	printf("Enter the isotope of Oxygen\n");
	scanf("%d",&number);
	
	if(number == 14) m_O = 14.008595;
	else if(number == 15) m_O = 15.003065;
	else if(number == 16) m_O = 15.994915;
	else if(number == 17) m_O = 16.999131;
	else if(number == 18) m_O = 17.999160;
	else if(number == 19) m_O = 19.003577;
	else if(number == 20) m_O = 20.004076;
	else m_O = (double)(abs(number)) + 1.0E-30;
	
	printf("Enter the initial vibrational state\n");
	scanf("%d",&number);
	v = abs(number);
	
	printf("Enter the initial rotational state\n");
	scanf("%d",&number);
	J = abs(number);
	
	E = Energy(v, J, m_C, m_O);
	
	printf("Enter the final vibrational state\n");
	scanf("%d",&number);
	v = abs(number);
	
	printf("Enter the final rotational state\n");
	scanf("%d",&number);
	J = abs(number);
	
	E = E - Energy(v, J, m_C, m_O);
	
	/* since constants in MHz, convert to Hz */
	E = E * 1.0E6;
	/* calculate the wavelength in um */
	wl = c / E * 1E+4;
	
	printf("\nThe frequency is %.10E Hz\n", E);
	printf("\nThe wavelength is %.10E um\n", wl);

	return 0;
}








double Energy(int v_i, int J_i, double massC, double massO)
{

	double m_e = 5.485797958676E-4;
	double v, J, mu, Y;
	int k, l;
	int max0 = 10, max1 = 8, max2 = 5, max3 = 3, max4 = 3, max5 = 2, max6 = 1;
	double E = 0.0;
	
	/*   Use U[l] to represent U(k,l), deltaC(k,l), and deltaO(k,l) */
	/*	in U[l][k][0], U[l][k][1], and U[l][k][2] respectively */
	double U[7][10][3] ={	0.0,				0.0, 			0.0,
							1.703231207129E8,	7.00930E-1,		-1.71324E-1,
							-2.73128342148E6,	4.2815E-1,		-9.1035E-1,
							5.61274901E3,		-1.22868E1,		-3.534E0,
							9.584848E1,			0.0,			0.0,
							9.98618E-1,			0.0,			0.0,
							-2.58278E-2,		0.0,			0.0,
							-1.81869E-2,		0.0,			0.0,
							6.92125E-4,			0.0,			0.0,
							-1.21496E-5,		0.0,			0.0,
							
							3.97029002897E5,	-2.05455E0, 	-2.09818E0,
							-9.422413539E3,		-1.2452E0,		-3.0206E0,
							1.0271995E0,		0.0,			0.0,
							-9.230681E-2,		0.0,			0.0,
							4.520158E-2,		0.0,			0.0,
							-3.577782E-3,		0.0,			0.0,
							8.680382E-5,		0.0,			0.0,
							-3.686248E-6,		0.0,			0.0,
							0.0,				0.0,			0.0,
							0.0,				0.0,			0.0,
							
							-8.62892397E0,		-7.353E0, 		1.350E0,
							1.686743E-3,		0.0,			0.0,
							-1.190827E-3,		0.0,			0.0,
							6.150632E-5,		0.0,			0.0,
							-6.912644E-6,		0.0,			0.0,
							0.0,				0.0,			0.0,
							0.0,				0.0,			0.0,
							0.0,				0.0,			0.0,
							0.0,				0.0,			0.0,
							0.0,				0.0,			0.0,
							
							5.242352E-5,		0.0, 			0.0,
							1.631888E-5,		0.0,			0.0,
							-8.11727E-8,		0.0,			0.0,
							0.0,				0.0,			0.0,
							0.0,				0.0,			0.0,
							0.0,				0.0,			0.0,
							0.0,				0.0,			0.0,
							0.0,				0.0,			0.0,
							0.0,				0.0,			0.0,
							0.0,				0.0,			0.0,
							
							1.445283E-8,		0.0, 			0.0,
							-7.867657E-8,		0.0,			0.0,
							-2.29581E-12,		0.0,			0.0,
							0.0,				0.0,			0.0,
							0.0,				0.0,			0.0,
							0.0,				0.0,			0.0,
							0.0,				0.0,			0.0,
							0.0,				0.0,			0.0,
							0.0,				0.0,			0.0,
							0.0,				0.0,			0.0,
							
							-2.068987E-14,		0.0, 			0.0,
							-7.04332E-15,		0.0,			0.0,
							0.0,				0.0,			0.0,
							0.0,				0.0,			0.0,
							0.0,				0.0,			0.0,
							0.0,				0.0,			0.0,
							0.0,				0.0,			0.0,
							0.0,				0.0,			0.0,
							0.0,				0.0,			0.0,
							0.0,				0.0,			0.0,
							
							-4.73222E-18,		0.0, 			0.0,
							0.0,				0.0,			0.0,
							0.0,				0.0,			0.0,
							0.0,				0.0,			0.0,
							0.0,				0.0,			0.0,
							0.0,				0.0,			0.0,
							0.0,				0.0,			0.0,
							0.0,				0.0,			0.0,
							0.0,				0.0,			0.0,
							0.0,				0.0,			0.0
					};
	
	v = (double)(v_i);
	J = (double)(J_i);
	if(J_i == 0) J = 1.0E-30;
	mu = 1.0 / (1.0/massC + 1.0/massO);
	
	for(l=0; l<7; l++) {
		for(k=0; k<max0; k++) {
			Y = U[l][k][0] * pow(mu,-0.5*(double)(k+2*l));
			Y = Y * (1.0 + m_e * (U[l][k][1]/massC + U[l][k][2]/massO));
			E = E + Y * pow(v + 0.5, (double)(k)) * pow(J*(J+1.0), (double)(l));
		}
	}
	
	return E;
}
