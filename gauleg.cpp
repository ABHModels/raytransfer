/* Calculate 100 radial values using same method as RELCONV */
void gauleg(long double rdisk_i, long double rdisk_f, long double rdisk[])
{
    long double auxr1, auxr2, auxxm, auxxl, auxz, auxp1, auxp2, auxp3, auxpp, auxz1, auxr[imax];

	auxr1 = 1.0/sqrt(rdisk_f);
	auxr2 = 1.0/sqrt(rdisk_i);

	auxxm = 0.5*(auxr2 + auxr1);
	auxxl = 0.5*(auxr2 - auxr1);

	for (int i = 1; i <= (imax + 1.)/2.; i++) {

		auxz = cos( Pi*(i - 0.25)/(imax + 0.5) );

		do {

			auxp1 = 1.0;
			auxp2 = 0.0;

			for (int j = 1; j <= imax; j++) {

				auxp3 = auxp2;
				auxp2 = auxp1;
				auxp1 = ((2.0*j - 1.0)*auxz*auxp2 - (j - 1.0)*auxp3)/j;

			}

			auxpp = imax*(auxz*auxp1 - auxp2)/(auxz*auxz - 1.0);
			auxz1 = auxz;
			auxz  = auxz1 - auxp1/auxpp;

		} while ( fabs(auxz - auxz1) > 3.0e-14 );

		auxr[i - 1]    = auxxm - auxxl*auxz;
		auxr[imax - i] = auxxm + auxxl*auxz;

	}

	for (int i = 1; i <= imax; i++)
		rdisk[i] = 1.0/(auxr[i-1]*auxr[i-1]);

	/* Add extra radius at beginning and end so that transfer function can be computing in Python script */
	rdisk[0] = rdisk[1] + 10.0;
	rdisk[imax+1] = rdisk[imax] - 0.0001;
}
