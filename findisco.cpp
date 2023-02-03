/* Calculate ISCO radius */
void find_isco()
{
    long double d2Veff, d2Veff2, E_var, Lz_var, Omega_var, denom, den, den_r, den_rr, num, num_r, num_rr;
    long double m[4][4],dmdr[4][4],dmdr2[4][4];
    long double l, j, r, d2Veff_last = 0, d2Veff_last2 = 1000, rin;
    long double jmin, jmax, lmin, lmax, rmin, factor, rstep;

    rstep = 0.1;
    rmin = 1.1;
    jmax = 15.0;
    jmin = rmin;
    factor = 1.0e2;

    for(j=jmax; j>jmin && j>rmin; j-=rstep/factor)
    {
        r = j;

		/* Calculate 2nd-derivative of the effective potential - d2Veff */
        metric(r, Pi/2., m);
        metric_rderivatives(r, Pi/2., dmdr);
        metric_r2derivatives(r, Pi/2., dmdr2);

        Omega_var = (-dmdr[0][3] + sqrt(dmdr[0][3]*dmdr[0][3] - dmdr[0][0]*dmdr[3][3])) / dmdr[3][3];
        denom = sqrt(-(m[0][0] + 2.0*m[0][3]*Omega_var + m[3][3]*Omega_var*Omega_var));
        E_var = -(m[0][0] + m[0][3]*Omega_var) / denom;
        Lz_var = (m[0][3] + m[3][3]*Omega_var) / denom;

        den = m[0][3]*m[0][3]-m[0][0]*m[3][3];
        den_r = 2.0*m[0][3]*dmdr[0][3]-dmdr[0][0]*m[3][3]-m[0][0]*dmdr[3][3];
        den_rr = 2.0*dmdr[0][3]*dmdr[0][3]+2.0*m[0][3]*dmdr2[0][3]-dmdr2[0][0]*m[3][3]-2.0*dmdr[0][0]*dmdr[3][3]-m[0][0]*dmdr2[3][3];
        num = E_var*E_var*m[3][3]+2.0*E_var*Lz_var*m[0][3]+Lz_var*Lz_var*m[0][0];
        num_r = E_var*E_var*dmdr[3][3]+2.0*E_var*Lz_var*dmdr[0][3]+Lz_var*Lz_var*dmdr[0][0];
        num_rr = E_var*E_var*dmdr2[3][3]+2.0*E_var*Lz_var*dmdr2[0][3]+Lz_var*Lz_var*dmdr2[0][0];

        d2Veff = (num_rr + (-2.0*num_r*den_r - num*den_rr + 2.0*num*den_r*den_r/den) / den) / den;

		/* Search for when d2Veff flips sign */
        if( (d2Veff < 0.0 && d2Veff_last > 0.0) || (d2Veff > 0.0 && d2Veff_last < 0.0) )
        {
            rin = j;
            lmax = rin + rstep/factor;
            lmin = rin - rstep/factor;
			jmin = lmin;
            factor *= 100.0;
            d2Veff_last2 = 1000;

			/* Zoom in on where d2Veff flips sign to increase ISCO accuracy */
            while(factor < 1.0e10)
            {
                for(l=lmax; l>lmin && l>rmin; l-=rstep/factor)
                {
                    r = l;

					/* Calculate 2nd-derivative of the effective potential - d2Veff */
                    metric(r, Pi/2., m);
                    metric_rderivatives(r, Pi/2., dmdr);
                    metric_r2derivatives(r, Pi/2., dmdr2);

                    Omega_var = (-dmdr[0][3] + sqrt(dmdr[0][3]*dmdr[0][3] - dmdr[0][0]*dmdr[3][3])) / dmdr[3][3];
                    denom = sqrt(-(m[0][0] + 2.0*m[0][3]*Omega_var + m[3][3]*Omega_var*Omega_var));
                    E_var = -(m[0][0] + m[0][3]*Omega_var) / denom;
                    Lz_var = (m[0][3] + m[3][3]*Omega_var) / denom;

                    den = m[0][3]*m[0][3]-m[0][0]*m[3][3];
                    den_r = 2.0*m[0][3]*dmdr[0][3]-dmdr[0][0]*m[3][3]-m[0][0]*dmdr[3][3];
                    den_rr = 2.0*dmdr[0][3]*dmdr[0][3]+2.0*m[0][3]*dmdr2[0][3]-dmdr2[0][0]*m[3][3]-2.0*dmdr[0][0]*dmdr[3][3]-m[0][0]*dmdr2[3][3];
                    num = E_var*E_var*m[3][3]+2.0*E_var*Lz_var*m[0][3]+Lz_var*Lz_var*m[0][0];
                    num_r = E_var*E_var*dmdr[3][3]+2.0*E_var*Lz_var*dmdr[0][3]+Lz_var*Lz_var*dmdr[0][0];
                    num_rr = E_var*E_var*dmdr2[3][3]+2.0*E_var*Lz_var*dmdr2[0][3]+Lz_var*Lz_var*dmdr2[0][0];

                    d2Veff2 = (num_rr + (-2.0*num_r*den_r - num*den_rr + 2.0*num*den_r*den_r/den) / den) / den;

                    if(fabs(d2Veff2)<fabs(d2Veff_last2))
                        rin = l;

                    d2Veff_last2 = d2Veff2;
                }

                lmax = rin + rstep/factor;
                lmin = rin - rstep/factor;
                factor *= 100.0;
                d2Veff_last2 = 1000;
            }

            factor = 1.0e2;
        }

        d2Veff_last = d2Veff;
  	}

    isco = rin;
}
