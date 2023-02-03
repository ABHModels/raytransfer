void rayprecise(long double rdisk, long double rerrtol, long double pscr, long double traced[4])
{
    long double rscr, rmiss, rhigh, rlow, rinit;

	//Set initial search values
    rinit = 5.0*rdisk;
    rscr = rinit;
    rhigh = rinit;
    rlow = 0.0;

    xyfromrphi(rhigh,pscr,rdisk);
    raytrace(xscr,yscr,traced,rdisk);

	//While missing the disk, shift the high value - rhigh
    while(traced[0]<rdisk)
    {
        if(rhigh>1.0e6) //Abort when radius gets too large, must be some error
            abort();

        rhigh *= 2.0;

        xyfromrphi(rhigh,pscr,rdisk);
        raytrace(xscr,yscr,traced,rdisk);
    }

	//Loop until radial value is found and satisfies error tolerance
    while(1)
    {
        rmiss = traced[0] - rdisk; //Variable to check error tolerance

        //printf("%Le %Le %Le %Le\n", traced[0], traced[1], traced[2], rscr);

        if(fabs(rmiss)<rerrtol) //Error tolerance satisfied
        {
            traced[3] = rscr;
            break;
        }
        else if(traced[0] < rdisk || traced[0] != traced[0]) //Value too low
            rlow = rscr;
        else if(traced[0] > rdisk) //Value too high
            rhigh = rscr;

        rscr = 0.5*(rlow+rhigh); //Next r_screen value is average of low and high values

        if((rhigh-rlow)<rerrtol/10.0) //Lower error tolerance if can't find correct value
            rerrtol *= 2.0;

        xyfromrphi(rscr,pscr,rdisk);
        raytrace(xscr,yscr,traced,rdisk);
	}
}

// Calculate x_screen and y_screen values from r_screen and phi_screen
void xyfromrphi(long double rscr, long double pscr, long double rdisk)
{
    long double h, h_eff = 0.0;

    h = 3.0*Mdl*(1.0 - sqrt(isco/rdisk))/eta; //height of disk at given radius

    if(rdisk >= isco)
        h_eff = h*sin(inc); //effective height adjusted by inclination angle

    xscr = rscr*cos(pscr);
	yscr = rscr*sin(pscr) + h_eff; //shift yscreen by effective height as the center of the image will shift
}
