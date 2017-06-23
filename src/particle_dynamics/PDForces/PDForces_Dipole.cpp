
# include "PDForces_Dipole.hpp"



// -------------------------------------------------------------------------
// Constructor...
// -------------------------------------------------------------------------

Dipole::Dipole(const GetPot& p, vector<double>& rin, vector<double>& vin,
        vector<double>& fin, vector<double>& radin) :
    r(rin), v(vin), f(fin), rad(radin)
{
    rcut = p("PDApp/inter_particle_forces/rcut",4.0);
    rcut2 = rcut*rcut;
    K = p("PDApp/inter_particle_forces/K",50.0);
    n = p("PDApp/inter_particle_forces/n",-13.0);
    Ex = p("PDApp/inter_particle_forces/Ex",0.0);
    Ey = p("PDApp/inter_particle_forces/Ey",0.0);
    Ez = p("PDApp/inter_particle_forces/Ez",1.0);
    box[0] = p("Domain/nx",5);
    box[1] = p("Domain/ny",5);
    box[2] = p("Domain/nz",5);
    box[0] *= p("Domain/dx",1.0);
    box[1] *= p("Domain/dy",1.0);
    box[2] *= p("Domain/dz",1.0);
    Emag = sqrt(Ex*Ex + Ey*Ey + Ez*Ez);
    edir[0] = Ex/Emag;
    edir[1] = Ey/Emag;
    edir[2] = Ez/Emag;
}



// -------------------------------------------------------------------------
// Destructor...
// -------------------------------------------------------------------------

Dipole::~Dipole()
{
}



// -------------------------------------------------------------------------
// Function to calculate fij:
// -------------------------------------------------------------------------

void Dipole::fijFunc(int i, int j)
{
    //compute the squared particle distance:
    rij2 = 0.0;
    for (int k=0; k<3; k++) 
    {
        rij[k] = r[i*3+k] - r[j*3+k];
        rij[k] -= round(rij[k]/box[k])*box[k];  // <-- pbc's
        rij2 += rij[k]*rij[k];
    }
    // if inside the cutoff radius then calculate ij pair interaction
    if (rij2 <= rcut2)
    {
        rijMag = sqrt(rij2);

        // add soft repulsive force
        radComp = 0.0;
        avgRad = (rad.at(i) + rad.at(j))/2;
        radComp = K*pow(rijMag-avgRad,n);

        if (abs(Ex)>0 || abs(Ey)>0 || abs(Ez)>0)
        {
            // add dipole-dipole interaction
            thetaComp = 0.0;
            EdotRij = rij[0]*Ex+rij[1]*Ey+rij[2]*Ez;
            theta = acos(EdotRij/(rijMag*Emag));
            ScaleFactor = pow(avgRad,6.0)*pow(Emag,2.0)/pow(rijMag,4.0);
            radComp -= ScaleFactor*(3*pow(cos(theta),2.0)-1.0);
            thetaComp = -ScaleFactor*sin(2.0*theta);

            // calculate theta hat direction
            cross(edir,rij,normal);
            cross(normal,rij,thetaHat);
            thetaMag = sqrt(pow(thetaHat[0],2.0)+pow(thetaHat[1],2.0)+pow(thetaHat[2],2.0));
            for (int k=0;k<3;k++) {
                if (thetaMag > 0.00001)
                {
                    thetaHat[k]*=1/thetaMag;
                }
                else
                {
                    thetaHat[k]*=0.0;
                }
            }

            // put all force components together
            for (int k=0;k<3;k++) 
            {
                f[3*i+k] += radComp*rij[k]/rijMag + thetaComp*thetaHat[k];
                f[3*j+k] -= radComp*rij[k]/rijMag + thetaComp*thetaHat[k];
            }
        }
        else
        {
            for (int k=0;k<3;k++) 
            {
                f[3*i+k] += radComp*rij[k]/rijMag;
                f[3*j+k] -= radComp*rij[k]/rijMag;
            }
        }
    }
}


void Dipole::cross(double (&a)[3], double(&b)[3],double (&crs)[3])
{
    crs[0] = a[1]*b[2]-b[1]*a[2];
    crs[1] = a[2]*b[0]-b[2]*a[0];
    crs[2] = a[0]*b[1]-b[0]*a[1];
}