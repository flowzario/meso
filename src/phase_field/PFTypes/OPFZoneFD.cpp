
# include "OPFZoneFD.hpp"
# include <iostream>
# include <fstream>


// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

OPFZoneFD::OPFZoneFD(const CommonParams& pin,
                   const GetPot& input_params) : p(pin), c(p)
{

    // ---------------------------------------
    // set needed parameters:
    // ---------------------------------------

    nxyz = p.nx*p.ny*p.nz;
    nx = p.nx;
    ny = p.ny;
    nz = p.nz;
    deli = (nz+2)*(ny+2);
	delj = (nz+2);
	delk = 1;
    co = input_params("PFApp/co",0.5);
    M = input_params("PFApp/M",1.0);
    noiseStr = input_params("PFApp/noiseStr",0.1);
    wzone = input_params("PFApp/wzone",10.0);
	vzone = input_params("PFApp/vzone",1.0);
	S = input_params("PFApp/S",0.05);
	chiN = input_params("PFApp/chiN",15);
	chiN2 = chiN*chiN;
	chiN3 = chiN2*chiN;
	chiN4 = chiN2*chiN2;

	c_2 = S*(-0.00368*chiN2 - 1.964*chiN +15.99); 
	c_4 = S*(0.01882*chiN2 + 3.176*chiN -26.9); 
	c_5 = S*(-0.000003452*chiN4 + 0.0004019*chiN3 - 0.01791*chiN2 + 0.3735*chiN -1.684); 
	c_7 = S*(0.000004339*chiN4 - 0.000446*chiN3 + 0.01729*chiN2 - 0.32*chiN + 11.59); 
	singleWellZero = sqrt(-c_2/(2*c_4))+0.5;
}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

OPFZoneFD::~OPFZoneFD()
{

}



// -------------------------------------------------------------------------
// Initialize phase-field method:
// -------------------------------------------------------------------------

void OPFZoneFD::initPhaseField()
{

    //	---------------------------------------
    // initialize the concentration field:
    //	---------------------------------------

    srand(time(NULL)*(p.rank+1));   // set the random seed
    for (int i=1; i<nx+1; i++) {
        for (int j=1; j<ny+1; j++) {
            for (int k=1; k<nz+1; k++) {
                int ndx = i*deli + j*delj + k*delk;
                double r = (double)rand()/RAND_MAX;
                double val = co + 0.1*(r-0.5);
                c.setValue(ndx,val);
            }
        }
    }

}



// -------------------------------------------------------------------------
// Step forward in time the phase-field method:
// -------------------------------------------------------------------------

void OPFZoneFD::updatePhaseField()
{

    // ---------------------------------------
    // calculate chemical potential & mobility
    // ---------------------------------------

    c.updatePBCFluxY();
    MPI::COMM_WORLD.Barrier();
	

    SfieldFD mu(p);
    SfieldFD mob(p);
    for (int i=1; i<nx+1; i++) {
		
		//calculate mobility
		double Mc;
		if (vzone != 0) //if zone annealing is active
		{
			//calculate time 
			double time = current_step*p.dt;
				
			//calculate offset
			int xf = int(time*vzone - 0.5*wzone - p.xOff);

			//calculate local mobility
			Mc = 0.5*(1 - tanh(6.0*double(i-xf)/wzone));
		}
		else 
		{ 
			//assign constant mobility
			Mc = M;
		}
        
		for (int j=1; j<ny+1; j++) {
            for (int k=1; k<nz+1; k++) {
                
				//get current position
				int ndx = i*deli + j*delj + k*delk;
				//get current concentration 
				double cc = c.getValue(ndx);
				//calculate chemical potential
                double dFdc = 2*c_2*(cc-0.5) + 4*c_4*(cc-0.5)*(cc-0.5)*(cc-0.5); 
				
				//write mobility values
				mob.setValue(ndx, Mc);
				//write potential values
				mu.setValue(ndx,dFdc);

            }//k
        }//j
    }//i
	
    // ---------------------------------------
    // update CH equation:
    // ---------------------------------------

	//update boundary conditions for mu and mobility
    
	MPI::COMM_WORLD.Barrier();
	mu.updatePBCFluxY();
    mob.updatePBCFluxY();

	//calculate first laplacian of CH right hand side
	SfieldFD inside = mu - c_5*c.Laplacian();
	inside.updatePBCFluxY();
	
	//calculate total right hand side of CH equation
	SfieldFD RHS = inside.Laplacian(mob); 	


	for (int i=1; i<nx+1; i++) {
        for (int j=1; j<ny+1; j++) {
            for (int k=1; k<nz+1; k++) {
                int ndx = i*deli + j*delj + k*delk;
				
				// Apply energy barrier for Di-block system:
				RHS.addValue(ndx,-c_7*(c.getValue(ndx)-co));    //comment to simulate binary fluid system
				
				// Add random fluctuations:
				double r = (double)rand()/RAND_MAX;
                double val = noiseStr*(r-0.5);
                RHS.addValue(ndx,val);
            }
        }
    }
	
	RHS.updatePBCFluxY();
	c += p.dt*RHS;
}



// -------------------------------------------------------------------------
// Write output for the phase-field method:
// -------------------------------------------------------------------------

void OPFZoneFD::outputPhaseField()
{
    int iskip = p.iskip;
    int jskip = p.jskip;
    int kskip = p.kskip;
    c.writeVTKFile("c",current_step,iskip,jskip,kskip);
}