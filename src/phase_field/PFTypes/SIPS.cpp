/************************************************************
 * DISCLAIMER: 
 * SIPS class development continued with SIPSternary
 *
 * This class models ternary phase separation with a binary
 * system of polymer (c) and solvent and the non-solvent/solvent
 * exchange is modeled with a 1D solution to fick's second law.
 * 
 * Focus has been directed to SIPSternary.
 *
************************************************************/
# include "SIPS.hpp"
# include <iostream>
# include <fstream>



// -------------------------------------------------------------------------
// Constructor:
// -------------------------------------------------------------------------

SIPS::SIPS(const CommonParams& pin,
                   const GetPot& input_params) : p(pin), c(p), w(p)
{

    // ---------------------------------------
    // set needed parameters:
    // ---------------------------------------

    nxyz = p.nx*p.ny*p.nz;
    nx = p.nx;
    ny = p.ny;
    nz = p.nz;
    dt = p.dt;
    bx = input_params("PFApp/bx",0);
    by = input_params("PFApp/by",1);
	bz = input_params("PFApp/bz",1);
    deli = (nz+2)*(ny+2);
	delj = (nz+2);
	delk = 1;
    co = input_params("PFApp/co",0.1);
    water_CB = input_params("PFApp/water_CB",1.0);
    chiPS = input_params("PFApp/chiPS",0.034);
    chiPN = input_params("PFApp/chiPN",2.5);
    phiCutoff = input_params("PFApp/phiCutoff",0.5);
    M = input_params("PFApp/M",1.0);
    kap = input_params("PFApp/kap",1.0);
    N = input_params("PFApp/N",100.0);
    A = input_params("PFApp/A",1.0);
    Tinit = input_params("PFApp/Tinit",273.15);
    noiseStr = input_params("PFApp/noiseStr",0.1);
    D0 = input_params("PFApp/D0", 1.0);
    DNS = input_params("PFApp/DNS",10.0);
    gamma = input_params("PFApp/gamma", 1.0);
    nu = input_params("PFApp/nu",1.0);
    gammaNS = input_params("PFApp/gammaNS",1.0);
    nuNS = input_params("PFApp/nuNS,",1.0);
    mobReSize = input_params("PFApp/mobReSize",1.0);
    Mweight = input_params("PFApp/Mweight",100.0);
    Mvolume = input_params("PFApp/Mvolume",0.1);
    
}



// -------------------------------------------------------------------------
// Destructor:
// -------------------------------------------------------------------------

SIPS::~SIPS()
{

}



// -------------------------------------------------------------------------
// Initialize phase-field method:
// -------------------------------------------------------------------------

void SIPS::initPhaseField()
{

    //	---------------------------------------
    // initialize the concentration field:
    //	---------------------------------------
    srand(time(NULL)*(p.rank+1));   // set the random seed
    for (int i=1; i<nx+1; i++) {
    	  //waterD = i + p.xOff - 1;    // NSdepth counter (water) (1)
        for (int j=1; j<ny+1; j++) {
            for (int k=1; k<nz+1; k++) {
                int ndx = i*deli + j*delj + k*delk;
                double r = (double)rand()/RAND_MAX;
                double val = 0.0;
                val = co + noiseStr*(r-0.5);  
                c.setValue(ndx,val);
                w.setValue(ndx,0.0);
            }
        } 
    } 

    //	---------------------------------------
    // Output the initial configuration:
    //	---------------------------------------

    current_step = 0;
    outputPhaseField();

}



// -------------------------------------------------------------------------
// Step forward in time the phase-field method:
// -------------------------------------------------------------------------

void SIPS::updatePhaseField()
{

    // ---------------------------------------
    // calculate chemical potential & mobility
    // ---------------------------------------

    c.updateBoundaries(bx,by,bz);
    w.updateBoundaries(bx,by,bz);
    MPI::COMM_WORLD.Barrier();
    // ----------------------------------------
    // initializing variables
    // ----------------------------------------
    SfieldFD mu_c(p);
    SfieldFD mob_c(p);
    SfieldFD mu_w(p);
    SfieldFD D_w(p);
    double chi = 0.0;
    for (int i=1; i<nx+1; i++) {
        double T = Tinit; 
 		double kT = T/273.0;
        for (int j=1; j<ny+1; j++) {
            for (int k=1; k<nz+1; k++) {
                int ndx = i*deli + j*delj + k*delk;
                if (i==1) w.setValue(ndx,water_CB);
                double cc = c.getValue(ndx);
                double ww = w.getValue(ndx);
                double cc_fh = 0.0;  
                cc_fh = getFHcc(cc);
                // calculate chi from linear weighted average
                chi = ww*chiPN + (1.0-ww)*chiPS;
                // ---------------------------------------
                // binary first derivative
                // ---------------------------------------
                double df_c = (log(cc_fh) + 1.0)/N - log(1.0-cc_fh) - 1.0 + chi*(1.0-2.0*cc_fh);
				// negative concentrations...
                if (cc <= 0.0) df_c = -1.5*A*sqrt(-cc);
                double lap_c = c.Laplacian(ndx);
                double lap_w = w.Laplacian(ndx);
                mu_c.setValue(ndx,df_c - kap*lap_c);
                mu_w.setValue(ndx,lap_w);
	 	        // --------------------------------------------------------
                // polymer self diffusion (Phillies) and 2nd Derivative FH
                // --------------------------------------------------------
		        double cc_phil = philliesDiffusion(cc);
                double ddf_c = secondDerFH(cc);	
	            ddf_c *= kT;
                double Dp = D0 * exp (- gamma * pow(cc_phil,nu));	
                // -------------------------
                // nonsolvent diffusion
                // --------------------------
                double Dpw = DNS * exp ( - gammaNS * pow(cc_phil,nuNS));
                D_w.setValue(ndx,Dpw);
                // ----------------------
                // mobility
                // ----------------------
                double Mc = Dp/ddf_c;
                if (Mc > 1.0) Mc = 1.0;     // making mobility max = 1
   	            else if (Mc < 0.0) Mc = 1e-12;
	            Mc *= mobReSize;
                if (cc >= phiCutoff) Mc *= 1e-6;
                mob_c.setValue(ndx,Mc);
            }
        }
    }

    // ---------------------------------------
    // update CH equation:
    // ---------------------------------------

    mu_c.updateBoundaries(bx,by,bz);
    mob_c.updateBoundaries(bx,by,bz);
    mu_w.updateBoundaries(bx,by,bz);
    D_w.updateBoundaries(bx,by,bz);

    MPI::COMM_WORLD.Barrier();

    c += p.dt*(mu_c.Laplacian(mob_c)); 
    w += p.dt*DNS*mu_w/*(mu_w.Laplacian(D_w)*/;

 

    // ---------------------------------------
    // Add random fluctuations: set coagulation bath composition
    // ---------------------------------------
    double val = 0;
    for (int i=1; i<nx+1; i++) {
        for (int j=1; j<ny+1; j++) {
            for (int k=1; k<nz+1; k++) {
                int ndx = i*deli + j*delj + k*delk;
                if (i == 1) w.setValue(ndx,water_CB);
                double r = (double)rand()/RAND_MAX;
                double c_check = c.getValue(ndx);
                if (c_check < 0) val =0;
                else if (c_check >= phiCutoff) val = 0;
                else val = noiseStr*(r-0.5);
                // dont want noise added to water bath...
                c.addValue(ndx,p.dt*val);
            }
        }
    }
    
}

// ------------------------------------------------------
// flory huggins -- creating stability
// ------------------------------------------------------
double SIPS::getFHcc(double cc)
{
	double cc_fh = 0.0;
	if (cc < 0.0) { cc_fh = 0.000001; }
	else if (cc > 1.0) { cc_fh = 0.999; }
	else { cc_fh = cc; }
	return cc_fh;
}


// -------------------------------------------------------
// phillies conversion
// -------------------------------------------------------
double SIPS::philliesDiffusion(double cc)
{
	double cc_phil = 0.0;
	if (cc >= 1.0) cc_phil = 1.0 * Mweight/Mvolume; // convert phi to g/L	
	else if (cc < 0.0) cc_phil = 0.000001 * Mweight/Mvolume; // convert phi to g/L 
	else { cc_phil = cc * Mweight/Mvolume; }// convert phi to g/L  
	return cc_phil;
}


// --------------------------------------------------------
// 2nd Derivative FH with phillies concentration
// --------------------------------------------------------
double SIPS::secondDerFH(double cc)
{
	double ddf_c = 0.0;
	if (cc < 0.0) ddf_c = 0.75*A*pow(-cc,-0.5);
	else if (cc == 0.0) ddf_c = 1e3;
	else ddf_c = 0.5* (1.0/(N*cc) + 1.0/(1.0-cc)); 
	return ddf_c;
}


// -------------------------------------------------------------------------
// Write output for the phase-field method:
// -------------------------------------------------------------------------

void SIPS::outputPhaseField()
{
    int iskip = p.iskip;
    int jskip = p.jskip;
    int kskip = p.kskip;
    c.writeVTKFile("c",current_step,iskip,jskip,kskip);
    w.writeVTKFile("w",current_step,iskip,jskip,kskip);
//    n.writeVTKFile("n",current_step,iskip,jskip,kskip);
}
