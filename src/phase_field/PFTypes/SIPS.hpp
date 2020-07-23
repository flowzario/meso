/************************************************************
 * DISCLAIMER:
 * SIPS class development continued with SIPSternary
 *
 * This class models ternary phase separation with a binary
 * system of polymer (c) and solvent and the non-solvent/solvent
 * exchange is modeled with a 1D solution to fick's second law.
 * 
 * Focus has been directed to SIPSternary
 *
************************************************************/

# ifndef SIPS_H
# define SIPS_H

# include "../PFBaseClass.hpp"
# include "../PFUtils/SfieldFD.hpp"
# include <mpi.h>


class SIPS : public PFBaseClass {

private:

    const CommonParams& p;
    int current_step;
    int nx,ny,nz;
    int NSdepth;
    int deli,delj,delk;
    int nxyz;
    int outAnalysisInterval;
    int numAnalysisOutputs;
    SfieldFD c;             // polymer concentration (3)
    SfieldFD w;             // non-solvent (water) concentration (1)
    double dt;
    double co;
    double water_CB;
    double phiCutoff;
    double chiPS;
    double chiPN;
    double M;
    double N;
    double kap;
    double A;
    double Tbath;
    double Tinit;
    double noiseStr;
    double gamma;
    double nu;
    double gammaNS;
    double nuNS;
    double D0;
    double DNS;
    double Mweight;
    double Mvolume;
    double mobReSize;
    bool bx,by,bz;
public:

    SIPS(const CommonParams&, const GetPot&);
    ~SIPS();
    void initPhaseField();
    void updatePhaseField();
    void outputPhaseField();
    void setTimeStep(int step) {current_step = step;}
    double secondDerFH(double cc);
    double philliesDiffusion(double cc);
    double getFHcc(double cc);
};

# endif  // SIPS_H
