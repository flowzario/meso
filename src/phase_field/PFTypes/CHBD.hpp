
# ifndef CHBD_H
# define CHBD_H

# include "../PFBaseClass.hpp"
# include "../PFUtils/Sfield.hpp"
# include "../PFUtils/Vfield.hpp"
# include "../PFUtils/ParticlesBDCH.hpp"
# include <complex.h>
# include <fftw3-mpi.h>

class CHBD: public PFBaseClass {

private:

    const CommonParams& p;
    int current_step;
    int nxyz;
    int part_step_skip;
    Sfield c1;
    Sfield c2;
    Sfield cp;
    Vfield k1;
    Sfield k2;
    Sfield k4;
    Sfield * kz; // for applying E-field in z-dir
    ParticlesBDCH particles;
    double co;
    double M;
    double w;
    double kap;
    double eCH;
    fftw_plan p_forward;
    fftw_plan p_backward;
    fftw_complex* dummy;

public:

    CHBD(const CommonParams&, const GetPot&);
    ~CHBD();
    void initPhaseField();
    void updatePhaseField();
    void outputPhaseField();
    void setTimeStep(int step) {current_step = step;}

private:

    void calculateKfields();

};

# endif  // CHBD_H
