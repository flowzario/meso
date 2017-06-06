
# include "MesoBase.hpp"

// -------------------------------------------------------------------------
// List of header files that need to be included...
// -------------------------------------------------------------------------

// Cahn-Hilliard classes:
# include "../cahn_hilliard/finite_difference/CHApp.hpp"
# include "../cahn_hilliard/spectral/CHFFTApp.hpp"

// Brownian Dynamics classes:
# include "../brownian_dynamics/BDApp.hpp"

// Lattice-Boltzmann classes:
# include "../lattice_boltzmann/LBApp.hpp"

// Phase-Field classes:
# include "../phase_field/PFApp.hpp"

// FIPI (Fast Interface-Particle Interactions) classes:
//# include "../fipi/FIPIApp.hpp"

// -------------------------------------------------------------------------
// Factory method: this function returns an object determined
// by the string 'specifier':
// {Note: all of the returnable objects inherent from 'MesoBase'}
// -------------------------------------------------------------------------

MesoBase* MesoBase::MesoObjectFactory(string specifier)
{

   // ------------------------------------------------
   // 'GetPot' object containing input parameters:
   // ------------------------------------------------

   GetPot InParams("input.dat");

   // -----------------------------------
   // return the requested object:
   // -----------------------------------

   if (specifier == "CHApp/") return new CHApp(InParams);
   if (specifier == "CHFFTApp/") return new CHFFTApp(InParams);
   if (specifier == "BDApp/") return new BDApp(InParams);
   if (specifier == "LBApp/") return new LBApp(InParams);
   if (specifier == "PFApp/") return new PFApp(InParams);
   //if (specifier == "FIPIApp/") return new FIPIApp(InParams);

}