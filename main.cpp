
#define HDF5

#include <iostream>
#include <cmath>
#include <uni10.hpp>
#include <string>
#include "tclap/CmdLine.h"

#include "timeEvolution2site.h"
#include "helper_compress2.h"
#include "hamiltonians.h"
#include "compress.h"
#include "grow.h"
#include "bigVec.h"


using namespace std;

typedef std::vector<double>     vd;
typedef std::vector<Bond>       vb;
typedef std::vector<UniTensor>  vt;

std::string globalHAMTYPE;
double globalLAM;
double globalDT;

void checkBETA( timeEvolveTwoSite &t_state, std::vector<UniTensor> &ham, int Nsteps, double tstep, int bondDim, std::string output_file, bool do_compress )
{
    if (output_file == "")
        output_file = "output_measurement_beta_dependence.csv";
    csvfile g(output_file,",");
    g.precision(6);
    g << "beta" << "H" << "HH" << "spHeat" << "norm" << "ham_type" << "lam" << "dt_program" << endrow;
    std::vector<UniTensor> mps = t_state.final_mps_state();
    double n  = normMPS(mps);
    double H  = expect(mps,ham);
    double HH = expect2(mps,ham,ham);
    double C  = 0;
    g << 0 << H/n << HH/n << C/n << n << globalHAMTYPE << globalLAM << globalDT << endrow;
    for(int i = 1; i <= Nsteps; ++i){
        t_state.sweep();
        mps = t_state.final_mps_state();
        n  = normMPS(mps);
        H  = expect(mps,ham)/n;
        HH = expect2(mps,ham,ham)/n;
        C  = i*i * (2*tstep)*(2*tstep) * (HH - H*H);
        g << i*(2*tstep) << H << HH << C << n << endrow;
        std::cout << " done sweep and calculation " << std::endl;
        std::cout << " max BondDim size: " << maximum_bond_size(t_state.final_mps_state()) << std::endl;

        if(do_compress){
            t_state.compress_mps_state(bondDim+1,bondDim);
            std::cout << " done compression " << std::endl;
            std::cout << " max BondDim size: " << maximum_bond_size(t_state.final_mps_state()) << std::endl;
        }

        t_state.save(output_file+"_beta="+doublestr(2*i*tstep)+"_STATE");
    }

}

int main(int argc, char** argv)
{
        try{
            TCLAP::CmdLine cmd("Command description message", ' ',  "0.9");
            TCLAP::ValueArg<double>      time_step_Arg ("t","timestep",      "", true,  0.1,          "double");
            TCLAP::ValueArg<double>            lam_Arg ("l","lambda",        "", true,  0.00,         "double");
            TCLAP::ValueArg<int>            length_Arg ("L","length",        "", true,  50,           "int"   );
            TCLAP::ValueArg<double>           kvec_Arg ("k","kvec",          "", false, 25.0/51.0,    "double");
            TCLAP::ValueArg<double>      maxHeight_Arg ("H","maxHeight",     "", true,  2000,         "int"   );
            TCLAP::ValueArg<double>      bandwidth_Arg ("","bandwidth",     "",  true,  100,          "int"   );
            TCLAP::ValueArg<bool>        Qcompress_Arg ("","Qcompress",     "",  true,  true,         "bool"  );
            TCLAP::ValueArg<double>        maxTemp_Arg ("","maxThermal",    "",  true,  10,           "int"   );
            TCLAP::ValueArg<std::string>   hamType_Arg ("","HamiltonianType","", true,  "J1J2aniso",  "string"   );
            TCLAP::ValueArg<double>         aniso1_Arg ("","AdditionalParam1","",true,  -1000,        "int"   );
            TCLAP::ValueArg<double>         aniso2_Arg ("","AdditionalParam2","",false, -1000,        "int"   );
            TCLAP::ValueArg<std::string>   outfile_Arg ("o","outfile",  "", true,  "soms", "string");
            TCLAP::ValueArg<std::string>   storage_Arg ("s","storage",  "", true,  "abc", "string");
            TCLAP::ValueArg<bool>           THorCH_Arg ("","thermalORcheb",  "", true,  true, "bool");
            TCLAP::ValueArg<bool>                D_Arg ("","maxBond",        "", true,  100, "int");
 
            cmd.add( time_step_Arg );
            cmd.add(       lam_Arg );
            cmd.add(    length_Arg );
            cmd.add(      kvec_Arg );
            cmd.add( maxHeight_Arg );
            cmd.add( bandwidth_Arg );
            cmd.add(   outfile_Arg );
            cmd.add(   storage_Arg );
            cmd.add( Qcompress_Arg );
            cmd.add(   maxTemp_Arg );
            cmd.add(   hamType_Arg );
            cmd.add(    THorCH_Arg );
            cmd.add(    aniso1_Arg );
            cmd.add(    aniso2_Arg );
            cmd.add(         D_Arg );
            cmd.parse(argc, argv);
 
        int    maxBond =         D_Arg.getValue();
        double    dt   = time_step_Arg.getValue(); globalDT  = dt;
        double    lam  =       lam_Arg.getValue(); globalLAM = lam;
        int       intL =    length_Arg.getValue(); double L  = intL;
        double    kvec =      kvec_Arg.getValue()*PI;
        int       maxH = maxHeight_Arg.getValue();
        double    band = bandwidth_Arg.getValue();
        bool qCompress = Qcompress_Arg.getValue();
        int       maxT =   maxTemp_Arg.getValue();
   std::string HamType =   hamType_Arg.getValue(); globalHAMTYPE = HamType;
        bool   thORch  =    THorCH_Arg.getValue(); 
        int     aniso1 =   aniso1_Arg.getValue();
        int     aniso2 =   aniso2_Arg.getValue();
 
        std::string outTemp = outfile_Arg.getValue();
        std::string outCheb = storage_Arg.getValue();
 
    /******************************************************************************************************/
    std::unique_ptr<HamiltonianClass> hamiltonian_class;
 
    if(HamType == "XXZ")
        hamiltonian_class = std::make_unique<hamXXZ>(hamXXZ(intL,lam));
    else if (HamType == "J1J2simple")
        hamiltonian_class = std::make_unique<hamJ1J2>(hamJ1J2(intL,lam,1.0));
    else if (HamType == "J1J2aniso")
        hamiltonian_class = std::make_unique<hamJ1J2>(hamJ1J2(intL,lam,aniso1));
    else if (HamType == "XXZext")
        hamiltonian_class = std::make_unique<hamXXZext>(hamXXZext(intL,lam));
    else if (HamType == "XXZspinOne")
        hamiltonian_class = std::make_unique<hamXXZspinOne>(hamXXZspinOne(intL,lam));
    else if (HamType == "XXZsia")
        hamiltonian_class = std::make_unique<hamXXZsia>(hamXXZsia(intL,lam,aniso1));
    else if (HamType == "XXZsiaSpinOne")
        hamiltonian_class = std::make_unique<hamXXZsiaSpinOne>(hamXXZsiaSpinOne(intL,lam,aniso1));
    else
        assert("hamiltonian class not found" == 0);
 
    std::vector<UniTensor> mps_dynamic_thermal_init = hamiltonian_class -> thermalMPS();
    std::vector<UniTensor> mpo_dynamic_thermal_heis = hamiltonian_class -> thermalMPO();
 
    assert(mps_dynamic_thermal_init.size() == 2*intL);
    assert(mps_dynamic_thermal_init.size() == mpo_dynamic_thermal_heis.size());
 
    /******************************************************************************************************/
 
    timeEvolveTwoSite t_state( // 2nd defn: timeEvolveTwoSite
                mps_dynamic_thermal_init.size(), // total length
                                             dt, // time step 
                       mps_dynamic_thermal_init,
                       mpo_dynamic_thermal_heis
            );
 
//   std::cout << "dt:              " <<  dt                       << std::endl;
//   std::cout << "lam:             " << lam                       << std::endl;
//   std::cout << "physical length: " << mps_dynamic_thermal_init.size()/2 << std::endl;
//
//   double Exp_Beta = 2*(maxT)*dt;
//   std::vector<UniTensor> tmp;
//   if (!thORch){
//       int N = mps_dynamic_thermal_init.size();
//       tmp.resize(N);
//       for(int i = 0; i < N; ++i){
//           tmp[i] = UniTensor(outTemp+"_beta="+doublestr(Exp_Beta)+"_STATE"+properfilename(0,0,i));
//           std::cout << "  - - -  Loaded tensor: " << tmp[i].elemNum() << " elements" << std::endl;
//       }
//       std::cout << "NORM:" << std::endl;
//       std::cout << normMPS(tmp) << std::endl;
//   }
//
//   double    Wprime = 1.000 - 0.025/2.0;
//
//   if(thORch){
//       checkBETA( t_state, mpo_dynamic_thermal_heis, maxT, dt, maxBond, outTemp, qCompress );
//   } else {
//       // cheb-evolution combine and test with thermal-state later
//       //      thermal MPS, maxBondDim, minBondDim, epsilon, littleLen,  
//       //      Delta,           beta,  kvec,   N,      W, Wprime, sweeps, storage_location, firstTime, saveOrNot
//       stem cheb_beta0 = stem( t_state.final_mps_state(), hamiltonian_class, maxBond+1, maxBond, 1e-5, 
//               intL, lam, Exp_Beta, kvec, maxH, band, Wprime, 30, outCheb, true, false);
//   }
//
//
    } catch (TCLAP::ArgException &e)  // catch any exceptions
    { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }



    return 0;
}

