#include <vector>
// #include "MParticle.hh"
using std::vector;

struct MParticle
{
    // int pdg;
    int PosInList;
    float Energy;
};
struct PCell
{
    vector <MParticle> prtcl;
    float chEnergy;
    float nuEnergy;
    float Noise; 
};