// #include <iostream>
#include <vector>
#include "Cell.hh"
#include "GlobalVariables.hh"
using std::vector;

class TopoCluster
{
    public:
        Cell CellArray[6][kMaxPixel][kMaxPixel];
        float NoiseArrayGran[6][kMaxPixel][kMaxPixel];
        vector < vector <float> > FinalListClust;
        // float NoiseValues[6] = {13., 34., 17., 54., 33., 54.};
        float NoiseValues[6] = {13., 34., 41., 75., 50., 25.};
        float LPixel[6];
        TopoCluster( const float (&TotalEnergyArray)[6][kMaxPixel][kMaxPixel], const float (&NoiseArray)[6][kMaxPixel][kMaxPixel], const float (&ChargeEnergyArray)[6][kMaxPixel][kMaxPixel], const float (&NeutralEnergyArray)[6][kMaxPixel][kMaxPixel], const int (&TopoLayersPix)[6]);
        void Run();
        void FillTopo(float (&TopoArray)[6][kMaxPixel][kMaxPixel]) ;
        void FillTotEn(float (&TotalEnergyArray)[6][kMaxPixel][kMaxPixel]);
        void FillChEn(float (&ChargeEnergyArray)[6][kMaxPixel][kMaxPixel]);
        void FillNuEn(float (&NeutralEnergyArray)[6][kMaxPixel][kMaxPixel]);
        void FillNoise(float (&NeutralEnergyArray)[6][kMaxPixel][kMaxPixel]);
};