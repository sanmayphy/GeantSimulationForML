#include "GlobalVariables.hh"
#include "CaloRConstants.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4NistManager.hh"
#include "cmath"
// using std::vector;
// using std::sort;
using namespace std;
inline bool sortcolf( const vector<float>& v1, const vector<float>& v2 ) 
{ 
    return v1[2] > v2[2]; 
} 
class PFlow
{
    public:
        PFlow(vector < vector <double> > &TrajList, vector < vector <float> > &TopoList, const float (&TotalEnergyArray)[6][kMaxPixel][kMaxPixel], const float (&LabelArray)[6][kMaxPixel][kMaxPixel]);
        void Fill(float (&PFlowA)[6][kMaxPixel][kMaxPixel]);
    private:
        void TrackMatch (vector <double>  &TrajList, vector < vector <float> > &TopoList);
        void RecoveringShowerSplit( vector <double>  &TrajList, vector < vector <float> > &TopoList);
        void LHED(vector <double>  &TrajList, const float (&TotalEnergyArray)[6][kMaxPixel][kMaxPixel], const float (&LabelArray)[6][kMaxPixel][kMaxPixel]) ;
        void CellByCellSubs( vector <double>  &TrajList, const float (&TotalEnergyArray)[6][kMaxPixel][kMaxPixel], const float (&LabelArray)[6][kMaxPixel][kMaxPixel]);
        float CellVoluem(int layer, int IndexEta, int IndexPhi);
        float PFlowArray[6][kMaxPixel][kMaxPixel];
        float EnergyDensityArray[6][kMaxPixel][kMaxPixel];
        vector <float> RateOfIncrease;
        float ErefTopref = 0.7;//0.625;//0.54
        float sigmaEref =  0.25;//2686;
        float SumEdep = 0;
};
