// #include <iostream>
#include "CaloRConstants.hh"
#include<cmath>
#include <algorithm>
#include <vector>


class Cell
{
    private:
        float TotalEnergy;
        float ChargeEnergy;
        float NeutralEnergy;
        float sigma;
        int GlobalLayer;
        int IndexEta;
        int IndexPhi;
        int GlobalLabel;
        bool shared;
        float FirstWeight;
        float SecondWeight;
        int FirstLocalMaxLabel;
        int SecondLocalMaxLabel;
        int ListSize4Sigam;
        std::vector<int> CellFirstLocmax;
        std::vector<int> CellSecondLocmax;
        float LPixel[6];
    public:
        Cell();
        Cell(const Cell& orig);
        Cell(float ener, float chener,float nuener, float LocalSigma, int lay, int ind_e, int ind_p, int LocalLabel, bool shar, const float (&TopoLayersPix)[6]);
        // ~Cell();
        void NeighborClusterMaker (std::vector<Cell> &ListSeedCells, Cell (&CellArray)[6][kMaxPixel][kMaxPixel], int SizeList);
        void Add_cellClusterMaker (std::vector<Cell> &ListSeedCells,Cell (&CellArray)[6][kMaxPixel][kMaxPixel], int lay, int eta , int phi);
        void NeighborLocalMax(std::vector<Cell> &LocalMaxCellsInClust, Cell (&CellArray)[6][kMaxPixel][kMaxPixel],     std::vector<Cell> &ShareCellsList, int ClusterLabel, int SizeList);
        void Add_cellLocalMax(std::vector<Cell> &LocalMaxCellsInClust, Cell (&CellArray)[6][kMaxPixel][kMaxPixel], std::vector<Cell> &ShareCellsList, int ClusterLabel, int lay, int eta , int phi);
        int modulo( int value, int m);
        bool IsLocalMax(Cell (&CellArray)[6][kMaxPixel][kMaxPixel], int LocalLaeyr, int LocalEta, int LocalPhi);
        void SetLabel(int LocalLabel);
        int GetLabel();
        void SetShar(bool shar);
        bool GetShar();
        float GetTotEnergy();
        void SetTotEnergy(float en);
        float GetChEnergy();
        void SetChEnergy(float en);
        float GetNuEnergy();
        void SetNuEnergy(float en);
        float GetSigma();
        void SetSigma(float LocalSigma);
        int GetLayer();
        int GetEta();
        int GetPhi();
        void update(float ener, float chener,float nuener, float LocalSigma, int lay, int ind_e, int ind_p, int LocalLabel, bool shar, const int (&TopoLayersPix)[6]);
        void update(Cell cell);
        float GetEta_min();
        float GetEta_max();
        float GetPhi_min();
        float GetPhi_max();
        void SetnewLMaxLabel(int LocalLabel);
        int GetnewLMaxLabel();
        void SetsecondLMaxLabel(int LocalLabel);
        int GetsecondLMaxLabel();
        void SetCellFirstLocmax(int lay, int eta, int phi);
        void SetCellFirstLocmax(std::vector<int> &coord);
        void GetCellFirstLocmax(std::vector<int> &coord);
        void SetCellSecondLocmax(int lay, int eta, int phi);
        void SetCellSecondLocmax(std::vector<int> &coord);
        void GetCellSecondLocmax(std::vector<int> &coord);
        void SetFirstWeight(float FWeight);
        void SetSecondWeight(float SWeight);
        float GetFirstWeight();
        float GetSecondWeight();
};