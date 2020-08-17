#include "Cell.hh"

Cell::Cell()
{
}
Cell::Cell(const Cell &orig)= default;
Cell::Cell(float ener, float chener,float nuener, float LocalSigma, int lay, int ind_e, int ind_p, int LocalLabel, bool shar, const float (&TopoLayersPix)[6])
{
    TotalEnergy = ener;
    ChargeEnergy = chener;
    NeutralEnergy = nuener;
    sigma = LocalSigma;
    GlobalLayer = lay;
    IndexEta = ind_e;
    IndexPhi = ind_p;
    GlobalLabel = LocalLabel;
    shared = shar;
    FirstWeight = 1;
    SecondWeight = 0;
    FirstLocalMaxLabel = 0;
    SecondLocalMaxLabel = 0;
    ListSize4Sigam = 0;
    for (int l = 0; l < 6; l++ )
    {
        LPixel[l] = TopoLayersPix[l];
    }
}
void Cell::update(float ener, float chener,float nuener, float LocalSigma, int lay, int ind_e, int ind_p, int LocalLabel, bool shar, const int (&TopoLayersPix)[6])
{
    TotalEnergy = ener;
    ChargeEnergy = chener;
    NeutralEnergy = nuener;
    sigma = LocalSigma;
    GlobalLayer = lay;
    IndexEta = ind_e;
    IndexPhi = ind_p;
    GlobalLabel = LocalLabel;
    shared = shar;
    FirstWeight = 1;
    SecondWeight = 0;
    FirstLocalMaxLabel = 0;
    SecondLocalMaxLabel= 0;
    ListSize4Sigam= 0;
    for (int l = 0; l < 6; l++ )
    {
        LPixel[l] = TopoLayersPix[l];
    }

}
void Cell::update(Cell cell)
{
    TotalEnergy = cell.GetTotEnergy();
    ChargeEnergy = cell.GetTotEnergy();
    NeutralEnergy = cell.GetNuEnergy();
    sigma = cell.GetSigma();
    GlobalLabel = cell.GetLabel();
    GlobalLayer = cell.GetLayer();
    IndexEta = cell.GetEta();
    IndexPhi = cell.GetPhi();
    shared = cell.GetShar();
    FirstWeight = cell.GetFirstWeight();
    SecondWeight = cell.GetSecondWeight();
    FirstLocalMaxLabel = cell.GetnewLMaxLabel();
    SecondLocalMaxLabel = cell.GetsecondLMaxLabel();
    cell.SetCellFirstLocmax(CellFirstLocmax);
    cell.SetCellSecondLocmax(CellSecondLocmax);

}



//? Step Cluster maker
void Cell::Add_cellClusterMaker(std::vector<Cell> &ListSeedCells,Cell (&CellArray)[6][kMaxPixel][kMaxPixel], int lay, int eta , int phi)
{
    int LocalListSize4Sigam=ListSize4Sigam;
    Cell &CurrentCell = CellArray[lay][eta][phi];
    float LocalSigma = CurrentCell.GetSigma();
    float LocalLabel = CurrentCell.GetLabel();
    if (LocalSigma>4.5)
    {
        if (LocalLabel > GlobalLabel)
        {
            for (int i = 0; i < 6; i++)
            {
                for (int j = 0; j < kMaxPixel; j++)
                {
                    for (int k = 0; k < kMaxPixel; k++)
                    {
                        Cell &LocalCell= CellArray[i][j][k];
                        int cell_lab = LocalCell.GetLabel();
                        if (cell_lab == LocalLabel)
                            LocalCell.SetLabel(GlobalLabel);
                        else if (cell_lab > LocalLabel)
                            LocalCell.SetLabel(cell_lab-1);
                    }
                }
            }

        }
        else if (LocalLabel!=0 && LocalLabel<GlobalLabel)
        {
            for (int i = 0; i < 6; i++)
            {
                for (int j = 0; j < kMaxPixel; j++)
                {
                    for (int k = 0; k < kMaxPixel; k++)
                    {
                        Cell &LocalCell= CellArray[i][j][k];
                        int cell_lab = LocalCell.GetLabel();
                        if (cell_lab == GlobalLabel)
                            LocalCell.SetLabel(LocalLabel);
                        else if (cell_lab > GlobalLabel)
                            LocalCell.SetLabel(cell_lab-1);
                    }
                }
            }
            GlobalLabel = LocalLabel;
        }

    }
    else if (LocalSigma>2)
    {
        if (LocalLabel > GlobalLabel)
        {
            for (int i = 0; i < 6; i++)
            {
                for (int j = 0; j < kMaxPixel; j++)
                {
                    for (int k = 0; k < kMaxPixel; k++)
                    {
                        Cell &LocalCell= CellArray[i][j][k];
                        int cell_lab = LocalCell.GetLabel();
                        if (cell_lab == LocalLabel)
                        {
                            LocalCell.SetLabel(GlobalLabel);
                        }
                        else if (cell_lab > LocalLabel)
                        {
                            LocalCell.SetLabel(cell_lab-1);
                        }
                    }
                }
            }

        }
        else if (LocalLabel!=0 && LocalLabel<GlobalLabel)
        {
            for (int i = 0; i < 6; i++)
            {
                for (int j = 0; j < kMaxPixel; j++)
                {
                    for (int k = 0; k < kMaxPixel; k++)
                    {
                        Cell &LocalCell= CellArray[i][j][k];
                        int cell_lab = LocalCell.GetLabel();
                        if (cell_lab == GlobalLabel)
                            LocalCell.SetLabel(LocalLabel);
                        else if (cell_lab > GlobalLabel)
                            LocalCell.SetLabel(cell_lab-1);
                    }
                }
            }
            GlobalLabel = LocalLabel;

        }
        else if (LocalLabel==0) //* =0
        {
            CurrentCell.SetLabel(GlobalLabel);
            int list_size =ListSeedCells.size();
            if (ListSeedCells.back().GetSigma()>LocalSigma)
            {
                ListSeedCells.emplace_back(CurrentCell);
            }
            else
            {
                if (list_size==ListSize4Sigam)
                {
                    ListSeedCells.emplace_back(CurrentCell);
                }
                else //* (list_size!=ListSize4Sigam)
                {
                    for (int ii = list_size-1; ii >= LocalListSize4Sigam; ii--)
                    {
                        if (ListSeedCells[ii].GetSigma()>LocalSigma)
                        {
                            ListSeedCells.emplace(ListSeedCells.begin() + ii+1, CurrentCell);
                            break;
                        }
                        if (ii==LocalListSize4Sigam)
                        {
                            ListSeedCells.emplace(ListSeedCells.begin() + ii, CurrentCell);
                        }

                    }
                }



            }

        }



    }
    else //* if (LocalSigma<2)
    {
        if (LocalLabel == 0)  CurrentCell.SetLabel(GlobalLabel);
    }
}
void Cell::NeighborClusterMaker(std::vector<Cell> &ListSeedCells, Cell (&CellArray)[6][kMaxPixel][kMaxPixel], int SizeList)
{
    // const int LPixel[6] = {PixelLayerEM1, PixelLayerEM2, PixelLayerEM3, PixelLayerHCAL1, PixelLayerHCAL2, PixelLayerHCAL3};
    // int up[5] = {(PixelLayerEM1/PixelLayerEM2), (PixelLayerEM2/PixelLayerEM3), (PixelLayerEM3/PixelLayerHCAL1),
    //                     (PixelLayerHCAL1/PixelLayerHCAL2), (PixelLayerHCAL2/PixelLayerHCAL3)};
    float up[5] = {(LPixel[0]/LPixel[1]), (LPixel[1]/LPixel[2]), (LPixel[2]/LPixel[3]),
                        (LPixel[3]/LPixel[4]), (LPixel[4]/LPixel[5])};
    ListSize4Sigam = SizeList;
    for (int i =-1; i < 2 ; i = i + 2)
    {
        int next_layer = GlobalLayer+i;
        if (next_layer>=0 && next_layer<6)
        {

            if (i == 1)
            {
                Add_cellClusterMaker(ListSeedCells, CellArray, next_layer ,(int) floor(((float)IndexEta)/up[next_layer-1]), (int) floor(((float)IndexPhi)/up[next_layer-1]) );
            }
            else//* i==-1
            {
                for (int UpStepEta = 0; UpStepEta < up[next_layer]; UpStepEta++)
                {
                    for (int UpStepPhi = 0; UpStepPhi < up[next_layer]; UpStepPhi++)
                    {
                    Add_cellClusterMaker( ListSeedCells, CellArray, next_layer, IndexEta*up[next_layer]+UpStepEta,   IndexPhi*up[next_layer]+UpStepPhi  );
                    }
                }
            }
        }
    }

    for (int i =-1; i < 2; i++)
    {
        for (int j =-1; j < 2; j++)
        {
            int step_eta = IndexEta+i;
            int step_phi = IndexPhi+j;
            if (step_eta>=0 && step_eta<LPixel[GlobalLayer] && (i!=0 || j!=0))
            {
                Add_cellClusterMaker( ListSeedCells, CellArray, GlobalLayer, step_eta, modulo(step_phi,LPixel[GlobalLayer]));
            }
        }
    }
    ListSeedCells.erase(ListSeedCells.begin());
    for (int it = 0; it < ListSeedCells.size(); it++ )
    {
        ListSeedCells[it].SetLabel(CellArray[ListSeedCells[it].GetLayer()][ListSeedCells[it].GetEta()][ListSeedCells[it].GetPhi()].GetLabel());
    }
}
//? Step Cluster maker end



//? Step Cluster spliter
void Cell::Add_cellLocalMax(std::vector<Cell> &LocalMaxCellsInClust, Cell (&CellArray)[6][kMaxPixel][kMaxPixel], std::vector<Cell> &ShareCellsList, int ClusterLabel, int lay, int eta , int phi)
{

    Cell &CurrentCell = CellArray[lay][eta][phi];
    float ener = CurrentCell.GetTotEnergy();
    float LocalLabel = CurrentCell.GetLabel();
    float newLabel = CurrentCell.GetnewLMaxLabel();
    float secondnewLabel = CurrentCell.GetsecondLMaxLabel();
    if (LocalLabel == ClusterLabel)
    {
        // cout<<"2layer "<<lay<<" eta "<<eta<<" phi "<<phi<<endl;
        if (newLabel==0)
        {
            CurrentCell.SetnewLMaxLabel(FirstLocalMaxLabel);
            CurrentCell.SetCellFirstLocmax(GlobalLayer, IndexEta, IndexPhi);
            // CurrentCell.SetShar(true);
            int list_size =LocalMaxCellsInClust.size();
            if (LocalMaxCellsInClust.back().GetTotEnergy()>ener)
            {
                LocalMaxCellsInClust.emplace_back(CurrentCell);
            }
            else
            {
                if (list_size==ListSize4Sigam)
                {
                    LocalMaxCellsInClust.emplace_back(CurrentCell);

                }
                else
                {
                    for (int ii = list_size-1; ii >= ListSize4Sigam; ii--)
                    {
                        if (LocalMaxCellsInClust[ii].GetTotEnergy()>ener)
                        {
                            LocalMaxCellsInClust.emplace(LocalMaxCellsInClust.begin() + ii+1, CurrentCell);
                            break;
                        }
                        if (ii==ListSize4Sigam)
                        {
                            LocalMaxCellsInClust.emplace(LocalMaxCellsInClust.begin() + ii, CurrentCell);
                        }
                    }
                }
            }
        }
        else // newLAbel!=0
        {
            if (newLabel != FirstLocalMaxLabel)
            {
                if (secondnewLabel == 0)
                {
                    CurrentCell.SetsecondLMaxLabel(FirstLocalMaxLabel);
                    CurrentCell.SetCellSecondLocmax(GlobalLayer, IndexEta, IndexPhi);
                    CurrentCell.SetShar(true);
                    int list_size =LocalMaxCellsInClust.size();
                    for (int ii= 0; ii < list_size; ii++)
                    {
                        if ( (LocalMaxCellsInClust[ii].GetLayer()== CurrentCell.GetLayer()) &&
                             (LocalMaxCellsInClust[ii].GetEta()== CurrentCell.GetEta())     &&
                             (LocalMaxCellsInClust[ii].GetPhi()== CurrentCell.GetPhi()))
                        {
                            LocalMaxCellsInClust.erase(LocalMaxCellsInClust.begin()+ii);
                            break;
                        }
                    }
                    ShareCellsList.push_back(CurrentCell);
                }
                else
                {
                    std::vector<int> indexfirstlocalmax;
                    CurrentCell.GetCellFirstLocmax(indexfirstlocalmax);
                    int layfirstlocmax = indexfirstlocalmax[0];
                    int etafirstlocmax = indexfirstlocalmax[1];
                    int phifirstlocmax = indexfirstlocalmax[2];
                    std::vector<int> indexsecondlocalmax;
                    CurrentCell.GetCellSecondLocmax(indexsecondlocalmax);
                    int laysecondlocmax = indexsecondlocalmax[0];
                    int etasecondlocmax = indexsecondlocalmax[1];
                    int phisecondlocmax = indexsecondlocalmax[2];
                    if (CellArray[layfirstlocmax][etafirstlocmax][phifirstlocmax].GetTotEnergy()<ener)//check TotalEnergy
                    {
                        CurrentCell.SetCellSecondLocmax(layfirstlocmax, etafirstlocmax, phifirstlocmax);
                        CurrentCell.SetsecondLMaxLabel(CurrentCell.GetnewLMaxLabel());
                        CurrentCell.SetnewLMaxLabel(FirstLocalMaxLabel);
                        CurrentCell.SetCellFirstLocmax(GlobalLayer,IndexEta, IndexPhi);
                    }
                    else if (CellArray[laysecondlocmax][etasecondlocmax][phisecondlocmax].GetTotEnergy()<ener)
                    {
                        CurrentCell.SetsecondLMaxLabel(FirstLocalMaxLabel);
                        CurrentCell.SetCellSecondLocmax(GlobalLayer,IndexEta, IndexPhi);
                    }
                }
            }
        }
    }
}
void Cell::NeighborLocalMax(std::vector<Cell> &LocalMaxCellsInClust, Cell (&CellArray)[6][kMaxPixel][kMaxPixel],     std::vector<Cell> &ShareCellsList, int ClusterLabel, int SizeList)
{
    // const int LPixel[6] = {PixelLayerEM1, PixelLayerEM2, PixelLayerEM3, PixelLayerHCAL1, PixelLayerHCAL2, PixelLayerHCAL3};
    // const int up[5] = {(PixelLayerEM1/PixelLayerEM2), (PixelLayerEM2/PixelLayerEM3), (PixelLayerEM3/PixelLayerHCAL1),
    //                      (PixelLayerHCAL1/PixelLayerHCAL2), (PixelLayerHCAL2/PixelLayerHCAL3)};
    float up[5] = {(LPixel[0]/LPixel[1]), (LPixel[1]/LPixel[2]), (LPixel[2]/LPixel[3]),
                        (LPixel[3]/LPixel[4]), (LPixel[4]/LPixel[5])};
    ListSize4Sigam = SizeList;
    for (int i =-1; i < 2 ; i = i + 2)
    {
        int next_layer = GlobalLayer+i;
        if (next_layer>=0 && next_layer<6)
        {

            if (i == 1)
            {
                Add_cellLocalMax(LocalMaxCellsInClust, CellArray, ShareCellsList, ClusterLabel, next_layer ,(int) floor(((float)IndexEta)/up[next_layer-1]), (int) floor(((float)IndexPhi)/up[next_layer-1]));
            }
            else//* i==-1
            {
                for (int UpStepEta = 0; UpStepEta < up[next_layer]; UpStepEta++)
                {
                    for (int UpStepPhi = 0; UpStepPhi < up[next_layer]; UpStepPhi++)
                    {
                    Add_cellLocalMax(LocalMaxCellsInClust, CellArray, ShareCellsList, ClusterLabel, next_layer, (int)(IndexEta*(up[next_layer]))+UpStepEta,   (int)(IndexPhi*(up[next_layer]))+UpStepPhi  );
                    }
                }
            }
        }
    }

    for (int i =-1; i < 2; i++)
    {
        for (int j =-1; j < 2; j++)
        {
            int step_eta = IndexEta+i;
            int step_phi = IndexPhi+j;

            if (step_eta >= 0 && step_eta < LPixel[GlobalLayer] && (i!=0 || j!=0))
            {
                Add_cellLocalMax( LocalMaxCellsInClust, CellArray, ShareCellsList, ClusterLabel, GlobalLayer, step_eta, modulo(step_phi,LPixel[GlobalLayer]));
            }
        }
    }
    LocalMaxCellsInClust.erase(LocalMaxCellsInClust.begin());
    int list_size =  LocalMaxCellsInClust.size();
    for (int it = 0; it < list_size; it++ )
    {
        LocalMaxCellsInClust[it].update(CellArray[LocalMaxCellsInClust[it].GetLayer()][LocalMaxCellsInClust[it].GetEta()][LocalMaxCellsInClust[it].GetPhi()]);
    }
    list_size =  ShareCellsList.size();
    for (int it = 0; it < list_size; it++ )
    {
        ShareCellsList[it].update(CellArray[ShareCellsList[it].GetLayer()][ShareCellsList[it].GetEta()][ShareCellsList[it].GetPhi()]);
    }
}



bool Cell::IsLocalMax(Cell (&CellArray)[6][kMaxPixel][kMaxPixel], int LocalLaeyr, int LocalEta, int LocalPhi)
{
    // const int LPixel[6] = {PixelLayerEM1, PixelLayerEM2, PixelLayerEM3, PixelLayerHCAL1, PixelLayerHCAL2, PixelLayerHCAL3};
    // const int up[5] = {(PixelLayerEM1/PixelLayerEM2), (PixelLayerEM2/PixelLayerEM3), (PixelLayerEM3/PixelLayerHCAL1),
                        //  (PixelLayerHCAL1/PixelLayerHCAL2), (PixelLayerHCAL2/PixelLayerHCAL3)};
    float up[5] = {(LPixel[0]/LPixel[1]), (LPixel[1]/LPixel[2]), (LPixel[2]/LPixel[3]),
                        (LPixel[3]/LPixel[4]), (LPixel[4]/LPixel[5])};
    int num_nei = 0;
    bool IsMax = true;
    float ener = CellArray[LocalLaeyr][LocalEta][LocalPhi].GetTotEnergy();
    if (CellArray[LocalLaeyr][LocalEta][LocalPhi].GetLabel()>0)
    {
        for (int i =-1; i < 2; i = i+2 )
        {
            int next_layer = LocalLaeyr+i;
            if (next_layer>=0 && next_layer<6)
            {
                if (i == 1)
                {
                    if (CellArray[LocalLaeyr+i][(int)floor(((float)LocalEta)/(up[next_layer-1]))][(int)floor(LocalPhi/(up[next_layer-1]))].GetLabel()>0) num_nei++;
                    IsMax = ((CellArray[LocalLaeyr+i][(int)floor(((float)LocalEta)/(up[next_layer-1]))][(int)floor(LocalPhi/(up[next_layer-1]))].GetTotEnergy() < ener) && IsMax);
                }
                else //i == -1
                {
                    for (int UpStepEta = 0; UpStepEta < up[next_layer]; UpStepEta++)
                    {
                        for (int UpStepPhi = 0; UpStepPhi < up[next_layer]; UpStepPhi++)
                        {
                            if (CellArray[LocalLaeyr+i][(int) (LocalEta*(up[next_layer])+UpStepEta)][(int)(LocalPhi*(up[next_layer])+UpStepPhi)].GetLabel()>0) num_nei++;
                            IsMax = ((CellArray[LocalLaeyr+i][(int) (LocalEta*(up[next_layer])+UpStepEta)][(int) (LocalPhi*(up[next_layer])+UpStepPhi)].GetTotEnergy() < ener)&& IsMax);
                        }
                    }
                }
            }
        }
        for (int i =-1; i < 2; i++)
        {
            for (int j =-1; j < 2; j++)
            {
                int step_eta = LocalEta+i;
                int step_phi = LocalPhi+j;
                if (step_eta>=0 && step_eta<LPixel[LocalLaeyr] && (i!=0 || j!=0))
                {
                    if (CellArray[LocalLaeyr][step_eta][modulo(step_phi,LPixel[LocalLaeyr]  )].GetLabel()>0) num_nei++;
                    IsMax = ((CellArray[LocalLaeyr][step_eta][modulo(step_phi,LPixel[LocalLaeyr])].GetTotEnergy() < ener)&& IsMax);
                }
            }
        }
    }
    bool res = false;
    if (num_nei>=4 && IsMax)
    {
        res = true;
    }
    return res;
}

int Cell::modulo( int value, int m)
{
    int mod = value % (int)m;
    if (value < 0) {
        mod += m;
    }
    return mod;
}

float Cell::GetTotEnergy()
{
    return TotalEnergy;
}

float Cell::GetSigma()
{
    return sigma;
}

int Cell::GetLayer()
{
    return GlobalLayer;
}

int Cell::GetEta()
{
    return IndexEta;
}

int Cell::GetPhi()
{
    return IndexPhi;
}

void Cell::SetLabel(int LocalLabel)
{
    GlobalLabel = LocalLabel;
}
int Cell::GetLabel()
{
//    cout<<"GlobalLabel "<<GlobalLabel<<endl;
    return GlobalLabel;
}
void Cell::SetShar(bool shar)
{
    shared = shar;
}
bool Cell::GetShar()
{
    return shared;
}
void Cell::SetFirstWeight(float FWeight)
{
    FirstWeight = FWeight;
}
void Cell::SetSecondWeight(float SWeight)
{
    SecondWeight = SWeight;
}
float Cell::GetFirstWeight()
{
    return FirstWeight;
}
float Cell::GetSecondWeight()
{
    return SecondWeight;
}
void Cell::SetnewLMaxLabel(int LocalLabel)
{
    FirstLocalMaxLabel = LocalLabel;
}
int Cell::GetnewLMaxLabel()
{
    return FirstLocalMaxLabel;
}
void Cell::SetsecondLMaxLabel(int LocalLabel)
{
    SecondLocalMaxLabel = LocalLabel;
}
int Cell::GetsecondLMaxLabel()
{
    return SecondLocalMaxLabel;
}
void Cell::SetCellFirstLocmax(int lay, int eta, int phi)
{
    CellFirstLocmax.clear();
    CellFirstLocmax.push_back(lay);
    CellFirstLocmax.push_back(eta);
    CellFirstLocmax.push_back(phi);
}
void Cell::SetCellFirstLocmax(std::vector<int> &coord)
{
    CellFirstLocmax=coord;
}
void Cell::GetCellFirstLocmax(std::vector<int> &coord)
{
    coord = CellFirstLocmax;
}
void Cell::SetCellSecondLocmax(int lay, int eta, int phi)
{
    CellSecondLocmax.clear();
    CellSecondLocmax.push_back(lay);
    CellSecondLocmax.push_back(eta);
    CellSecondLocmax.push_back(phi);
}
void Cell::SetCellSecondLocmax(std::vector<int> &coord)
{
    CellSecondLocmax=coord;
}
void Cell::GetCellSecondLocmax(std::vector<int> &coord)
{
    coord = CellSecondLocmax;
}
void Cell::SetTotEnergy(float en)
{
    TotalEnergy=en;
}
void Cell::SetChEnergy(float en)
{
    ChargeEnergy=en;
}
float Cell::GetChEnergy()
{
    return ChargeEnergy;
}
void Cell::SetNuEnergy(float en)
{
    NeutralEnergy=en;
}
float Cell::GetNuEnergy()
{
    return NeutralEnergy;
}
void Cell::SetSigma(float LocalSigma)
{
    sigma=LocalSigma;
}