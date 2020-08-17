// #include "Cell.hh"
#include "TopoCluster.hh"
//* Constructor creates arrays that imitate diffrent granularity
TopoCluster::TopoCluster(const float (&TotalEnergyArray)[6][kMaxPixel][kMaxPixel], const float (&NoiseArray)[6][kMaxPixel][kMaxPixel], const float (&ChargeEnergyArray)[6][kMaxPixel][kMaxPixel], const float (&NeutralEnergyArray)[6][kMaxPixel][kMaxPixel])
{
    for (int l = 0; l < 6; l++)
    {
        SumPix(l, TotalEnergyArray[l], NoiseArray[l], ChargeEnergyArray[l], NeutralEnergyArray[l]);
    }
}

void TopoCluster::SumPix(int layer,  const float (&TotalEnergyArray)[kMaxPixel][kMaxPixel], const float (&NoiseArray)[kMaxPixel][kMaxPixel],const float (&ChargeEnergyArray)[kMaxPixel][kMaxPixel],const float (&NeutralEnergyArray)[kMaxPixel][kMaxPixel])
{
    int size = LayersPix[layer];
    //* If we have some energy dep. from particle (False)
    //* If event is empty and only noise was generated (True).
    //* Do not do unnecessary calculations
    bool OnlyNoiseCell = true;
    //* Fill NoiseArrayGran and CellArray
    for (int e = 0; e < kMaxPixel; e++)
    {
        for (int p = 0; p < kMaxPixel; p++)
        {
            NoiseArrayGran[layer][e][p]= 0;
            OnlyNoiseCell = true && (TotalEnergyArray[e][p]==NoiseArray[e][p]);
            CellArray[layer][e][p].update(0, 0, 0, 0, layer, e, p, 0, false); //Todo Could be done when cell is initialized
        }
    }
    int scale = kMaxPixel/size;
    //* Fill CellArray with in case of diffrent granularity
    for (int it_x = 0; it_x < size; it_x++)
    {
        for (int it_y = 0; it_y < size; it_y++)
        {
            float CellTotSignal    = NoiseArray[it_x][it_y];
            float CellChSignal = 0;//CellTotSignal;
            float CellNuSignal = 0;//CellTotSignal;
            NoiseArrayGran[layer][it_x][it_y] = CellTotSignal;
            // if ( !OnlyNoiseCell)
            // {
                for (int ii = scale*it_x; ii < scale*(it_x+1); ii++ )
                {
                    for (int jj = scale*it_y; jj < scale*(it_y+1); jj++ )
                    {
                        CellTotSignal = CellTotSignal + TotalEnergyArray[ii][jj];
                        CellChSignal = CellChSignal + ChargeEnergyArray[ii][jj];
                        CellNuSignal = CellNuSignal + NeutralEnergyArray[ii][jj];
                    }
                }
            // }
            CellArray[layer][it_x][it_y].SetSigma(CellTotSignal/NoiseValues[layer]);
            CellArray[layer][it_x][it_y].SetTotEnergy(CellTotSignal);
            CellArray[layer][it_x][it_y].SetNuEnergy(CellNuSignal);
            CellArray[layer][it_x][it_y].SetChEnergy(CellChSignal);
        }
    }
}

void TopoCluster::FillTopo(float (&TopoArray)[6][kMaxPixel][kMaxPixel])
{
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < kMaxPixel; j++)
        {
            for (int k = 0; k < kMaxPixel; k++)
            {
                TopoArray[i][j][k] = CellArray[i][j][k].GetLabel();
            }
        }
    }

}
void TopoCluster::FillTotEn(float (&TotalEnergyArray)[6][kMaxPixel][kMaxPixel])
{
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < kMaxPixel; j++)
        {
            for (int k = 0; k < kMaxPixel; k++)
            {
                if (CellArray[i][j][k].GetLabel()!=0) TotalEnergyArray[i][j][k] = CellArray[i][j][k].GetTotEnergy();
                else TotalEnergyArray[i][j][k] = 0;
            }
        }
    }

}
void TopoCluster::FillChEn(float (&ChargeEnergyArray)[6][kMaxPixel][kMaxPixel])
{
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < kMaxPixel; j++)
        {
            for (int k = 0; k < kMaxPixel; k++)
            {
                if (CellArray[i][j][k].GetLabel()!=0) ChargeEnergyArray[i][j][k] = CellArray[i][j][k].GetChEnergy();
                else ChargeEnergyArray[i][j][k] = 0;
            }
        }
    }

}
void TopoCluster::FillNuEn(float (&NeutralEnergyArray)[6][kMaxPixel][kMaxPixel])
{
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < kMaxPixel; j++)
        {
            for (int k = 0; k < kMaxPixel; k++)
            {
                if (CellArray[i][j][k].GetLabel()!=0) NeutralEnergyArray[i][j][k] = CellArray[i][j][k].GetNuEnergy();
                else NeutralEnergyArray[i][j][k] = 0;
            }
        }
    }

}
void TopoCluster::FillNoise(float (&Noise)[6][kMaxPixel][kMaxPixel])
{
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < kMaxPixel; j++)
        {
            for (int k = 0; k < kMaxPixel; k++)
            {
                if (CellArray[i][j][k].GetLabel()!=0) Noise[i][j][k] = NoiseArrayGran[i][j][k];
                else Noise[i][j][k] = 0;
            }
        }
    }

}
void TopoCluster::Run()
{
    // * ListSeedCells - list of 4 sigma Cell
    vector< Cell > ListSeedCells;
    //* There was problem of relocation of vectors
    ListSeedCells.reserve(5000);
    //* find all Seed cells with above 4 sigma signal
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < kMaxPixel; j++)
        {
            for (int k = 0; k < kMaxPixel; k++)
            {

                Cell &RefToCell = CellArray[i][j][k];
                float sigma = RefToCell.GetSigma();

                if (sigma > 4 )
                {
                    int SizeListSeedCells =ListSeedCells.size();
                    if (SizeListSeedCells==0)
                    {
                        RefToCell.SetLabel(1);
                        ListSeedCells.push_back(RefToCell);
                    }
                    else//* if (SizeListSeedCells!=0)
                    {

                        Cell &LastCellInList = ListSeedCells.back();

                        if (LastCellInList.GetSigma()>sigma)
                        {
                            RefToCell.SetLabel(LastCellInList.GetLabel()+1);
                            ListSeedCells.push_back(RefToCell);
                        }
                        else//* if (LastCellInList.GetSigma()<=sigma)
                        {
                            for (int ii = SizeListSeedCells-1; ii >= 0 ; ii--)
                            {
                                if (ListSeedCells[ii].GetSigma()>sigma)
                                {
                                    int label = ListSeedCells[ii].GetLabel();
                                    RefToCell.SetLabel(label+1);
                                    ListSeedCells.insert(ListSeedCells.begin() + ii+1, RefToCell);
                                    break;
                                }
                                else//* (ListSeedCells[ii].GetSigma()<=sigma)
                                {
                                    int label = ListSeedCells[ii].GetLabel();
                                    ListSeedCells[ii].SetLabel(label+1);
                                }
                                if (ii == 0)
                                {
                                    RefToCell.SetLabel(1);
                                    ListSeedCells.insert(ListSeedCells.begin(),RefToCell);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    // cout<<"Point 2 "<<endl;

    //* collecting Cells
    int SizeListSeedCells = ListSeedCells.size();
    for (int i = 0; i < SizeListSeedCells;i++)
    {
        CellArray[ListSeedCells[i].GetLayer()][ListSeedCells[i].GetEta()][ListSeedCells[i].GetPhi()].SetLabel(ListSeedCells[i].GetLabel());
    }

    //* Cluster Maker Step
    while (ListSeedCells.size()>0)
    {
        ListSeedCells.front().NeighborClusterMaker(ListSeedCells,CellArray,SizeListSeedCells);
        SizeListSeedCells--;
        if (SizeListSeedCells==0)
        {
            SizeListSeedCells=ListSeedCells.size();
        }
    }
    //* Get number of laybels
    int NumberOfLabels = 0;
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < kMaxPixel; j++)
        {
            for (int k = 0; k < kMaxPixel; k++)
            {
                int label = CellArray[i][j][k].GetLabel();
                if (label > NumberOfLabels) NumberOfLabels = label;
            }
        }
    }
    // cout<<"Point 4 "<<endl;


    // * ListLocalMaxCells - list of Cells with energy bigger then 500 MeV
    vector< Cell> ListLocalMaxCells;
    //* There was problem of relocation of vectors
    ListLocalMaxCells.reserve(5000);
    //* Find all Local Max. cells
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < kMaxPixel; j++)
        {
            for (int k = 0; k < kMaxPixel; k++)
            {
                Cell &RefToCell = CellArray[i][j][k];
                float CellEner = RefToCell.GetTotEnergy();
                //* 4 sigma
                if (CellEner > 500 && (i >= 0 && i < 3))
                {
                    if (RefToCell.IsLocalMax(CellArray,i,j,k))
                    {
                        int SizeListLocalMaxCells =ListLocalMaxCells.size();
                        if (SizeListLocalMaxCells==0)
                        {
                            ListLocalMaxCells.push_back(RefToCell);
                        }
                        else // if (SizeListLocalMaxCells!=0)
                        {

                            Cell &LastCellInList = ListLocalMaxCells.back();

                            if (LastCellInList.GetTotEnergy()>CellEner)
                            {
                                ListLocalMaxCells.push_back(RefToCell);
                            }
                            else//if (LastCellInList.GetTotEnergy()<=CellEner)
                            {
                                for (int ii = SizeListLocalMaxCells-1; ii >= 0 ; ii--)
                                {
                                    if (ListLocalMaxCells[ii].GetTotEnergy()>CellEner)
                                    {
                                        ListLocalMaxCells.insert(ListLocalMaxCells.begin() + ii+1, RefToCell);
                                        break;
                                    }
                                    if (ii == 0)
                                    {
                                        ListLocalMaxCells.insert(ListLocalMaxCells.begin(),RefToCell);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }


    // * ShareCellsList - list of Share cells between two local max.
    vector< Cell> ShareCellsList;
    //* There was problem of relocation of vectors
    ShareCellsList.reserve(5000);
    //* NewClustCount -- cluster number (instead of i), since we split cluster total number of clusters change
    int NewClustCount = 0;
    //* SplitClustCount -- amount of new clusters after one iteration of spliting
    int SplitClustCount = 0;

    for (int ClustCount = 1; ClustCount <= NumberOfLabels; ClustCount++)
    {
        NewClustCount = NewClustCount + SplitClustCount + 1;
        int SplitClustCount = 0;
        vector< Cell> LocalMaxCellsInClust;
        LocalMaxCellsInClust.reserve(5000);

        int SizeListLocalMaxCells = ListLocalMaxCells.size();
        int LabelLocalMax = 0;
        for (int CellinListLocalMaxCells = 0; CellinListLocalMaxCells < SizeListLocalMaxCells; CellinListLocalMaxCells++)
        {
            if(ListLocalMaxCells[CellinListLocalMaxCells].GetLabel()==NewClustCount)
            {
                ListLocalMaxCells[CellinListLocalMaxCells].SetnewLMaxLabel(NewClustCount+LabelLocalMax);
                ListLocalMaxCells[CellinListLocalMaxCells].SetCellFirstLocmax(ListLocalMaxCells[CellinListLocalMaxCells].GetLayer(),ListLocalMaxCells[CellinListLocalMaxCells].GetEta(),ListLocalMaxCells[CellinListLocalMaxCells].GetPhi());
                LocalMaxCellsInClust.push_back(ListLocalMaxCells[CellinListLocalMaxCells]);
                CellArray[ListLocalMaxCells[CellinListLocalMaxCells].GetLayer()][ListLocalMaxCells[CellinListLocalMaxCells].GetEta()][ListLocalMaxCells[CellinListLocalMaxCells].GetPhi()].SetnewLMaxLabel(NewClustCount+LabelLocalMax);
                CellArray[ListLocalMaxCells[CellinListLocalMaxCells].GetLayer()][ListLocalMaxCells[CellinListLocalMaxCells].GetEta()][ListLocalMaxCells[CellinListLocalMaxCells].GetPhi()].SetCellFirstLocmax(ListLocalMaxCells[CellinListLocalMaxCells].GetLayer(),ListLocalMaxCells[CellinListLocalMaxCells].GetEta(),ListLocalMaxCells[CellinListLocalMaxCells].GetPhi());
                LabelLocalMax++;
            }
        }

        if (LocalMaxCellsInClust.size()>1)
        {
            SplitClustCount = LabelLocalMax-1;

            int SizeListSeedCells = LocalMaxCellsInClust.size();
            while (LocalMaxCellsInClust.size()>0)
            {
                LocalMaxCellsInClust.front().NeighborLocalMax(LocalMaxCellsInClust, CellArray, ShareCellsList, NewClustCount, SizeListSeedCells);
                SizeListSeedCells--;
                if (SizeListSeedCells==0)
                {
                    SizeListSeedCells=LocalMaxCellsInClust.size();
                }
            }
            for (int CellinListLocalMaxCells = 0; CellinListLocalMaxCells < SizeListLocalMaxCells;CellinListLocalMaxCells++)
            {
                if(ListLocalMaxCells[CellinListLocalMaxCells].GetLabel()>NewClustCount)
                {
                    ListLocalMaxCells[CellinListLocalMaxCells].SetLabel(ListLocalMaxCells[CellinListLocalMaxCells].GetLabel()+LabelLocalMax-1);
                }
            }
            for (int i = 0; i < 6; i++)
            {
                for (int j = 0; j < kMaxPixel; j++)
                {
                    for (int k = 0; k < kMaxPixel; k++)
                    {
                        int label = CellArray[i][j][k].GetLabel();
                        if (label > NewClustCount)
                        {
                            CellArray[i][j][k].SetLabel(label + LabelLocalMax-1);
                        }
                        else if (NewClustCount==label)
                        {
                                CellArray[i][j][k].SetLabel(CellArray[i][j][k].GetnewLMaxLabel());

                        }
                    }
                }
            }
        }
    }


    // int NumberOfLabels = 0;
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < kMaxPixel; j++)
        {
            for (int k = 0; k < kMaxPixel; k++)
            {
                if (CellArray[i][j][k].GetLabel() > NumberOfLabels )
                {
                    NumberOfLabels = CellArray[i][j][k].GetLabel();
                }
            }
        }
    }
    int ShareCellsListsize = ShareCellsList.size();
    for (int CellinShareCellsList = 0; CellinShareCellsList < ShareCellsListsize; CellinShareCellsList++)
    {
        CellArray[ShareCellsList[CellinShareCellsList].GetLayer()][ShareCellsList[CellinShareCellsList].GetEta()][ShareCellsList[CellinShareCellsList].GetPhi()].SetLabel(0);
        ShareCellsList[CellinShareCellsList].SetLabel(0);
    }
    //* TopoClusters list [Energy, eta, phi, z, Energy_ch, Energy_nu, NoiseArray]
    vector < vector <float> > ListClust (NumberOfLabels, vector <float> (9,0));//before (4,0);
    //* We save abs energy for COM calculation
    vector <float> EnerList(NumberOfLabels,0);

    //* List for calculation of COM when we have periodic boundary condition
    vector < vector <float> > PhiList(NumberOfLabels, vector <float> (2,0));
    int count  = 0;
    float dEta[6] = {d_eta*(kMaxPixel/PixelLayerEM1),   d_eta*(kMaxPixel/PixelLayerEM2),   d_eta*(kMaxPixel/PixelLayerEM3),
                    d_eta*(kMaxPixel/PixelLayerHCAL1), d_eta*(kMaxPixel/PixelLayerHCAL2), d_eta*(kMaxPixel/PixelLayerHCAL3)};
    float dPhi[6] = {divided_tube_dPhi*(kMaxPixel/PixelLayerEM1),   divided_tube_dPhi*(kMaxPixel/PixelLayerEM2),   divided_tube_dPhi*(kMaxPixel/PixelLayerEM3),
                    divided_tube_dPhi*(kMaxPixel/PixelLayerHCAL1), divided_tube_dPhi*(kMaxPixel/PixelLayerHCAL2), divided_tube_dPhi*(kMaxPixel/PixelLayerHCAL3)};
    float dDeep[6] = {(r_inn_ECAL1+(r_out_ECAL1-r_inn_ECAL1)/2)/50.0, (r_inn_ECAL2+(r_out_ECAL2-r_inn_ECAL2)/2)/50.0, (r_inn_ECAL3+(r_out_ECAL3-r_inn_ECAL3)/2)/50.0, 
                    (r_inn_HCAL1+(r_out_HCAL1-r_inn_HCAL1)/2)/50.0, (r_inn_HCAL2+(r_out_HCAL2-r_inn_HCAL2)/2)/50.0, (r_inn_HCAL3+(r_out_HCAL3-r_inn_HCAL3)/2)/50.0};
    //Fill TopoCluster List without share cells
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < kMaxPixel; j++)
        {
            for (int k = 0; k < kMaxPixel; k++)
            {
                int label = CellArray[i][j][k].GetLabel()-1;
                if (label >= 0 && CellArray[i][j][k].GetShar()==false)
                {
                    float CellEnergy =  CellArray[i][j][k].GetTotEnergy();

                    ListClust[label][0] += CellEnergy;
                    ListClust[label][6] += CellArray[i][j][k].GetChEnergy();
                    ListClust[label][7] += CellArray[i][j][k].GetNuEnergy();
                    ListClust[label][8] += NoiseArrayGran[i][j][k];
                    EnerList[label] += abs(CellEnergy);
                    float ListEta = ListClust[label][1];
                    float ListDeep = ListClust[label][3];
                    float CosPhi = PhiList[label][0];
                    float SinPhi = PhiList[label][1];
                    float hdeep = dDeep[i];//ListClust[label][3];
                    float eta = -1 * eta_max + dEta[i] * (j+0.5);
                    float phi = dPhi[i] * (k+0.5);

                    float TopoTotEn = EnerList[label];
                    if (TopoTotEn!=0)
                    {
                        ListClust[label][1]+=abs(CellEnergy)*(eta-ListEta)/TopoTotEn;
                        ListClust[label][3]+=abs(CellEnergy)*(hdeep-ListDeep)/TopoTotEn;
                        PhiList[label][0] += abs(CellEnergy)*(cos(phi)-CosPhi)/TopoTotEn;
                        PhiList[label][1] += abs(CellEnergy)*(sin(phi)-SinPhi)/TopoTotEn;
                    }
                }
            }
        }
    }

    ShareCellsListsize = ShareCellsList.size();
    //* Add to TopoClust share cells
    for (int kk = 0; kk<ShareCellsListsize; kk++)
    {
        // G4cout<<"Orig"<<ShareCellsList[kk].GetLayer()<<" "<<ShareCellsList[kk].GetEta()<<" "<<ShareCellsList[kk].GetPhi()<<G4endl;
        float ShareCellEnergy = ShareCellsList[kk].GetTotEnergy();
        //* FirstLocalMax (SecondLocalMax) about TopoClusters associated with this local max
        int LabelOfFirstLocalMax = ShareCellsList[kk].GetnewLMaxLabel();
        int LabelOfSecondLocalMax= ShareCellsList[kk].GetnewLMaxLabel();
        float EnergyOfFirstLocalMax = ListClust[LabelOfFirstLocalMax-1][0];
        float EnergyOfSecondLocalMax = ListClust[LabelOfSecondLocalMax-1][0];
        int i = ShareCellsList[kk].GetLayer();
        int j = ShareCellsList[kk].GetEta();
        int k = ShareCellsList[kk].GetPhi();
        float eta = -1 * eta_max + dEta[i] * (j+0.5);
        float phi = dPhi[i] * (k+0.5);

        float  DistancePhiOfFirstLocalMax = abs(phi-ListClust[LabelOfFirstLocalMax-1][2]);
        if (DistancePhiOfFirstLocalMax > M_PI)
        {
            DistancePhiOfFirstLocalMax = abs(DistancePhiOfFirstLocalMax - 2*M_PI);
        }
        // float DistanceOfFirstLocalMax  = pow(
        //                                         pow(ShareCellsList[kk].GetEta()-ListClust[LabelOfFirstLocalMax-1][1],2) +
        //                                         pow(DistancePhiOfFirstLocalMax, 2) +
        //                                         pow(ShareCellsList[kk].GetLayer()-ListClust[LabelOfFirstLocalMax-1][3],2),
        //                                     0.5);
        float thetaCell = 2 * atan(exp(-1*eta));
        float thetaClus = 2 * atan(exp(-1*ListClust[LabelOfFirstLocalMax-1][1])); 
        float sinCell = sin(thetaCell);
        float sinClus = sin(thetaClus);
        float HCell = dDeep[ShareCellsList[kk].GetLayer()];
        float HClus = ListClust[LabelOfFirstLocalMax-1][3];
        float RCell = HCell/sinCell;
        float RClus = HClus/sinClus; 
        float DistanceOfFirstLocalMax  = pow( pow(RCell,2) + pow(RClus,2) - 2*RCell*RClus*
                                            (sinCell * sinClus * cos(DistancePhiOfFirstLocalMax) + cos(thetaCell)*cos(thetaClus)), 0.5);
        float  DistancePhiOfSecondLocalMax = abs(phi-ListClust[LabelOfFirstLocalMax-1][2]);
        if (DistancePhiOfSecondLocalMax > M_PI)
        {
            DistancePhiOfSecondLocalMax = abs(DistancePhiOfSecondLocalMax - 2*M_PI);
        }
        // float DistanceOfSecondLocalMax  = pow(
        //                                         pow(ShareCellsList[kk].GetEta()-ListClust[LabelOfSecondLocalMax-1][1],2) +
        //                                         pow(DistancePhiOfSecondLocalMax, 2) +
        //                                         pow(ShareCellsList[kk].GetLayer()-ListClust[LabelOfSecondLocalMax-1][3],2),
        //                                     0.5);
        thetaClus = 2 * atan(exp(-1*ListClust[LabelOfSecondLocalMax-1][1])); 
        sinClus = sin(thetaClus);
        HClus = ListClust[LabelOfSecondLocalMax-1][3];
        RClus = HClus/sinClus; 
        float DistanceOfSecondLocalMax  = pow( pow(RCell,2) + pow(RClus,2) - 2*RCell*RClus*
                                            (sinCell * sinClus * cos(DistancePhiOfSecondLocalMax) + cos(thetaCell)*cos(thetaClus)), 0.5);

        float r = exp(DistanceOfFirstLocalMax - DistanceOfSecondLocalMax);
        float WeightForFirstLocalMax = EnergyOfFirstLocalMax/(EnergyOfFirstLocalMax + r*EnergyOfSecondLocalMax);
        float WeightForSecondLocalMax = 1 - WeightForFirstLocalMax;
        ListClust[LabelOfFirstLocalMax -1][0] += WeightForFirstLocalMax * ShareCellEnergy;
        ListClust[LabelOfSecondLocalMax-1][0] += WeightForSecondLocalMax * ShareCellEnergy;



        float IndexEtaOfFirstLocalMax  = ListClust[LabelOfFirstLocalMax  - 1][1];
        float IndexEtaOfSecondLocalMax = ListClust[LabelOfSecondLocalMax - 1][1];
        float CosPhiOfFirstLocalMax  = PhiList[LabelOfFirstLocalMax  - 1][0];
        float SinPhiOfFirstLocalMax  = PhiList[LabelOfFirstLocalMax  - 1][1];
        float CosPhiOfSecondLocalMax = PhiList[LabelOfSecondLocalMax - 1][0];
        float SinPhiOfSecondLocalMax = PhiList[LabelOfSecondLocalMax - 1][1];
        float ZofFirstLocalMax = ListClust[LabelOfFirstLocalMax - 1][3];
        float ZofSecondLocalMax = ListClust[LabelOfSecondLocalMax - 1][3];
        EnerList[LabelOfFirstLocalMax - 1] += abs(WeightForFirstLocalMax * ShareCellEnergy); 
        EnerList[LabelOfSecondLocalMax - 1]+= abs(WeightForSecondLocalMax * ShareCellEnergy);
        float TotEnergyOfFirstLocalMax = EnerList[LabelOfFirstLocalMax - 1];
        float TotEnergyOfSecondLocalMax = EnerList[LabelOfSecondLocalMax - 1];

        if (TotEnergyOfFirstLocalMax!=0)
        {
            ListClust[LabelOfFirstLocalMax -1][1]+=abs(WeightForFirstLocalMax * ShareCellEnergy)*(eta-IndexEtaOfFirstLocalMax)/TotEnergyOfFirstLocalMax;
            PhiList[LabelOfFirstLocalMax -1][0] += abs(WeightForFirstLocalMax * ShareCellEnergy)*(cos(phi)-CosPhiOfFirstLocalMax)/TotEnergyOfFirstLocalMax;
            PhiList[LabelOfFirstLocalMax -1][1] += abs(WeightForFirstLocalMax * ShareCellEnergy)*(sin(phi)-SinPhiOfFirstLocalMax)/TotEnergyOfFirstLocalMax;
            ListClust[LabelOfFirstLocalMax -1][3]+=abs(WeightForFirstLocalMax * ShareCellEnergy)*(HCell-ZofFirstLocalMax)/TotEnergyOfFirstLocalMax;
        }
        if (TotEnergyOfSecondLocalMax!=0)
        {
            ListClust[LabelOfSecondLocalMax-1][1]+=abs(WeightForSecondLocalMax * ShareCellEnergy)*(eta-IndexEtaOfSecondLocalMax)/TotEnergyOfSecondLocalMax;
            PhiList[LabelOfSecondLocalMax-1][0] += abs(WeightForSecondLocalMax * ShareCellEnergy)*(cos(phi)-CosPhiOfSecondLocalMax)/TotEnergyOfSecondLocalMax;
            PhiList[LabelOfSecondLocalMax-1][1] += abs(WeightForSecondLocalMax * ShareCellEnergy)*(sin(phi)-SinPhiOfSecondLocalMax)/TotEnergyOfSecondLocalMax;
            ListClust[LabelOfSecondLocalMax-1][3] += abs(WeightForSecondLocalMax * ShareCellEnergy)*(HCell-ZofSecondLocalMax)/TotEnergyOfSecondLocalMax;
        }

        CellArray[ShareCellsList[kk].GetLayer()][ShareCellsList[kk].GetEta()][ShareCellsList[kk].GetPhi()].SetLabel(LabelOfFirstLocalMax);
        CellArray[ShareCellsList[kk].GetLayer()][ShareCellsList[kk].GetEta()][ShareCellsList[kk].GetPhi()].SetFirstWeight(WeightForFirstLocalMax);
        CellArray[ShareCellsList[kk].GetLayer()][ShareCellsList[kk].GetEta()][ShareCellsList[kk].GetPhi()].SetSecondWeight(WeightForSecondLocalMax);
        ListClust[LabelOfFirstLocalMax -1][6] += WeightForFirstLocalMax * CellArray[i][j][k].GetChEnergy();
        ListClust[LabelOfSecondLocalMax-1][6] += WeightForSecondLocalMax * CellArray[i][j][k].GetChEnergy();
        ListClust[LabelOfFirstLocalMax -1][7] += WeightForFirstLocalMax * CellArray[i][j][k].GetNuEnergy();
        ListClust[LabelOfSecondLocalMax-1][7] += WeightForSecondLocalMax * CellArray[i][j][k].GetNuEnergy();
        ListClust[LabelOfFirstLocalMax -1][8] += WeightForFirstLocalMax * NoiseArrayGran[i][j][k];
        ListClust[LabelOfSecondLocalMax-1][8] += WeightForSecondLocalMax * NoiseArrayGran[i][j][k];
    }

    //* Find COM in phi axis
    // G4cout<<"A0"<<G4endl;
    for (int index = 0; index < NumberOfLabels; index++)
    {
        ListClust[index][2] = atan2f(-1*PhiList[index][1], -1*PhiList[index][0]) + M_PI;
    }
    // G4cout<<"A1"<<G4endl;
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < kMaxPixel; j++)
        {
            for (int k = 0; k < kMaxPixel; k++)
            {
                int label = CellArray[i][j][k].GetLabel();
                // G4cout<<"A2"<<G4endl;
                if (label > 0)
                {
                    // G4cout<<"A3"<<G4endl;
                    float eta = -1 * eta_max + dEta[i] * (j+0.5);
                    float phi = dPhi[i] * (k+0.5);
                    float CellEnergy =  CellArray[i][j][k].GetTotEnergy();
                    // G4cout<<"A4"<<G4endl;
                    if (CellArray[i][j][k].GetShar())
                    {
                        // G4cout<<i<<" "<<j<<" "<<k<<G4endl;
                        // G4cout<<"A5"<<G4endl;
                        int FirstLabel = CellArray[i][j][k].GetnewLMaxLabel();
                        float etaClus = ListClust[FirstLabel-1][1];
                        float argumEta = pow(eta-etaClus,2);
                        float FirstWeight =  CellArray[i][j][k].GetFirstWeight();
                        ListClust[FirstLabel-1][4] += FirstWeight*abs(CellEnergy)*argumEta;
                        float phiClus = ListClust[FirstLabel-1][2];
                        float DistancePhi = abs(phi-phiClus);
                        if (DistancePhi > M_PI)
                        {
                            DistancePhi = abs(DistancePhi - 2*M_PI);
                        }
                        float argumPhi = pow(DistancePhi,2);
                        ListClust[FirstLabel-1][5] += FirstWeight*abs(CellEnergy)*argumPhi;
                        // G4cout<<"A6"<<G4endl;
                        int SecondLabel = CellArray[i][j][k].GetsecondLMaxLabel();
                        // G4cout<<"A6.0"<<G4endl;
                        // G4cout<<SecondLabel-1<<G4endl;
                        etaClus = ListClust[SecondLabel-1][1];
                        // G4cout<<"A6.1"<<G4endl;
                        argumEta = pow(eta-etaClus,2);
                        // G4cout<<"A6.2"<<G4endl;
                        float SecondWeight =  CellArray[i][j][k].GetFirstWeight();
                        // G4cout<<"A6.3"<<G4endl;
                        ListClust[SecondLabel-1][4] += FirstWeight*abs(CellEnergy)*argumEta;
                        // G4cout<<"A6.4"<<G4endl;
                        phiClus = ListClust[SecondLabel-1][2];
                        // G4cout<<"A6.5"<<G4endl;
                        DistancePhi = abs(phi-phiClus);
                        // G4cout<<"A6.6"<<G4endl;
                        if (DistancePhi > M_PI)
                        {
                            DistancePhi = abs(DistancePhi - 2*M_PI);
                        }
                        // G4cout<<"A6.7"<<G4endl;
                        argumPhi = pow(DistancePhi,2);
                        // G4cout<<"A6.8"<<G4endl;
                        ListClust[SecondLabel-1][5] += SecondWeight*abs(CellEnergy)*argumPhi;
                        // G4cout<<"A7"<<G4endl;
                    }
                    else
                    {
                        // G4cout<<"A8"<<G4endl;
                        float etaClus = ListClust[label-1][1];
                        float argumEta = pow(eta-etaClus,2);
                        // G4cout<<"A9"<<G4endl;
                        ListClust[label-1][4] += abs(CellEnergy)*argumEta;
                        float phiClus = ListClust[label-1][2];
                        float  DistancePhi = abs(phi-phiClus);
                        // G4cout<<"A10"<<G4endl;
                        if (DistancePhi > M_PI)
                        {
                            DistancePhi = abs(DistancePhi - 2*M_PI);
                        }
                        // G4cout<<"A11"<<G4endl;
                        float argumPhi = pow(DistancePhi,2);
                        ListClust[label-1][5] += abs(CellEnergy)*argumPhi;
                        // G4cout<<"A12"<<G4endl;
                    }
                }
            }
        }
    }

    for (int topoclust = 0; topoclust < ListClust.size(); topoclust++)
    {
        if (ListClust[topoclust][0]>0)
        {
            ListClust[topoclust][4]/=EnerList[topoclust];
            ListClust[topoclust][5]/=EnerList[topoclust];
            ListClust[topoclust][4]=sqrt(ListClust[topoclust][4]);
            ListClust[topoclust][5]=sqrt(ListClust[topoclust][5]);
            if (ListClust[topoclust][5]<0.05) ListClust[topoclust][5]=0.05;
            if (ListClust[topoclust][4]<0.05) ListClust[topoclust][4]=0.05;
        }
        else
        {
            ListClust.erase(ListClust.begin()+(topoclust));
            EnerList.erase(EnerList.begin()+(topoclust));
            for (int i = 0; i < 6; i++)
            {
                for (int j = 0; j < kMaxPixel; j++)
                {
                    for (int k = 0; k < kMaxPixel; k++)
                    {
                        int label = CellArray[i][j][k].GetLabel();
                        if (label==(topoclust+1))
                        {
                            CellArray[i][j][k].SetLabel(0);
                        }

                    }
                }
            }
            for (int i = 0; i < 6; i++)
            {
                for (int j = 0; j < kMaxPixel; j++)
                {
                    for (int k = 0; k < kMaxPixel; k++)
                    {
                        int label = CellArray[i][j][k].GetLabel();
                        if (label>(topoclust+1))
                        {
                            CellArray[i][j][k].SetLabel(label-1);
                        }

                    }
                }
            }
            topoclust--;
        }
    } 

    FinalListClust = ListClust;

}