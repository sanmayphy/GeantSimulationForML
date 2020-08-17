#include "H02PFlow.hh"

PFlow::PFlow(vector < vector <double> > &TrajList, vector < vector <float> > &TopoList, const float (&TotalEnergyArray)[6][kMaxPixel][kMaxPixel], const float (&LabelArray)[6][kMaxPixel][kMaxPixel])
{
    // sort(TrajList.begin(), TrajList.end(),sortcol);
    // sort(TrajList.begin(), TrajList.end(),sortcol);
    // G4cout<<"C0"<<G4endl;
    // RecoveringShowerSplit(TrajList,TopoList);
    // G4cout<<"C1"<<G4endl;
    
    int numberOfTrekcs = TrajList.size();
    for (int i = 0; i < 6; i++)
    {
        
        for (int j = 0; j < LayersPix[0]; j++)
        {
            for (int k = 0; k < LayersPix[0]; k++)
            {
                PFlowArray[i][j][k] =  TotalEnergyArray[i][j][k];
            }
        }
    }
    G4cout<<"C2"<<G4endl;
    for (int track = 0; track < numberOfTrekcs; track++)
    {
        TrackMatch(TrajList[track],TopoList);
        RecoveringShowerSplit(TrajList[track],TopoList);
        G4cout<<"C2.1"<<G4endl;
        LHED(TrajList[track], TotalEnergyArray, LabelArray);
        G4cout<<"C2.2"<<G4endl;
        CellByCellSubs(TrajList[track], TotalEnergyArray, LabelArray);
    }
    // G4cout<<"C3"<<G4endl; 

}
void PFlow::Fill(float (&PFlowA)[6][kMaxPixel][kMaxPixel])
{
    for (int i = 0; i < 6; i++)
    {
        
        for (int j = 0; j < LayersPix[0]; j++)
        {
            for (int k = 0; k < LayersPix[0]; k++)
            {
                PFlowA[i][j][k] =  PFlowArray[i][j][k];
            }
        }
    }

}

void PFlow::TrackMatch(vector <double>  &TrajList, vector < vector <float> > &TopoList) 
{
    int numberOfTopocLust = TopoList.size();
    float trackMom = TrajList[2];
      // G4cout<<"A9.1"<<G4endl;
    for (int clust = 0;clust < numberOfTopocLust; clust++ )
    {
        if ((TopoList[clust][0]/TrajList[2])>0.1)
        {
            float deta = TrajList[8] - TopoList[clust][1];
            float dphi = abs(TrajList[9] - TopoList[clust][2]);
            if (dphi > M_PI)
            {
                dphi -= 2*M_PI;
            }
            float R = sqrt(sqr(dphi/TopoList[clust][5])+sqr(deta/TopoList[clust][4]));
            if (R < TrajList[28])
            {
            TrajList[29] = clust+1;
            TrajList[28] = R;            
            }
        }
    }
      // G4cout<<"A9.8"<<G4endl;
      // if (TrajList[28]>2.0)
      // {
      //   runAction->ChargeParticlTraj.erase(runAction->ChargeParticlTraj.begin()+track);
      //   track--;
      // }
      // G4cout<<"A9.9"<<G4endl;
}
void PFlow::RecoveringShowerSplit( vector <double>  &TrajList, vector < vector <float> > &TopoList)
{
    // G4cout<<"F0"<<G4endl; 
    int numberOfTopocLust = TopoList.size();
    // int numberOfTrekcs = TrajList.size();
    float threshold = -0.35;
    // G4cout<<"F1"<<G4endl; 
    // for (int track = 0; track < numberOfTrekcs; track++)
    // {
    G4cout<<"F2"<<G4endl; 
    float mom = TrajList[2];
    // G4cout<<"TrajList[track][2] "<<TrajList[track][2]<<G4endl;
    // G4cout<<"((int) TrajList[track][29])-1 "<<((int) TrajList[track][29])-1<<G4endl;
    float en_topo = TopoList[((int) TrajList[29])-1][0];
    
    float S = (en_topo - mom*ErefTopref)/(sigmaEref*mom);
    G4cout<<"F3"<<G4endl; 
    if (S<threshold)
    {
        G4cout<<"!mom*ErefTopref! "<< mom*ErefTopref <<G4endl; 
        for (int clust = 0; clust < numberOfTopocLust; clust++)
        {
            G4cout<<"F4"<<G4endl; 
            float deta = TrajList[4] - TopoList[clust][1];
            float dphi = abs(TrajList[5] - TopoList[clust][2]);
            G4cout<<"F5"<<G4endl; 
            if (dphi > M_PI)
            {
                dphi -= 2*M_PI;
            }
            float R = sqrt(sqr(dphi)+sqr(deta));
            G4cout<<"R "<<R<<G4endl;
            if ((((int) TrajList[29])-1!=clust) && (R<0.4) )
            {
                TrajList.push_back(clust+1);
                G4cout<<"clust "<<clust<<G4endl; 
            }
            
        }
    }
    G4cout<<"F6"<<G4endl; 
    // }

}
void PFlow::LHED( vector <double>  &Traj, const float (&TotalEnergyArray)[6][kMaxPixel][kMaxPixel], const float (&LabelArray)[6][kMaxPixel][kMaxPixel])
{
    float constant = 0.035;//0.00236672;
    int track_size = Traj.size();
    float EnergyDensity[7] = {0,0,0,0,0,0,0};
    long double Lambda_int_ECAL = 37.578;
    long double Deepth[7] ={0, (r_inn_ECAL1-r_out_ECAL1)/Lambda_int_ECAL, (r_inn_ECAL2-r_out_ECAL2)/Lambda_int_ECAL, (r_inn_ECAL3-r_out_ECAL3)/Lambda_int_ECAL, 1.5, 4.1,1.8 }; 
    // float RateOfIncrease[6] = {0,0,0,0,0,0};
    // vector <float> RateOfIncrease (6,0);
    for (int clustmatch = 29; clustmatch < track_size; clustmatch++)
    {
        for (int i = 0; i < 6; i++)
        {
            
            for (int j = 0; j < LayersPix[i]; j++)
            {
                for (int k = 0; k < LayersPix[i]; k++)
                {
                    if (LabelArray[i][j][k]==Traj[clustmatch])
                    {
                        float dEta = d_eta*(kMaxPixel/LayersPix[i]);
                        float dPhi = divided_tube_dPhi*(kMaxPixel/LayersPix[i]);
                        float Eta = j*dEta-eta_max;
                        float Phi = k*dPhi;
                        float cellvol = CellVoluem(i,j,k);
                        float deta = Traj[4+i*4] - Eta;
                        float dphi = abs(Traj[5+i*4] - Phi);// Which value Phi
                        if (dphi > M_PI)
                        {
                            dphi -= 2*M_PI;
                        }
                        float R = (sqr(dphi)+sqr(deta));
                        float LHED_var= TotalEnergyArray[i][j][k]/cellvol;
                        float Weight_var = exp(-1*(R)/(constant*sqr(LayersPix[i]))) ;
                        EnergyDensityArray[i][j][k] = LHED_var*Weight_var;
                        EnergyDensity[i+1]+=EnergyDensityArray[i][j][k];
                        // G4cout<<G4endl;
                        // // G4cout<<"R "<<R<<G4endl;
                        // // G4cout<<"LHED_var "<<LHED_var<<G4endl;
                        // // G4cout<<"Weight_var "<<Weight_var<<G4endl;
                        
                        // G4cout<<G4endl;
                    }
                }
            }
        }

    }
    for (int i = 0; i < 6; i++)
    {
        RateOfIncrease.push_back((EnergyDensity[i+1]-EnergyDensity[i])/(Deepth[i+1]-Deepth[i]));

    }

}
float PFlow::CellVoluem(int layer, int IndexEta, int IndexPhi)
{
    float X0_CAL[6] = {3.897,3.897,3.897,2.357,2.357,2.357};
    float R_inn[6] = {r_inn_ECAL1,r_inn_ECAL2,r_inn_ECAL3,r_inn_HCAL1,r_inn_HCAL2,r_inn_HCAL3};
    float R_out[6] = {r_out_ECAL1,r_out_ECAL2,r_out_ECAL3,r_out_HCAL1,r_out_HCAL2,r_out_HCAL3};
    float dEta = d_eta*(kMaxPixel/LayersPix[layer]);
    
    float Eta = IndexEta*dEta-eta_max;
    float voluem = 0;
    if (Eta>=0)
    {
        voluem = ((2.0/3.0)*M_PI*abs(sinh(Eta+dEta)-sinh(Eta))*(pow(R_out[layer]/X0_CAL[layer],3)-pow(R_inn[layer]/X0_CAL[layer],3)))/LayersPix[layer];
    }
    else
    {
        voluem = ((2.0/3.0)*M_PI*abs(sinh(Eta-dEta)-sinh(Eta))*(pow(R_out[layer]/X0_CAL[layer],3)-pow(R_inn[layer]/X0_CAL[layer],3)))/LayersPix[layer];
    }

    return voluem;
}

void PFlow::CellByCellSubs(vector <double> &Traj, const float (&TotalEnergyArray)[6][kMaxPixel][kMaxPixel], const float (&LabelArray)[6][kMaxPixel][kMaxPixel])
{
    G4cout<<"D0"<<G4endl;


    float pflowneutral = 0;
    for (int i = 0; i < 6; i++)
    {
        
        for (int j = 0; j < LayersPix[i]; j++)
        {
            for (int k = 0; k < LayersPix[i]; k++)
            {
                pflowneutral += PFlowArray[i][j][k];
            }
        }
    }
    G4cout<<"pflowneutralInPF "<<pflowneutral<<G4endl;
    float  max_ = 0;
    int layer_max = 0;
    for (int i = 0; i<6; i++)
    {
        if (RateOfIncrease[i] > max_)
        {
            max_ = RateOfIncrease[i];
            layer_max = i; 
        }
    }
    G4cout<<"D1"<<G4endl;
    float trk_pos_Eta = Traj[4+layer_max*4]; 
    float trk_pos_Phi = Traj[5+layer_max*4]; 
    float mom = Traj[2];
    int track_size = Traj.size();
    float E_clust_tot = 0;
    vector < vector <float> > CellList;
    G4cout<<"D2"<<G4endl;
    float Rest_energy = 0;
    for (int i = 0; i < 6; i++)
    {
        float dEta = d_eta*(kMaxPixel/LayersPix[i]);
        float dPhi = divided_tube_dPhi*(kMaxPixel/LayersPix[i]);
        
        for (int j = 0; j < LayersPix[i]; j++)
        {
            for (int k = 0; k < LayersPix[i]; k++)
            {
                if (LabelArray[i][j][k]!=0)
                {
                    for (int clustmatch = 29; clustmatch < track_size; clustmatch++)
                    {
                        if (LabelArray[i][j][k]==Traj[clustmatch])
                        {
                            float cell_Eta = j*dEta-eta_max;
                            float cell_Phi = k*dPhi;
                            float deta = trk_pos_Eta - cell_Eta;
                            float dphi = abs(trk_pos_Phi - cell_Phi);
                            if (dphi > M_PI)
                            {
                                dphi -= 2*M_PI;
                            }
                            E_clust_tot+=PFlowArray[i][j][k];
                            float R = sqrt(sqr(dphi)+sqr(deta));
                            vector <float> Cell;
                            Cell.push_back(i);
                            Cell.push_back(j);
                            Cell.push_back(k);
                            Cell.push_back(R);
                            Cell.push_back(EnergyDensityArray[i][j][k]);
                            Cell.push_back(PFlowArray[i][j][k]);
                            CellList.push_back(Cell);
                            Rest_energy +=PFlowArray[i][j][k];
                        }
                    }
                    
                }
            }
        }
    }
    
    float E_dep = Traj[2]*ErefTopref;
    SumEdep+=E_dep;
    G4cout<<"E_dep "<< E_dep <<G4endl;
    G4cout<<"E_clust_tot "<< E_clust_tot <<G4endl;
    G4cout<<"Rest_energy "<< Rest_energy <<G4endl;
    if (E_dep>=E_clust_tot)
    {
        G4cout<<"No rings needed! "<<G4endl;
        Rest_energy = 0;
        
        for (int i = 0; i < 6; i++)
        {
            
            for (int j = 0; j < LayersPix[i]; j++)
            {
                for (int k = 0; k < LayersPix[i]; k++)
                {
                    if (LabelArray[i][j][k]!=0)
                    {
                        for (int clustmatch = 29; clustmatch < track_size; clustmatch++)
                        {
                            if (LabelArray[i][j][k]==Traj[clustmatch])  PFlowArray[i][j][k]=0;
                        }
                    }
                }
            }
        }
    }
    else
    {
        int numRing = 1;
        vector < vector <float> > Rings;
        vector <float> Ring;
        float sizeCellList= CellList.size();
        for (int i = 0; i < 6; i++)
        {
            float dEta = d_eta*(kMaxPixel/LayersPix[i]);
            float factor = 1.2*dEta;
            // G4cout<<"LayersPix[i]/2+1 "<<LayersPix[i]/2+1<<G4endl;
            for (int r = 0; r < LayersPix[i]/2+1; r++)
            { 
                // G4cout<<"sizeCellList "<<sizeCellList<<G4endl;
                Ring.clear();
                // vector <float> Ring;
                Ring.push_back(i);
                Ring.push_back(0);
                Ring.push_back(0);
                Ring.push_back(numRing);
                for (int cell = 0; cell < sizeCellList; cell++)
                {
                    int i_cell = CellList[cell][0];
                    // G4cout<<"sqr(CellList[cell][3]) "<<sqr(CellList[cell][3])<<G4endl;
                    // G4cout<<"sqr((r+1)*(factor)) "<<sqr((r+1)*(factor))<<G4endl;
                    // G4cout<<"sqr((r)*(factor)) "<<sqr((r)*(factor))<<G4endl;
                    // G4cout<<"sqr(CellList[cell][5]) "<<sqr(CellList[cell][5])<<G4endl;
                    // G4cout<<"layer "<<(int) CellList[cell][0]<<" eta "<<(int) CellList[cell][1]<<" phi "<<(int) CellList[cell][2]<<G4endl;
                    // G4cout<<G4endl;
                    if ((i_cell==i) && ((sqr((r+1)*(factor)) >= sqr(CellList[cell][3])) && (sqr(r*factor) <= sqr(CellList[cell][3]))))
                    {
                        // G4cout<<"Attention!!!"<<G4endl;
                        Ring[2]+=CellList[cell][4];
                        Ring[1]+=CellList[cell][5];
                        if (CellList[cell].size()==6)
                        {
                            CellList[cell].push_back(numRing);
                        }
                        else
                        {
                            G4cout<<"Attention!!Something wrong!!"<<G4endl;
                        }
                        
                    }
                        
                }
                numRing++;
                Rings.push_back(Ring);
            }
        }
        // Rings.erase(Rings.begin());
        G4cout<<"D6"<<G4endl;
        int sizeRings = Rings.size();
        sort(Rings.begin(), Rings.end(),sortcolf);
        // for (int ring = 0; ring < sizeRings; ring++)
        // {
        //     G4cout<<G4endl;
        //     G4cout<<"Ring "<<ring<<G4endl;
        //     G4cout<<"RingN "<<Rings[ring][3]<<G4endl;
        //     G4cout<<"EnergyDensityArray "<<  Rings[ring][2]<<G4endl;
        //     G4cout<<"TotalEnergyArray "  <<  Rings[ring][1]<<G4endl;
        //     G4cout<<G4endl;
        // }
        
        for (int ring = 0; ring < sizeRings; ring++)
        {
            
            if (E_dep <= Rings[ring][1])
            {
                G4cout<<"Substruct not all ring "<<E_dep<<" "<<Rings[ring][1]<<G4endl;
                float NuFrac = 1-E_dep/Rings[ring][1];
                E_dep = 0 ;
                for (int cell = 0; cell < sizeCellList; cell++)
                {
                    if (CellList[cell][6]==Rings[ring][3])
                    {
                        PFlowArray[(int) CellList[cell][0]][(int) CellList[cell][1]][(int) CellList[cell][2]]=NuFrac*PFlowArray[(int) CellList[cell][0]][(int) CellList[cell][1]][(int) CellList[cell][2]];
                        Rest_energy = Rest_energy - ( (1 - NuFrac)*CellList[cell][5]);
                        CellList[cell][5] = NuFrac*CellList[cell][5];
                        
                        CellList[cell][4] = NuFrac*CellList[cell][4]; 
                    }
                }
                break;
            }
            else
            {
                G4cout<<"Substruct all ring "<<E_dep<<" "<<Rings[ring][1]<<G4endl;
                E_dep = E_dep - Rings[ring][1];
                for (int cell = 0; cell < sizeCellList; cell++)
                {
                    if (CellList[cell][6]==Rings[ring][3])
                    {
                        
                        Rest_energy = Rest_energy - CellList[cell][5];
                        PFlowArray[(int) CellList[cell][0]][(int) CellList[cell][1]][(int) CellList[cell][2]]=0;
                        CellList[cell][5] = 0.0;
                        CellList[cell][4] = 0.0;

                    }
                }
            }
            pflowneutral = 0;
            for (int i = 0; i < 6; i++)
            {
                for (int j = 0; j < LayersPix[i]; j++)
                {
                    for (int k = 0; k < LayersPix[i]; k++)
                    {
                        pflowneutral += PFlowArray[i][j][k];
                    }
                }
            }
            G4cout<<"pflowneutralInRing "<<pflowneutral<<G4endl;
            
        } 
        if (Rest_energy<=1.5*sigmaEref*mom)
        {
            G4cout<<"All energy substructed"<<G4endl;
            for (int cell = 0; cell < sizeCellList; cell++)
            {  
                PFlowArray[(int) CellList[cell][0]][(int) CellList[cell][1]][(int) CellList[cell][2]]=0;
                CellList[cell][5] = 0.0;
                CellList[cell][4] = 0.0;
            }
        }
        G4cout<<"Rest_energy "<<Rest_energy<<" 1.5*sigmaEref "<<1.5*sigmaEref*mom<<G4endl;
         
        G4cout<<"D8"<<G4endl;
    }
    G4cout<<"SumEdep "<<SumEdep<<G4endl;
    pflowneutral = 0;
    for (int i = 0; i < 6; i++)
    {
        
        for (int j = 0; j < LayersPix[i]; j++)
        {
            for (int k = 0; k < LayersPix[i]; k++)
            {
                pflowneutral += PFlowArray[i][j][k];
            }
        }
    }
    G4cout<<"pflowneutralInPFinal "<<pflowneutral<<G4endl;
    G4cout<<G4endl;


}
