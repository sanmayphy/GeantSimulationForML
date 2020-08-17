#include "DownResolution.hh"


DownResolution::DownResolution(const float(&TotalEnergyArray)[6][kMaxPixel][kMaxPixel], const float (&NoiseArray)[6][kMaxPixel][kMaxPixel], const float (&ChargeEnergyArray)[6][kMaxPixel][kMaxPixel], const float (&NeutralEnergyArray)[6][kMaxPixel][kMaxPixel]) 
{
    for (int l = 0; l < 6; l++)
    {
        SumPix(l, TotalEnergyArray[l], NoiseArray[l], ChargeEnergyArray[l], NeutralEnergyArray[l]);
    }
}

void DownResolution::SumPix(int layer, const float (&TotalEnergyArray)[kMaxPixel][kMaxPixel], const float (&NoiseArray)[kMaxPixel][kMaxPixel], const float (&ChargeEnergyArray)[kMaxPixel][kMaxPixel], const float (&NeutralEnergyArray)[kMaxPixel][kMaxPixel]) 
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
            NoiseArrayGran[layer][e][p] = 0;
            TotalEnergyGran[layer][e][p] = 0;
            NeutralEnergyGran[layer][e][p] = 0;
            ChargeEnergyGran[layer][e][p] = 0;
            OnlyNoiseCell = true && (TotalEnergyArray[e][p]==NoiseArray[e][p]);
        }
    }
    int scale = kMaxPixel/size;
    //* Fill CellArray with in case of diffrent granularity
    for (int it_x = 0; it_x < size; it_x++)
    {
        for (int it_y = 0; it_y < size; it_y++)
        {
            float CellTotSignal=  0;//NoiseArray[it_x][it_y];
            float CellChSignal = 0;//CellTotSignal;
            float CellNuSignal = 0;//CellTotSignal;


            NoiseArrayGran[layer][it_x][it_y] =  NoiseArray[it_x][it_y];
            // if ( !OnlyNoiseCell)
            // {
                for (int ii = scale*it_x; ii < scale*(it_x+1); ii++ )
                {
                    for (int jj = scale*it_y; jj < scale*(it_y+1); jj++ )
                    {
                        CellTotSignal += TotalEnergyArray[ii][jj];
                        CellChSignal += ChargeEnergyArray[ii][jj];
                        CellNuSignal += NeutralEnergyArray[ii][jj];
                    }
                }
            // }
            // CellArray[layer][it_x][it_y].SetSigma(CellTotSignal/NoiseValues[layer]);

            TotalEnergyGran[layer][it_x][it_y] = CellTotSignal;
            NeutralEnergyGran[layer][it_x][it_y] = CellNuSignal;
            ChargeEnergyGran[layer][it_x][it_y] = CellChSignal;
        }
    }
}

void DownResolution::FillTotEn(float (&TotalEnergyArray)[6][kMaxPixel][kMaxPixel])
{
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < kMaxPixel; j++)
        {
            for (int k = 0; k < kMaxPixel; k++)
            {
                TotalEnergyArray[i][j][k] = TotalEnergyGran[i][j][k];
            }
        }
    }
}

void DownResolution::FillChEn(float (&ChargeEnergyArray)[6][kMaxPixel][kMaxPixel])
{
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < kMaxPixel; j++)
        {
            for (int k = 0; k < kMaxPixel; k++)
            {
                ChargeEnergyArray[i][j][k] = ChargeEnergyGran[i][j][k];
            }
        }
    }
}

void DownResolution::FillNuEn(float (&NeutralEnergyArray)[6][kMaxPixel][kMaxPixel])
{
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < kMaxPixel; j++)
        {
            for (int k = 0; k < kMaxPixel; k++)
            {
                NeutralEnergyArray[i][j][k] = NeutralEnergyGran[i][j][k];
            }
        }
    }
}

void DownResolution::FillNoise(float (&Noise)[6][kMaxPixel][kMaxPixel])
{
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < kMaxPixel; j++)
        {
            for (int k = 0; k < kMaxPixel; k++)
            {
                Noise[i][j][k] = NoiseArrayGran[i][j][k];
            }
        }
    }

}
