#include "GlobalVariables.hh"


class DownResolution
{
    public:
        DownResolution(const float (&TotalEnergyArray)[6][kMaxPixel][kMaxPixel], const float (&NoiseArray)[6][kMaxPixel][kMaxPixel], const float (&ChargeEnergyArray)[6][kMaxPixel][kMaxPixel], const float (&NeutralEnergyArray)[6][kMaxPixel][kMaxPixel]);
        void FillTotEn(float (&TotalEnergyArray)[6][kMaxPixel][kMaxPixel]);
        void FillChEn(float (&ChargeEnergyArray)[6][kMaxPixel][kMaxPixel]);
        void FillNuEn(float (&NeutralEnergyArray)[6][kMaxPixel][kMaxPixel]);
        void FillNoise(float (&NoiseArray)[6][kMaxPixel][kMaxPixel]);
    private:
        float TotalEnergyGran[6][kMaxPixel][kMaxPixel];
        float ChargeEnergyGran[6][kMaxPixel][kMaxPixel];
        float NeutralEnergyGran[6][kMaxPixel][kMaxPixel];
        float NoiseArrayGran[6][kMaxPixel][kMaxPixel];
        void SumPix(int layer, const float (&TotalEnergyArray)[kMaxPixel][kMaxPixel], const float (&NoiseArray)[kMaxPixel][kMaxPixel], const float (&ChargeEnergyArray)[kMaxPixel][kMaxPixel], const float (&NeutralEnergyArray)[kMaxPixel][kMaxPixel]);


};