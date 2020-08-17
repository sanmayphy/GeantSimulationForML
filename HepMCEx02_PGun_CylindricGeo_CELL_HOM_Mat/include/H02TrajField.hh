#include "GlobalVariables.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4NistManager.hh"
#include "cmath"
using std::vector;
using std::sort;

class ChargeParticlTraj
{
    public:
        ChargeParticlTraj(vector < vector <double> > &RealParticalInfo,vector < vector <double> > &TrajList, float (&Track_Cell)[kMaxPixel][kMaxPixel]);
    private:
        void TrajFinder(vector <double> &Particle, vector < vector <double> > &TrajList, float (&Track_Cell)[kMaxPixel][kMaxPixel]);
        template <typename T> int sgn(T val);
};
template <typename T> int ChargeParticlTraj::sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}
bool sortcol( const vector<double>& v1, const vector<double>& v2 ) 
{ 
    return v1[1] > v2[1]; 
} 
ChargeParticlTraj::ChargeParticlTraj(vector < vector <double> > &RealParticalInfo,vector < vector <double> > &TrajList, float (&Track_Cell)[kMaxPixel][kMaxPixel])
{
    int SizeRealParticalInfo = RealParticalInfo.size();
    for (int nParticle = 0; nParticle < SizeRealParticalInfo; nParticle++)
    {
        TrajFinder(RealParticalInfo[nParticle], TrajList, Track_Cell);
    }

    sort(TrajList.begin(), TrajList.end(),sortcol);

}

void ChargeParticlTraj::TrajFinder(vector <double> &Particle,  vector < vector <double> > &TrajList, float (&Track_Cell)[kMaxPixel][kMaxPixel])
{
    auto DefParticle =  G4ParticleTable::GetParticleTable()->FindParticle(Particle[0]);
    double Hmagnetic = fieldValue;//3.8*tesla;//* tesla
    double mass = DefParticle->GetPDGMass();//* Mev
    double charge = DefParticle->GetPDGCharge();//* 1 elementary charge
    
    if (charge!=0)
    {
        // G4cout<<"GGGGGGGG0"<<G4endl;
        if (Hmagnetic!=0)
        {
            double energy = Particle[4];
            double LighSpeed = 299792458 * m/s;//* mm/ ns
            double px = Particle[1];//* Mev
            double py = Particle[2];//* Mev
            double pz = Particle[3];//* Mev
            double x0 = py/(charge*Hmagnetic*LighSpeed);//gamma
            double y0 = -1*px/(charge*Hmagnetic*LighSpeed);//gamma
            double alpha = sgn(charge*Hmagnetic)*r_inn_ECAL1/(2*pow(pow(x0,2)+pow(y0,2),0.5));//r_inn_ECAL1
            if (abs(alpha)<1)
            {
                double pt = sqrt(sqr(px)+sqr(py));
                double px_f = (px*(1-2*sqr(alpha))+2*py*alpha*sqrt(1-sqr(alpha)));
                double py_f = (py*(1-2*sqr(alpha))-2*px*alpha*sqrt(1-sqr(alpha)));
                // double pz_f = ;

                // G4cout<<"Energy= "<<energy<<G4endl;
                // G4cout<<"Alpha= "<<alpha<<G4endl;
                int IndexEta = -1;
                int IndexPhi = -1;

                vector <float> phi(6,0);
                vector <float> eta(6,0);
                vector <float> IndexPhi_v(6,-1);
                vector <float> IndexEta_v(6,-1);
                float R_inn[6] = {r_inn_ECAL1,r_inn_ECAL2,r_inn_ECAL3,r_inn_HCAL1,r_inn_HCAL2,r_inn_HCAL3};

                double x = 2*(x0*pow(alpha,2)-y0*alpha*pow(1-pow(alpha,2),0.5));
                double y = 2*(y0*pow(alpha,2)+x0*alpha*pow(1-pow(alpha,2),0.5));
                double z = (pz)/(abs(charge*Hmagnetic)*LighSpeed)*abs(2*asin(alpha));//gamma
                double r2 = pow(x,2)+pow(y,2);
                // G4cout<<"r2 "<<r2<<" r_inn_ECAL1 "<<sqr(r_inn_ECAL1)<<G4endl;
                double a = sqr(pt/mass);
                double b = 2*(px_f*x+py_f*y)/mass; 
                // G4cout<<"G1"<<G4endl;
                for (int lay = 0; lay < 6; lay++)
                {
                    // G4cout<<"G2"<<G4endl;
                    float c = sqr(R_inn[lay])-r2;
                    float t = (-b+sqrt(4*a*c+sqr(b)))/(2*a);
                    double x_f = x+(px/mass)*t;
                    double y_f = y+(py/mass)*t;
                    double z_f = z+(pz/mass)*t;
                    double phi_ = atan2(y_f,x_f);
                    // G4cout<<"G2"<<G4endl;
                    if (phi_<0) phi_+=tube_dPhi;
                    phi[lay]=phi_;
                    double eta_ = -1*log(tan(0.5*acos(z_f/sqrt(sqr(R_inn[lay])+sqr(z_f)))));
                    // G4cout<<"G3"<<G4endl;
                    eta[lay]=eta_;
                    float deta = d_eta*(kMaxPixel/LayersPix[lay]);
                    float dphi = divided_tube_dPhi*(kMaxPixel/LayersPix[lay]);
                    IndexPhi = (float)floor(phi_/dphi);
                    IndexEta = (float)floor((eta_max+eta_)/deta);
                    if ((IndexPhi>LayersPix[lay]) || (IndexEta>LayersPix[lay]) ) 
                    {
                        IndexPhi=-1;
                        IndexEta=-1;
                    }
                    IndexPhi_v[lay] = IndexPhi;
                    IndexEta_v[lay] = IndexEta;
                }
                if ((IndexEta_v[0]!=-1 && IndexPhi_v[0]!=-1) && pt<40e+03)
                {
                    // G4cout<<"G5"<<G4endl;
                    vector <double> TrajOfParticl;
                    float momentum = sqrt(sqr(px)+sqr(py)+sqr(pz));
                    // G4cout<<"G6"<<G4endl;
                    TrajOfParticl.push_back(Particle[0]);
                    TrajOfParticl.push_back(energy);
                    TrajOfParticl.push_back(momentum);//momentum
                    TrajOfParticl.push_back(pt);

                    for (int lay = 0; lay<6;lay++)
                    {
                        TrajOfParticl.push_back(eta[lay]);
                        TrajOfParticl.push_back(phi[lay]);
                        TrajOfParticl.push_back(IndexEta_v[lay]);
                        TrajOfParticl.push_back(IndexPhi_v[lay]);
                    }
                    // G4cout<<"G7"<<G4endl;
                    TrajOfParticl.push_back(100000);
                    TrajOfParticl.push_back(0);
                    

                    // G4cout<<"IndexEta "<<IndexEta<<" IndexPhi "<<IndexPhi<<G4endl;
                    // G4cout<<"G8"<<G4endl;
                    TrajList.push_back(TrajOfParticl);
                    // G4cout<<"G9"<<G4endl;
                    // G4cout<<(int)IndexEta_v[1]<<" "<<(int)IndexPhi_v[1]<<G4endl;
                    Track_Cell[(int)IndexEta_v[1]][(int)IndexPhi_v[1]] += energy;
                    // G4cout<<"G10"<<G4endl;
                }
            }
        
        }
        else
        {
            // G4cout<<"GGGGGGGG1"<<G4endl;
            vector <float> phi(6,0);
            vector <float> eta(6,0);
            vector <float> IndexPhi_v(6,-1);
            vector <float> IndexEta_v(6,-1);
            int IndexEta = -1;
            int IndexPhi = -1;
            double px = Particle[1];//* Mev
            double py = Particle[2];//* Mev
            double pz = Particle[3];//* Mev
            double pt = sqrt(sqr(px)+sqr(py));
            float R_inn[6] = {r_inn_ECAL1,r_inn_ECAL2,r_inn_ECAL3,r_inn_HCAL1,r_inn_HCAL2,r_inn_HCAL3};
            for (int lay = 0; lay < 6; lay++)
            {
                double t = R_inn[lay]/(pt/mass);
                float z = pz/mass*t;
                float x = px/mass*t;
                float y = py/mass*t;
                eta[lay] = -1*log(tan(0.5*acos(z/sqrt(sqr(R_inn[lay])+sqr(z)))));
                phi[lay] = atan2(y,x);
                if (phi[lay]<0) phi[lay]+=tube_dPhi;
                float deta = d_eta*(kMaxPixel/LayersPix[lay]);
                float dphi = divided_tube_dPhi*(kMaxPixel/LayersPix[lay]);
                IndexPhi = (float)floor(phi[lay]/dphi);
                IndexEta = (float)floor((eta_max+eta[lay])/deta);
                if ((IndexPhi>LayersPix[lay]) || (IndexEta>LayersPix[lay]) ) 
                {
                    IndexPhi=-1;
                    IndexEta=-1;
                }
                IndexPhi_v[lay] = IndexPhi;
                IndexEta_v[lay] = IndexEta;
            }
            if ((IndexEta_v[5]!=-1 && IndexPhi_v[5]!=-1) && sqrt(sqr(px)+sqr(py))<40e+03)
            {
                // G4cout<<"GGGGGGGG2"<<G4endl;
                double energy = Particle[4];
                vector <double> TrajOfParticl;
                float momentum = sqrt(sqr(px)+sqr(py)+sqr(pz));
                TrajOfParticl.push_back(Particle[0]);
                TrajOfParticl.push_back(energy);
                TrajOfParticl.push_back(momentum);
                TrajOfParticl.push_back(pt);
                for (int lay = 0; lay<6;lay++)
                {
                    TrajOfParticl.push_back(eta[lay]);
                    TrajOfParticl.push_back(phi[lay]);
                    TrajOfParticl.push_back(IndexEta_v[lay]);
                    TrajOfParticl.push_back(IndexPhi_v[lay]);
                }
                TrajOfParticl.push_back(100000);
                TrajOfParticl.push_back(0);
                TrajList.push_back(TrajOfParticl);
                Track_Cell[(int)IndexEta_v[1]][(int)IndexPhi_v[1]] += energy;
            }
        }
        
    }
    




}