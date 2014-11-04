#include "TBRecHit.h"
#include "Mapper.h"
#include <stdio.h>
#include <stdlib.h>

static const TString gainNameMap[4] = {"Low","Mid","High","VHigh"};
static const Int_t CUTAMP_PLV_LOWLOW100  =  300;  static const Int_t CUTChi2_PLV_LOWLOW100  =   80;
static const Int_t CUTAMP_PLV_LOWLOW200  = 9999;  static const Int_t CUTChi2_PLV_LOWLOW200  = 9999;
static const Int_t CUTAMP_PLV_LOWLOW300  =  650;  static const Int_t CUTChi2_PLV_LOWLOW300  =  200;
static const Int_t CUTAMP_PLV_LOWHIGH100 =  300;  static const Int_t CUTChi2_PLV_LOWHIGH100 =  300;
static const Int_t CUTAMP_PLV_LOWHIGH200 = 9999;  static const Int_t CUTChi2_PLV_LOWHIGH200 = 9999;
static const Int_t CUTAMP_PLV_LOWHIGH300 =  500;  static const Int_t CUTChi2_PLV_LOWHIGH300 =  500;
static const Int_t CUTAMP_PLV_MIDLOW100  =  500;  static const Int_t CUTChi2_PLV_MIDLOW100  =   80;
static const Int_t CUTAMP_PLV_MIDLOW200  = 9999;  static const Int_t CUTChi2_PLV_MIDLOW200  = 9999;
static const Int_t CUTAMP_PLV_MIDLOW300  =  800;  static const Int_t CUTChi2_PLV_MIDLOW300  =  200;
static const Int_t CUTAMP_PLV_MIDLOW500  = 1700;  static const Int_t CUTChi2_PLV_MIDLOW500  =  200;
static const Int_t CUTAMP_PLV_MIDHIGH100 =  450;  static const Int_t CUTChi2_PLV_MIDHIGH100 =  300;
static const Int_t CUTAMP_PLV_MIDHIGH200 =  700;  static const Int_t CUTChi2_PLV_MIDHIGH200 = 80;
static const Int_t CUTAMP_PLV_MIDHIGH300 = 1200;  static const Int_t CUTChi2_PLV_MIDHIGH300 = 2000;
static const Int_t CUTAMP_PLV_MIDHIGH500 = 1800;  static const Int_t CUTChi2_PLV_MIDHIGH500 = 2000;

TBRecHit::TBRecHit(PadeChannel *pc, Float_t zsp, UInt_t options){
  nzsp=zsp;
  if (options&kNoFit) {
    status|=kNoFit;
    Init(pc,false);
  }
  else Init(pc);
}

TBRecHit::TBRecHit(const TBRecHit &hit, UShort_t idx, UInt_t newstatus) :
  TObject(hit),
  maxADC(hit.maxADC), pedestal(hit.pedestal),
  noise(hit.noise), aMaxValue(hit.aMaxValue), tRiseValue(hit.tRiseValue),
  chi2(hit.chi2), ndof(hit.ndof), nzsp(hit.nzsp), status(hit.status|newstatus),
  cfactor(hit.cfactor), ts(hit.ts)
{
  if (idx<=127) channelIndex=idx;
}

void TBRecHit::Init(PadeChannel *pc,  Float_t zsp){
  channelIndex=-1;
  maxADC=-1;
  pedestal=-999;
  noise=-999;
  aMaxValue=-999;
  tRiseValue=-1;
  aMaxError=0;
  tRiseError=0;
  ndof=0;
  status=0;
  nzsp=zsp;
  cfactor=1;
  if (!pc) return;
  channelIndex=pc->GetChannelIndex();
  ts=pc->GetTimeStamp();
  maxADC=pc->GetMax();
  double ped,sig;
  pc->GetPedestal(ped,sig);
  pedestal=ped;
  noise=sig;
  if (status&kNoFit) return;
  FitPulse(pc);
}


void TBRecHit::GetXYZ(double &x, double &y, double &z) const {
  Mapper *mapper=Mapper::Instance(ts);
  mapper->ChannelIdxXYZ(channelIndex,x,y,z); // does not depend on run epoch
}

void TBRecHit::GetXYZ(float &x, float &y, float &z) const {
  double p[3];
  GetXYZ(p[0],p[1],p[2]);
  x=p[0];
  y=p[1];
  z=p[2];
}


Int_t TBRecHit::GetChannelID() const{  // depends on run epoch b/c of PADE id's
  Mapper *mapper=Mapper::Instance(ts);
  return mapper->ChannelIndex2ChannelID(channelIndex);
}

void TBRecHit::GetModuleFiber(int &moduleID, int &fiberID) const{
  Mapper *mapper=Mapper::Instance(ts);  
  mapper->ChannelIndex2ModuleFiber(channelIndex,moduleID,fiberID); // result independent of run epoch
}

void TBRecHit::FitPulse(PadeChannel *pc){
  if ( TMath::Abs(maxADC-pedestal) / (noise+0.001) < nzsp ) { // avoid div by 0
    status|=kZSP;
    return;
  }
  PulseFit fit=PadeChannel::FitPulse(pc);
  pedestal=fit.pedestal;
  noise=fit.noise;
  aMaxValue=fit.aMaxValue;
  tRiseValue=fit.tRiseValue;
  aMaxError=fit.aMaxError;
  tRiseError=fit.tRiseError;
  chi2=fit.chi2Peak;
  ndof=fit.ndofPeak;
  if (fit.status>0) status|=kPoorFit;
}

std::ostream& operator<<(std::ostream& s, const TBRecHit& hit) {
  double x,y,z;
  hit.GetXYZ(x,y,z);
  return s << "TBRecHit (index,x,y,z,noise,aMax,tRise,Chi2) " 
	   << hit.ChannelIndex() << "," << x << "," << y << "," << z<< ","
	   << hit.NoiseRMS() << "," << hit.AMax() << "," 
	   << hit.TRise() << "," << hit.Chi2();
}

void TBRecHit::Calibrate(float *calconstants){
  float c=calconstants[channelIndex];
  if (IsCalibrated()) c/=cfactor;  // undo pevious calibration
  aMaxValue*=c;
  aMaxError*=c;
  noise*=c;
  cfactor=calconstants[channelIndex];
}
void TBRecHit::Calibrate(vector<TBRecHit> *rechits, float *calconstants){
  for (unsigned i=0; i<rechits->size(); i++)
    (rechits->at(i)).Calibrate(calconstants);
}



Bool_t TBRecHit::GoodPulse(PadeChannel* pc, UShort_t pga, UShort_t lna, ULong_t vga) {

  bool isGood_ = true;
  double max = 9999; double chi2 = 99999; double ndof = 1;

  int nSigmaCut=30;
  Init(pc, nSigmaCut);
  if (!(Status() && kZSP))
    {
      max = pc->GetMaxCalib();
      chi2= Chi2();
      ndof= Ndof();
    }

  //  PadeChannel TBRecHit::GetPade(){
  //  return pc

  TString Pga = gainNameMap[pga];  TString Lna = gainNameMap[lna]; 
  TString Vga=TString::Itoa(vga,16);
  TString pga_lna_vga = Pga+"_"+Lna+"_"+Vga;
  cout << "pga="<<pga<< " lna="<<lna<< " vga="<<vga<<endl;
  cout << "pga_lna_vga is "<< pga_lna_vga<<endl;
  if      (pga_lna_vga=="LOW_LOW_100")
    isGood_ = max < CUTAMP_PLV_LOWLOW100 && 1.0*chi2/ndof < CUTChi2_PLV_LOWLOW100;
  else if (pga_lna_vga=="LOW_LOW_200")
    isGood_ = max < CUTAMP_PLV_LOWLOW200 && 1.0*chi2/ndof < CUTChi2_PLV_LOWLOW200;
  else if (pga_lna_vga=="LOW_LOW_300")
    isGood_ = max < CUTAMP_PLV_LOWLOW300 && 1.0*chi2/ndof < CUTChi2_PLV_LOWLOW300;
  else if (pga_lna_vga=="LOW_HIGH_100")
    isGood_ = max < CUTAMP_PLV_LOWHIGH100 && 1.0*chi2/ndof < CUTChi2_PLV_LOWHIGH100;
  else if (pga_lna_vga=="LOW_HIGH_200")
    isGood_ = max < CUTAMP_PLV_LOWHIGH200 && 1.0*chi2/ndof < CUTChi2_PLV_LOWHIGH200;
  else if (pga_lna_vga=="LOW_HIGH_300")
    isGood_ = max < CUTAMP_PLV_LOWHIGH300 && 1.0*chi2/ndof < CUTChi2_PLV_LOWHIGH300;
  else if (pga_lna_vga=="MID_LOW_100")
    isGood_ = max < CUTAMP_PLV_MIDLOW100 && 1.0*chi2/ndof < CUTChi2_PLV_MIDLOW100;
  else if (pga_lna_vga=="MID_LOW_200")
    isGood_ = max < CUTAMP_PLV_MIDLOW200 && 1.0*chi2/ndof < CUTChi2_PLV_MIDLOW200;
  else if (pga_lna_vga=="MID_LOW_300")
    isGood_ = max < CUTAMP_PLV_MIDLOW300 && 1.0*chi2/ndof < CUTChi2_PLV_MIDLOW300;
  else if (pga_lna_vga=="MID_LOW_500")
    isGood_ = max < CUTAMP_PLV_MIDLOW500 && 1.0*chi2/ndof < CUTChi2_PLV_MIDLOW500;
  else if (pga_lna_vga=="MID_HIGH_100")
    isGood_ = max < CUTAMP_PLV_MIDHIGH100 && 1.0*chi2/ndof < CUTChi2_PLV_MIDHIGH100;
  else if (pga_lna_vga=="MID_HIGH_200")
    isGood_ = max < CUTAMP_PLV_MIDHIGH200 && 1.0*chi2/ndof < CUTChi2_PLV_MIDHIGH200;
  else if (pga_lna_vga=="MID_HIGH_300")
    isGood_ = max < CUTAMP_PLV_MIDHIGH300 && 1.0*chi2/ndof < CUTChi2_PLV_MIDHIGH300;
  else if (pga_lna_vga=="MID_HIGH_500")
    isGood_ = max < CUTAMP_PLV_MIDHIGH500 && 1.0*chi2/ndof < CUTChi2_PLV_MIDHIGH500;
  return isGood_;
}
