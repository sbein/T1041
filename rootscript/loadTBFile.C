// simple tool to load a TB file and library

loadTBFile(TString name="latest.root"){
  if (!TClassTable::GetDict("TBEvent")) {
    gSystem->Load("$TBHOME/build/lib/libTB.so");  // n.b. make sure to compile if changed
  }
  TFile *f=new TFile(name);
  f->ls();
  TBrowser *tb=new TBrowser();
}
