#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1I.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <random>
#include <iostream>

using namespace std;

struct Event {
  int nevent;
  int nparton;
  int ID[100000]; // [nparton]
  float Px[100000]; // [nparton]
  float Py[100000]; // [nparton]
  float Pz[100000]; // [nparton]
  float X[100000]; // [nparton]
  float Y[100000]; // [nparton]
  float Z[100000]; // [nparton]
};

enum EventType {
  kPureRandom ,
  kSampleFromAMPT ,
};

int main(int argc , char* argv[]) {
  int nEventsType = kPureRandom;
  int nPartonsYouWant = 0;
  TString amptPartonDistFileName;
  TString randPartonsFileName;

  if (argc == 5) {
    nPartonsYouWant = atoi(argv[1]);
    nEventsType = atoi(argv[2]);
    amptPartonDistFileName = argv[3];
    randPartonsFileName = argv[4];
  }
  else {
    cerr << "Usage: " << argv[0] << " <nPartonsYouWant> <AMPT parton distribution file> <output file name>" << endl;
    return 1;
  }

  // PDG code map bin 0 -> u, 1 -> ubar, 2 -> d, 3 -> dbar, 4 -> s, 5 -> sbar, 6 -> c, 7 -> cbar, 8 -> b, 9 -> bbar, 10 -> other positive, 11 -> other negative
  const int nEvents = 500;
  const int pdgMap[12] = { 2, -2, 1, -1, 3, -3, 4, -4, 5, -5, 6, -6 };
  const float R = 10.; // 半径
  const float pT = 1.5; // 横向动量
  gRandom = new TRandom3(time(nullptr)); // 使用时间作为随机数种子

  cout << "Reading AMPT parton distribution from " << amptPartonDistFileName << endl;
  auto inputFile = new TFile(amptPartonDistFileName , "READ");
  if (!inputFile || inputFile->IsZombie()) {
    cerr << "Cannot open file " << amptPartonDistFileName << endl;
    return 1;
  }
  auto h2_xy = dynamic_cast<TH2D*>(inputFile->Get("h2_xy"));
  if (!h2_xy) {
    cerr << "Cannot find h2_xy" << endl;
    return 1;
  }
  auto h_z = dynamic_cast<TH1D*>(inputFile->Get("h_z"));
  if (!h_z) {
    cerr << "Cannot find h_z" << endl;
    return 1;
  }
  auto h2_pxpy = dynamic_cast<TH2D*>(inputFile->Get("h2_pxpy"));
  if (!h2_pxpy) {
    cerr << "Cannot find h2_pxpy" << endl;
    return 1;
  }
  auto h_pz = dynamic_cast<TH1D*>(inputFile->Get("h_pz"));
  if (!h_pz) {
    cerr << "Cannot find h_pz" << endl;
    return 1;
  }
  auto h_pdg = dynamic_cast<TH1I*>(inputFile->Get("h_pdg"));
  if (!h_pdg) {
    cerr << "Cannot find h_pdg" << endl;
    return 1;
  }
  auto h_nparton = dynamic_cast<TH1I*>(inputFile->Get("h_nparton"));
  if (!h_nparton) {
    cerr << "Cannot find h_nparton" << endl;
    return 1;
  }

  auto outputFile = new TFile(randPartonsFileName , "RECREATE");
  auto outputTree = new TTree("RandPartons" , "Random Generator according to AMPT Parton Distribution");

  Event event;
  outputTree->Branch("Event" , &event , "nevent/I:nparton/I");
  outputTree->Branch("ID" , event.ID , "ID[nparton]/I");
  outputTree->Branch("Px" , event.Px , "Px[nparton]/F");
  outputTree->Branch("Py" , event.Py , "Py[nparton]/F");
  outputTree->Branch("Pz" , event.Pz , "Pz[nparton]/F");
  outputTree->Branch("X" , event.X , "X[nparton]/F");
  outputTree->Branch("Y" , event.Y , "Y[nparton]/F");
  outputTree->Branch("Z" , event.Z , "Z[nparton]/F");

  cout << "Generating " << nEvents << " events with " << nPartonsYouWant << " partons each." << endl;

  for (int iEvt = 0; iEvt < nEvents; iEvt++) {
    if (iEvt % 100 == 0) cout << "Event " << iEvt << endl;
    event.nevent = iEvt + 1;
    event.nparton = nPartonsYouWant;
    double x{ 0 } , y{ 0 } , z{ 0 };
    double px{ 0 } , py{ 0 } , pz{ 0 };
    int pdg{ 0 };

    for (int iTrk = 0; iTrk < event.nparton; iTrk++) {

      switch (nEventsType) {
      case kPureRandom: {
        // 均匀球体分布
        double u1 = gRandom->Rndm();
        double u2 = gRandom->Rndm();
        double u3 = gRandom->Rndm();
        double r = R * std::cbrt(u1);  // 半径
        double theta = 2 * M_PI * u2;  // 方位角
        double phi = std::acos(2 * u3 - 1);  // 极角
        x = r * std::sin(phi) * std::cos(theta);
        y = r * std::sin(phi) * std::sin(theta);
        z = r * std::cos(phi);

        // 均匀圆环分布，为了简化, pz取0GeV/c,pT取1.5GeV/c
        double phi_p = 2 * M_PI * gRandom->Rndm();
        px = pT * std::cos(phi_p);
        py = pT * std::sin(phi_p);
        pz = 0;
        break;
      }

      case kSampleFromAMPT: {
        h2_xy->GetRandom2(x , y);
        z = h_z->GetRandom();
        h2_pxpy->GetRandom2(px , py);
        pz = h_pz->GetRandom();
        break;
      }

      default:
        cerr << "Unknown event type" << endl;
      }


      pdg = pdgMap[(int)floor(h_pdg->GetRandom())];
      event.ID[iTrk] = pdg;
      event.Px[iTrk] = px;
      event.Py[iTrk] = py;
      event.Pz[iTrk] = pz;
      event.X[iTrk] = x;
      event.Y[iTrk] = y;
      event.Z[iTrk] = z;
    }
    outputTree->Fill();
  }


  outputFile->cd();
  cout << "Writing random partons to " << randPartonsFileName << endl;
  outputTree->Write();
  outputFile->Close();
  inputFile->Close(); // Close input file
  delete inputFile; // Ensure input file is properly deleted
  delete outputFile; // Ensure output file is properly deleted

  return 0;
}