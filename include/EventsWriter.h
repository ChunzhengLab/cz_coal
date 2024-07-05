#ifndef EVENTS_WRITER_H
#define EVENTS_WRITER_H

#include <string>
#include <TFile.h>
#include <TTree.h>
#include "Event.h"
#include "Par.h"

class EventsWriter {
  //Write到一个root文件的TTree里面

  private:
    std::unique_ptr<TFile> outputFile;
    std::unique_ptr<TTree> tree;
    HadronEventStruct hadronEventStruct;

  public:
    explicit EventsWriter(const TString& filename) {
      outputFile = std::unique_ptr<TFile>(new TFile(filename, "RECREATE"));
      tree = std::unique_ptr<TTree>(new TTree("CoalHadrons", "CoalHadrons"));
      // Initialize tree branches once
      tree->Branch("nSeries", &hadronEventStruct.nSeries, "nSeries/I");
      tree->Branch("nTracks", &hadronEventStruct.nTracks, "nTracks/I");
      tree->Branch("PDG", &hadronEventStruct.PDG);
      tree->Branch("Pt", &hadronEventStruct.Pt);
      tree->Branch("Eta", &hadronEventStruct.Eta);
      tree->Branch("Phi", &hadronEventStruct.Phi);
      tree->Branch("x", &hadronEventStruct.x);
      tree->Branch("y", &hadronEventStruct.y);
      tree->Branch("z", &hadronEventStruct.z);
      tree->Branch("time", &hadronEventStruct.time);
      tree->Branch("dis", &hadronEventStruct.dis);
      tree->Branch("quark0", &hadronEventStruct.quark0);
      tree->Branch("quark1", &hadronEventStruct.quark1);
      tree->Branch("quark2", &hadronEventStruct.quark2);
      if (par::isDebug) {
        tree->Branch("quark0_x", &hadronEventStruct.quark0_x);
        tree->Branch("quark0_y", &hadronEventStruct.quark0_y);
        tree->Branch("quark0_z", &hadronEventStruct.quark0_z);
        tree->Branch("quark1_x", &hadronEventStruct.quark1_x);
        tree->Branch("quark1_y", &hadronEventStruct.quark1_y);
        tree->Branch("quark1_z", &hadronEventStruct.quark1_z);
        tree->Branch("quark2_x", &hadronEventStruct.quark2_x);
        tree->Branch("quark2_y", &hadronEventStruct.quark2_y);
        tree->Branch("quark2_z", &hadronEventStruct.quark2_z);
      }
    }
    ~EventsWriter() {
      outputFile->Write();
      // outputFile->Close();
    }

    void Event2HadronEventStruct(const Event<Hadron>& event);
    void WriteEvent(const Event<Hadron>& event);

    void Print() const {
      std::cout <<"--------------------------" << std::endl;
      std::cout << "EventsWriter:" << std::endl;
      std::cout << "EventsWriter with output file " << outputFile->GetName() << std::endl;
      std::cout <<"--------------------------" << std::endl;
    }
};


#endif // EVENTS_WRITER_H