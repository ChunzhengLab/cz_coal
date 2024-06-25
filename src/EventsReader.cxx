#include <fstream>
#include <iostream>
#include "EventsReader.h"
#include "Particle.h"

EventsReader::EventsReader(EventType eventType, TString nameFile) : eventType(eventType) {
  InitTree(nameFile);
  nEvents = chain->GetEntries();
}

Event<Parton>& EventsReader::GetEvent(int iEvent) {
  UpdateEvent(iEvent);
  return currentEvent;
}

void EventsReader::UpdateEvent(int iEvent) {
  chain->GetEntry(iEvent);
  currentEvent = std::move(Event<Parton>(partonEventStruct));
  if (par::isDebug) {
    std::cout << ">>>>>>>>Move to Entry(iEvent)>>>>>>>" << iEvent << std::endl;
    std::cout << "Event ID: " << partonEventStruct.nevent << std::endl;
    std::cout << "First parton ID: " << partonEventStruct.ID[0] << std::endl;
    std::cout << "First parton Z: " << partonEventStruct.Z[0] << std::endl;
    std::cout << "Second parton ID: " << partonEventStruct.ID[1] << std::endl;
    std::cout << "Second parton Z: " << partonEventStruct.Z[1] << std::endl;
    std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
  }
}

// AMPT event reader
bool EventsReader::InitTree(TString nameFile) {
  // 确定tree的名称
  std::string treeName = "AMPT";
  switch (eventType) {
    case EventType::kAMPT: treeName = "AMPT";
    break;
    case EventType::kRandom: treeName = "RandPartons";
    break;
  }

  chain = std::unique_ptr<TChain>(new TChain(treeName.c_str()));
  int fileNum = 0;
  std::string fileList;
  if (nameFile.EndsWith(".list")) {
    std::ifstream file(nameFile.Data());
    if (!file.is_open()) {
        std::cerr << "Error: cannot open file " << nameFile << std::endl;
        return false;
    }
    while (std::getline(file, fileList)) {
      if (fileList.empty()) continue; // 跳过空行
      TFile *f_tmp = TFile::Open(fileList.c_str());
      if (!f_tmp || !f_tmp->IsOpen() || f_tmp->GetNkeys() == 0) {
          std::cerr << "Error: cannot open list file " << fileList << std::endl;
          delete f_tmp; // 清理临时文件对象
          return false;
      } else {
          chain->Add(fileList.c_str());
          fileNum++;
      }
      delete f_tmp; // 关闭并删除临时文件对象
    }
    std::cout << fileNum << " files read" << std::endl;
  } else if (nameFile.EndsWith(".root")) {
    chain->Add(nameFile.Data());
  } else {
    std::cerr << "Error: unknown file extension for " << nameFile << std::endl;
    return false;
  }

  // 设置分支地址
  chain->SetBranchAddress("Event", &partonEventStruct);
  chain->SetBranchAddress("ID", partonEventStruct.ID);
  chain->SetBranchAddress("Px", partonEventStruct.Px);
  chain->SetBranchAddress("Py", partonEventStruct.Py);
  chain->SetBranchAddress("Pz", partonEventStruct.Pz);
  chain->SetBranchAddress("X", partonEventStruct.X);
  chain->SetBranchAddress("Y", partonEventStruct.Y);
  chain->SetBranchAddress("Z", partonEventStruct.Z);

  // 调试信息
  if (par::isDebug) {
    std::cout << "-----------------" << std::endl;
    std::cout << "Initializing readin Tree" << std::endl;
    std::cout << "readin tree linked to partonEventStruct, but has not into any entry" << std::endl;
    std::cout << "readin tree has " << chain->GetEntries() << " entries" << std::endl;
    std::cout << "NOTE: Event ID starts from 1, but iEvent(entry) starts from 0" << std::endl;
    std::cout << "-----------------" << std::endl;
  }
  return true;
}






