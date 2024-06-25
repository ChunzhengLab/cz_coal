#ifndef EVENTSREADER_H
#define EVENTSREADER_H

#include <fstream>
#include "Par.h"
#include "Event.h"
#include "Particle.h"
#include "TFile.h"
#include "TChain.h"

class EventsReader {
  //这个类里面要实现读取事件的功能，分成两种模式，第一种是读取来随机抽样的事件
  private:
  //事件类型
  int nEvents;
  EventType eventType;
  PartonEventStruct partonEventStruct;
  std::unique_ptr<TChain> chain;
  Event<Parton> currentEvent;

  public:
  EventsReader(EventType eventType, TString inputFile);
  bool InitTree(TString inputFile);
  int GetNEvents() { return nEvents; }

  Event<Parton>& GetEvent(int iEvent);
  void UpdateEvent(int iEvent);
  void Print() const {
    std::cout <<"--------------------------" << std::endl;
    std::cout << "EventsReader:" << std::endl;
    std::cout << "EventsReader with " << nEvents << " events" << std::endl;
    std::cout <<"--------------------------" << std::endl;
  }

  EventsReader(const EventsReader&) = delete;
  EventsReader& operator=(const EventsReader&) = delete;
};

#endif





