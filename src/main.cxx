#include "EventsReader.h"
#include "EventsWriter.h"
#include "Event.h"
#include "Coalescence.h"
#include "CalculateObvs.h"
#include "Par.h"


int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <config-file>" << std::endl;
    return 1;
  }
  std::ifstream configFile(argv[1]);
  if (!configFile) {
    std::cerr << "Unable to open config file: " << argv[1] << std::endl;
    return 1;
  }
  std::string line;
  while (std::getline(configFile, line)) {
    try {
      parseConfig(line);
    } catch (const std::exception& e) {
      std::cerr << e.what() << std::endl;
      return 1;
    }
  }

  printConfig();

  //如果inputfile无法打开，抛出异常
  try {
    TFile file(par::inputFile.c_str(), "READ");
    if (!file.IsOpen() || file.IsZombie()) {
      throw std::runtime_error("Error: cannot open input file: " + par::inputFile);
    }
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
    return 1;
  }

  // 初始化读取器，写入器，聚合器以及观测量计算器
  EventsReader reader(par::eventType, par::inputFile);
  reader.Print();
  Coalescence coalescence(par::r_bm, par::coalescenceAlgorithm);
  coalescence.Print();
  CalculateObvs calculateObvs(par::obvsFile);
  calculateObvs.Print();
  EventsWriter writer(par::outputFile);
  writer.Print();

  int nEvents = reader.GetNEvents();
  if (par::isDebug) {
    std::cout << "Total number of events: " << nEvents << std::endl;
    std::cout << "For debug mode, only process the first 20 events" << std::endl;
    nEvents = 20;
  }

  std::cout << "Start processing events" << std::endl;
  for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
    // 进度条
    if (iEvent % 50 == 0) {
      std::cout << "Processing event " << iEvent << " / " << reader.GetNEvents() << std::endl;
    }
    
    std::cout << "======================" << std::endl;
    // 读取PartonEvent
    Event<Parton>& partonEvent = reader.GetEvent(iEvent);

    // 执行Coalescence
    std::vector<Hadron> hadrons;
    coalescence.Process(partonEvent.GetParticles(), hadrons);
    if (par::isDebug) {
      std::cout << "Hadrons after coalescence: " << hadrons.size() << std::endl;
    }

    //执行观测量计算
    if(par::isCalculateObvs) {
      calculateObvs.Process(hadrons);
    }

    //写入HadronEvent
    if(par::isWriteEvents) {
      writer.WriteEvent(Event<Hadron>(iEvent + 1, hadrons.size(), std::move(hadrons)));
    }
    std::cout << "======================" << std::endl;
  }
  std::cout << "All events processed" << std::endl;

  return 0;
}


// 考虑是否要将最后剩下的quark进行一次尝试聚合，如果聚合成功，那么就将其加入到hadrons中，否则就不加入。但是目前似乎没有那么紧张



// Insights:
// 检查一下AMPT 和 Random 的事件中 px和py与x和y的关系，提取所谓的 Radial Flow
// 我猜测AMPT事件有Radial Flow，而Random事件没有Radial Flow