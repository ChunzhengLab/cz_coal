#include "Par.h"
#include <stdexcept>

namespace par {
  bool isDebug = true;
  bool isWriteEvents = true;
  bool isCalculateObvs = true;
  bool isRemoveHFQuarks = true;

  std::string inputFile = "zpc-1.root";
  std::string outputFile = "output.root";
  std::string obvsFile = "obvs.root";

  //事件类型
  EventType eventType = EventType::kRandom;
  //聚合算法
  CoalescenceAlgorithm coalescenceAlgorithm = CoalescenceAlgorithm::kFromParton;
  //r_bm的默认值
  float r_bm = 1.0;
  float flavourBreakTolerance = 0.;

  //维护一个全局的随机数生成器
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> zero_or_one(0, 1);
}

void parseConfig(const std::string& line) {
  std::istringstream iss(line);
  std::string key, value;
  size_t eq_pos = line.find('=');
  if (eq_pos == std::string::npos) {
      throw std::runtime_error("Invalid configuration line: " + line);
  }
  key = line.substr(0, eq_pos);
  value = line.substr(eq_pos + 1);
  // 删除key和value周围的空格
  key.erase(key.find_last_not_of(" \t\n\r\f\v") + 1);
  value.erase(0, value.find_first_not_of(" \t\n\r\f\v"));

  // 下面是原来的判断逻辑
  if (key == "isDebug") par::isDebug = (value == "true");
  else if (key == "isWriteEvents") par::isWriteEvents = (value == "true");

  if (key == "isDebug") par::isDebug = (value == "true");
  else if (key == "isWriteEvents") par::isWriteEvents = (value == "true");
  else if (key == "isCalculateObvs") par::isCalculateObvs = (value == "true");
  else if (key == "isRemoveHFQuarks") par::isRemoveHFQuarks = (value == "true");
  else if (key == "eventType") {
      if (value == "kAMPT") par::eventType = EventType::kAMPT;
      else if (value == "kRandom") par::eventType = EventType::kRandom;
      else throw std::runtime_error("Unknown event type: " + value);
  }
  else if (key == "r_bm") par::r_bm = std::stof(value);
  else if (key == "flavourBreakTolerance") par::flavourBreakTolerance = std::stof(value);
  else if (key == "coalescenceAlgorithm") {
      if (value == "kClassic") par::coalescenceAlgorithm = CoalescenceAlgorithm::kClassic;
      else if (value == "kFromParton") par::coalescenceAlgorithm = CoalescenceAlgorithm::kFromParton;
      else throw std::runtime_error("Unknown coalescence algorithm: " + value);
  }
  else if (key == "inputFile") par::inputFile = value;
  else if (key == "outputFile") par::outputFile = value;
  else if (key == "obvsFile") par::obvsFile = value;
  else throw std::runtime_error("Unknown key: " + key);

  // for (int i = 0; i < 10; ++i) {
  //   int sample = par::zero_or_one(par::gen);
  //   std::cout << sample << " ";
  // }
}

void printConfig () {
  std::cout <<"--------------------------" << std::endl;
  std::cout << "Configuration: " << std::endl;
  std::cout << "isDebug = " << par::isDebug << std::endl;
  std::cout << "isWriteEvents = " << par::isWriteEvents << std::endl;
  std::cout << "isCalculateObvs = " << par::isCalculateObvs << std::endl;
  std::cout << "eventType = " << (par::eventType == EventType::kAMPT ? "kAMPT" : "kRandom") << std::endl;
  std::cout << "r_bm = " << par::r_bm << std::endl;
  std::cout << "flavourBreakTolerance = " << par::flavourBreakTolerance << std::endl;
  std::cout << "coalescenceAlgorithm = " << (par::coalescenceAlgorithm == CoalescenceAlgorithm::kClassic ? "kClassic" : "kFromParton") << std::endl;
  std::cout << "inputFile = " << par::inputFile << std::endl;
  std::cout << "outputFile = " << par::outputFile << std::endl;
  std::cout << "obvsFile = " << par::obvsFile << std::endl;
  std::cout << "Configuration loaded successfully." << std::endl;
  std::cout <<"--------------------------" << std::endl;
}