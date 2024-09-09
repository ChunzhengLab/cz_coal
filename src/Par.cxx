#include <iostream>
#include <string>
#include <stdexcept>
#include "Par.h"
#include <boost/program_options.hpp> // 包含 Boost.Program_options 头文件

namespace po = boost::program_options;

// CoalescenceAlgorithm 输出流重载
std::ostream& operator<<(std::ostream& os, const CoalescenceAlgorithm& algo) {
    switch (algo) {
        case CoalescenceAlgorithm::kClassic:
            os << "kClassic";
            break;
        case CoalescenceAlgorithm::kFromParton:
            os << "kFromParton";
            break;
        default:
            throw std::runtime_error("Unknown CoalescenceAlgorithm");
    }
    return os;
}

// CoalescenceAlgorithm 输入流重载
std::istream& operator>>(std::istream& is, CoalescenceAlgorithm& algo) {
    std::string token;
    is >> token;
    if (token == "kClassic") {
        algo = CoalescenceAlgorithm::kClassic;
    } else if (token == "kFromParton") {
        algo = CoalescenceAlgorithm::kFromParton;
    } else {
        throw std::runtime_error("Invalid CoalescenceAlgorithm value: " + token);
    }
    return is;
}

// EventType 输出流重载
std::ostream& operator<<(std::ostream& os, const EventType& event) {
    switch (event) {
        case EventType::kAMPT:
            os << "kAMPT";
            break;
        case EventType::kRandom:
            os << "kRandom";
            break;
        default:
            throw std::runtime_error("Unknown EventType");
    }
    return os;
}

// EventType 输入流重载
std::istream& operator>>(std::istream& is, EventType& event) {
    std::string token;
    is >> token;
    if (token == "kAMPT") {
        event = EventType::kAMPT;
    } else if (token == "kRandom") {
        event = EventType::kRandom;
    } else {
        throw std::runtime_error("Invalid EventType value: " + token);
    }
    return is;
}

namespace par {
  // 全局变量的定义
  bool isDebug = false;
  bool isCreateAnimation = false;
  bool isWriteCoalQA = false;
  bool isWriteEvents = false;
  bool isCalculateObvs = false;
  bool isRemoveHFQuarks = false;
  bool isEnableMassVarify = false;
  bool isEnableQuarkMoveOn = false;
  bool isBalanceQuarkNumber = false;
  bool isTurnOffRecursion = false;

  bool isResampleQuarkXY = false;
  bool isPiesampleXY = false;
  bool isRandomRotateQuarkXY = false;
  bool isForgetZ = false;

  bool isResampleQuarkPxPy = false;
  bool isRandomRotateQuarkPxPy = false;
  bool isPiesamplePxPy = false;

  bool isRandomCoal = false;

  float r_bm = 1.0;
  float flavourBreakTolerance = 0.0;
  EventType eventType = EventType::kAMPT;
  CoalescenceAlgorithm coalescenceAlgorithm = CoalescenceAlgorithm::kFromParton;

  std::string inputFile = "zpc-1.root";
  std::string outputFile = "output.root";
  std::string obvsFile = "obvs.root";

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> zero_or_one(0, 1);

  //质量表
  std::map<int, float> mass = {
    //quark mass
    {2, 0.00216}, //u
    {-2, 0.00216}, //ubar
    {1, 0.00470}, //d
    {-1, 0.00470}, //dbar
    {3, 0.09350}, //s
    {-3, 0.09350}, //sbar

    //meson mass
    {211, 0.13957}, //π^+
    {-211, 0.13957}, //π^-
    {111, 0.13498}, //π^0

    {321, 0.49368}, //K^+
    {-321, 0.49368}, //K^-
    {311, 0.49761}, //K^0
    {-311, 0.49761}, //anti-K^0

    {221, 0.54786}, //η

    {213, 0.77526}, //ρ^+
    {-213, 0.77526}, //ρ^-
    {113, 0.77526}, //ρ^0

    {223, 0.78266}, //ω

    {333, 1.01946}, //φ

    //baryon mass
    {2212, 0.93827}, // p
    {-2212, 0.93827}, // 反p

    {2112, 0.93957}, // n
    {-2112, 0.93957}, // 反n

    {2224, 1.232}, // Δ^++
    {-2224, 1.232}, // 反Δ^++
    {2214, 1.232}, // Δ^+
    {-2214, 1.232}, // 反Δ^+
    {2114, 1.232}, // Δ^0
    {-2114, 1.232}, // 反Δ^0
    {1114, 1.232}, // Δ^-
    {-1114, 1.232}, // 反Δ^-


    {3122, 1.11568}, // Λ^0
    {-3122, 1.11568}, // 反Λ^0

    {3222, 1.18937}, // Σ^+
    {-3222, 1.18937}, // 反Σ^+
    {3212, 1.19255}, // Σ^0
    {-3212, 1.19255}, // 反Σ^0
    {3112, 1.19745}, // Σ^-
    {-3112, 1.19745}, // 反Σ^-

    {3322, 1.31486}, // Ξ^0
    {-3322, 1.31486}, // 反Ξ^0
    {3312, 1.32131}, // Ξ^-
    {-3312, 1.32131}, // 反Ξ^+

    {3334, 1.67245}, // Ω^-
    {-3334, 1.67245}, // 反Ω^+
  };

  // 实现 printConfig 函数
  void printConfig() {
    std::cout << "--------------------------" << std::endl;
    std::cout << "Configuration: " << std::endl;
    std::cout << "isDebug = " << isDebug << std::endl;
    std::cout << "isCreateAnimation = " << isCreateAnimation << std::endl;
    std::cout << "isWriteCoalQA = " << isWriteCoalQA << std::endl;
    std::cout << "isWriteEvents = " << isWriteEvents << std::endl;
    std::cout << "isCalculateObvs = " << isCalculateObvs << std::endl;
    std::cout << "isRemoveHFQuarks = " << isRemoveHFQuarks << std::endl;
    std::cout << "isEnableMassVarify = " << isEnableMassVarify << std::endl;
    std::cout << "isEnableQuarkMoveOn = " << isEnableQuarkMoveOn << std::endl;
    std::cout << "isBalanceQuarkNumber = " << isBalanceQuarkNumber << std::endl;
    std::cout << "isResampleQuarkXY = " << isResampleQuarkXY << std::endl;
    std::cout << "isRandomRotateQuarkXY = " << isRandomRotateQuarkXY << std::endl;
    std::cout << "isPiesampleXY = " << isPiesampleXY << std::endl;
    std::cout << "isResampleQuarkPxPy = " << isResampleQuarkPxPy << std::endl;
    std::cout << "isRandomRotateQuarkPxPy = " << isRandomRotateQuarkPxPy << std::endl;
    std::cout << "isPiesamplePxPy = " << isPiesamplePxPy << std::endl;
    std::cout << "isTurnOffRecursion = " << isTurnOffRecursion << std::endl;
    std::cout << "isForgetZ = " << isForgetZ << std::endl;
    std::cout << "isRandomCoal = " << isRandomCoal << std::endl;
    std::cout << "flavourBreakTolerance = " << flavourBreakTolerance << std::endl;
    std::cout << "eventType = " << eventType << std::endl;
    std::cout << "coalescenceAlgorithm = " << coalescenceAlgorithm << std::endl;
    std::cout << "r_bm = " << r_bm << std::endl;
    std::cout << "inputFile = " << inputFile << std::endl;
    std::cout << "outputFile = " << outputFile << std::endl;
    std::cout << "obvsFile = " << obvsFile << std::endl;
    std::cout << "--------------------------" << std::endl << std::endl;
  }

  // 实现 initConfig 函数，负责处理命令行和配置文件
  void initConfig(int argc, char** argv) {
    // 定义命令行和配置文件的选项
    po::options_description config("Configuration options");
    config.add_options()
        ("isDebug", po::value<bool>(&isDebug)->default_value(false), "Enable debug mode")
        ("isCreateAnimation", po::value<bool>(&isCreateAnimation)->default_value(false), "Enable local draw")
        ("isWriteCoalQA", po::value<bool>(&isWriteCoalQA)->default_value(false), "Write coalescence QA")
        ("isWriteEvents", po::value<bool>(&isWriteEvents)->default_value(false), "Write events")
        ("isCalculateObvs", po::value<bool>(&isCalculateObvs)->default_value(false), "Calculate observables")
        ("isRemoveHFQuarks", po::value<bool>(&isRemoveHFQuarks)->default_value(true), "Remove HF quarks")
        ("isEnableMassVarify", po::value<bool>(&isEnableMassVarify)->default_value(false), "Enable mass verify")
        ("isEnableQuarkMoveOn", po::value<bool>(&isEnableQuarkMoveOn)->default_value(false), "Enable quark move on")
        ("isBalanceQuarkNumber", po::value<bool>(&isBalanceQuarkNumber)->default_value(false), "Balance quark number")
        ("isResampleQuarkXY", po::value<bool>(&isResampleQuarkXY)->default_value(false), "Resample quark XY")
        ("isPiesampleXY", po::value<bool>(&isPiesampleXY)->default_value(false), "Piesample XY")
        ("isRandomRotateQuarkXY", po::value<bool>(&isRandomRotateQuarkXY)->default_value(false), "Randomly rotate quark XY")
        ("isResampleQuarkPxPy", po::value<bool>(&isResampleQuarkPxPy)->default_value(false), "Resample quark PxPy")
        ("isPiesamplePxPy", po::value<bool>(&isPiesamplePxPy)->default_value(false), "Piesample PxPy")
        ("isRandomRotateQuarkPxPy", po::value<bool>(&isRandomRotateQuarkPxPy)->default_value(false), "Randomly rotate quark PxPy")
        ("isTurnOffRecursion", po::value<bool>(&isTurnOffRecursion)->default_value(false), "Turn off recursion")
        ("isForgetZ", po::value<bool>(&isForgetZ)->default_value(false), "Forget Z")
        ("isRandomCoal", po::value<bool>(&isRandomCoal)->default_value(false), "Random coalescence")
        ("flavourBreakTolerance", po::value<float>(&flavourBreakTolerance)->default_value(0.0), "Set flavour break tolerance")
        ("eventType", po::value<EventType>(&eventType)->default_value(EventType::kAMPT), "Set event type")
        ("coalescenceAlgorithm", po::value<CoalescenceAlgorithm>(&coalescenceAlgorithm)->default_value(CoalescenceAlgorithm::kFromParton), "Set coalescence algorithm")
        ("r_bm", po::value<float>(&r_bm)->default_value(0.1), "Set r_bm value")
        ("inputFile", po::value<std::string>(&inputFile)->default_value("zpc-1.root"), "Input file")
        ("outputFile", po::value<std::string>(&outputFile)->default_value("output.root"), "Output file")
        ("obvsFile", po::value<std::string>(&obvsFile)->default_value("obvs.root"), "Observables file")
        ("configFile", po::value<std::string>(), "configuration file");

    po::variables_map vm;

    // 解析命令行参数
    po::store(po::parse_command_line(argc, argv, config), vm);

    // 如果有配置文件，加载配置文件
    if (vm.count("configFile")) {
        std::ifstream ifs(vm["configFile"].as<std::string>());
        if (ifs) {
            po::store(po::parse_config_file(ifs, config), vm);
        } else {
            std::cerr << "Unable to open config file: " << vm["configFile"].as<std::string>() << std::endl;
            exit(1);
        }
    }

    // 应用命令行参数，覆盖配置文件中的设置
    po::notify(vm);

    // 打印最终的配置信息（可选）
    printConfig();
  }
}