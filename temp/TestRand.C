#include <iostream>
#include <random>

void TestRand() {
    // 创建随机数生成器
    std::random_device rd; // 用于生成种子
    std::mt19937 gen(rd()); // 使用 mt19937 算法的随机数生成器

    // 创建均匀分布对象，范围为 [0, 1]
    std::uniform_int_distribution<int> dis(0, 1);

    // 抽样生成 0 和 1
    for (int i = 0; i < 10; ++i) {
        int sample = dis(gen); // 生成一个随机的 0 或 1
        std::cout << sample << " ";
    }
}