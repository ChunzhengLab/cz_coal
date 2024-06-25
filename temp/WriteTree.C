#include <iostream>
#include <vector>
#include <TFile.h>
#include <TTree.h>

using namespace std;

int WriteTree() {
    // 创建新文件
    TFile *file = new TFile("vector_tree.root", "RECREATE");
    TTree *tree = new TTree("tree", "Tree with std::vector");

    // 定义存储实际数据的vector
    vector<float> data;
    vector<int> indices;

    // 设置tree的分支
    tree->Branch("data", &data);
    tree->Branch("indices", &indices);

    // 填充数据
    for (int i = 0; i < 10; ++i) {
        data.clear();
        indices.clear();
        
        for (int j = 0; j < 5; ++j) {
            data.push_back(i * 2.3 + j);  // 填充一些演示数据
            indices.push_back(i * 5 + j);
        }
        
        tree->Fill();  // 填充tree
    }

    // 写入文件并关闭
    tree->Write();
    file->Close();
    delete file;

    cout << "Data written successfully to 'vector_tree.root'" << endl;
    return 0;
}