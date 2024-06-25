#include <iostream>
#include <vector>
#include <TFile.h>
#include <TTree.h>

using namespace std;

int ReadTree() {
    // 打开文件
    TFile *file = new TFile("vector_tree.root", "READ");
    TTree *tree = (TTree*)file->Get("tree");

    // 定义存储读取数据的vector
    vector<float> *data = nullptr;
    vector<int> *indices = nullptr;

    // 设置读取tree的分支
    tree->SetBranchAddress("data", &data);
    tree->SetBranchAddress("indices", &indices);

    // 读取数据
    int nEntries = tree->GetEntries();
    for (int i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        
        cout << "Entry " << i << ":" << endl;
        for (size_t j = 0; j < data->size(); ++j) {
            cout << "  data[" << j << "] = " << data->at(j);
            cout << ", indices[" << j << "] = " << indices->at(j) << endl;
        }
    }

    // 关闭文件
    file->Close();
    delete file;

    return 0;
}