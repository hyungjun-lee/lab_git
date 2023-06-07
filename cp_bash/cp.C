
#include <iostream>
#include <fstream>
#include <vector>
#include <TROOT.h>

int cp() {
    std::ifstream inputFile("list.txt"); // Replace with your file name or path
    
    if (!inputFile.is_open()) {
        std::cout << "Error opening file." << std::endl;
        return 1;
    }

    std::vector<char*> lines; // Vector to store the lines

    std::string line;
    while (std::getline(inputFile, line)) {
        std::string prefixedLine = "alien://" + line; // Add the "alien://" prefix to the line
        char* charLine = new char[prefixedLine.length() + 1]; // Allocate memory for the char array
        strcpy(charLine, prefixedLine.c_str()); // Copy the line to the char array
        lines.push_back(charLine); // Store the char array in the vector
    }

    inputFile.close();
    int i = 1;

    for (const auto& line : lines) {
        TGrid::Connect("alien:///");
        gSystem->Exec(Form("mkdir -p ./AOD/%i", i));
        TFile::Cp(line, Form("./AOD/%i/AO2D.root", i));
        i++;
        delete[] line; // Free the allocated memory for each line
    }

    return 0;
}