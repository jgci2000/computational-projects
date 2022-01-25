#ifndef READ_NN_H
#define READ_NN_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>

std::vector<std::string> split(const std::string& s, char seperator) {
    std::vector<std::string> output;
    std::string::size_type prev_pos = 0, pos = 0;
    while((pos = s.find(seperator, pos)) != std::string::npos) {
        std::string substring( s.substr(prev_pos, pos-prev_pos) );
        output.push_back(substring);
        prev_pos = ++pos;
    }
    output.push_back(s.substr(prev_pos, pos-prev_pos)); // Last word
    return output;
}

void read_NN_talbe(std::string file_name, int *NN_table) {
    std::ifstream neighbour_tables_file(file_name);
    std::string line;
    if (neighbour_tables_file.is_open()) {
        int i = 0;
        while (std::getline(neighbour_tables_file, line)) {
            std::vector<std::string> a = split(line, ' ');
            for (int idx = 0; idx < a.size(); idx++)
                NN_table[i++] = std::stold(a.at(idx));
        }
        neighbour_tables_file.close();
    }
    else
        std::cout << "Unable to open neighbour table file. Invalid lattice size or lattice type." << std::endl;
}

#endif 
