#include "shape_data.h"

#include <fstream>
#include <iostream>

shape_info shape = shape_info();

// prepares filename to be stored in shape_file
void shape_info::set_shape_file(std::string filename) {
    // cut off first 7 characters ("-shape=")
    filename = filename.substr(7,filename.length()-6);

    auto front = filename.front();
    auto back = filename.back();

    // cut quote from beginning and end if they exist
    if (front == '\"')
        filename.erase(front);
    if (back == '\"')
        filename.erase(back);

    shape_file_ = filename;
    use_shape_data_ = true;

    set_data(shape_file_);
}

void shape_info::set_data(std::string filename) {
    // open the file
    std::ifstream infile;
    infile.open(filename, std::ifstream::in);

    // first line may be the sequence
    std::string input;
    char first_char;

    // continue until there is a line that is a number
    while (infile >> input && !is_number(input)) {}

    if (!is_number(input)) {
        printf("SHAPE data file error: file has no numbers in it\n");
        exit(-1);
    }

    // the first number is in input
    // it starts from 1
    data_.push_back(stof(input));
    // go through where each word is a shape data number
    int i = 1;
    while (infile >> input) {
        if (i > sequence_length()) {
            printf("SHAPE data file error: length greater than sequence length (sequence length = %d)\n",sequence_length());
            exit(-1);
        }
        if (!is_number(input)) {
            printf("SHAPE data file error: line after start is not a number\n",i,sequence_length());
            exit(-1);
        }

        data_.push_back(stof(input));
        ++i;
    }

    if (i != sequence_length()) {
        printf("SHAPE data file error: length less than sequence length (%d compared to %d)\n",i-1,sequence_length());
        exit(-1);
    }

    infile.close();
}
