# DinoKnot

#### Description:
This software is still under development

#### Supported OS: 
Linux 
macOS 

### Installation:  
Requirements: A compiler that supports C++11 standard (tested with g++ version 4.7.2 or higher)  and CMake version 3.1 or greater.    

[CMake](https://cmake.org/install/) version 3.1 or greater must be installed in a way that HFold can find it.    
To test if your Mac or Linux system already has CMake, you can type into a terminal:      
```
cmake --version
```
If it does not print a cmake version greater than or equal to 3.1, you will have to install CMake depending on your operating system.

#### Mac:    
Easiest way is to install homebrew and use that to install CMake.    
To do so, run the following from a terminal to install homebrew:      
```  
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"   
```    
When that finishes, run the following from a terminal to install CMake.     
```   
brew install cmake   
``` 
#### Linux:    
Run from a terminal     
```
wget http://www.cmake.org/files/v3.8/cmake-3.8.2.tar.gz
tar xzf cmake-3.8.2.tar.gz
cd cmake-3.8.2
./configure
make
make install
```
[Linux instructions source](https://geeksww.com/tutorials/operating_systems/linux/installation/downloading_compiling_and_installing_cmake_on_linux.php)

#### Steps for installation   
1. [Download the repository](https://github.com/HosnaJabbari/DinoKnot/archive/master.zip) and extract the files onto your system.
2. From a command line in the root directory (where this README.md is) run
```
cmake -H. -Bbuild
cmake --build build
```   
If you need to specify a specific compiler, such as g++, you can instead run something like   
```
cmake -H. -Bbuild -DCMAKE_CXX_COMPILER=g++
cmake --build build
```   
This can be useful if you are getting errors about your compiler not having C++11 features.

#### How to use:
    Arguments:
        DinoKnot:
            --s1 <sequence1>
            --r1 <restricted_structure1>
            --s2 <sequence2>
            --r2 <restricted_structure2>
            --t1 <type_for_sequence1>
            --t2 <type_for_sequence2>
            --pen <penalty of base pairs to interact>
            --n <max_number of suboptimal structure> (default is hotspot_num * hotspot_num)
            --o_dir </path/to/folder/to/contain/all/output/files>
            -i </path/to/file>
            -o </path/to/file>
            --hotspot_num <max_number_of_restricted_structures> (default is 20)
            --hotspot_only </path/to/file>
            --method_num <run_method_number> (range from 0 to 4, default is 0 which runs all of them)

        Remarks:
            Required arguments: 
            1. --s1 <sequence1>, --s2 <sequence2>, --t1 <type_for_sequence1>, --t2 <type_for_sequence2>
            or
            2. -i </path/to/file>, --t1 <type_for_sequence1>, --t2 <type_for_sequence2>
            

            make sure the <arguments> are enclosed in "", for example --r1 "..().." instead of --r1 ..()..
            input file for -i must be .txt

            if -i is provided with just a file name without a path, it is assuming the file is in the diretory where the executable is called

            if -o is provided with just a file name without a path, the output file will be generated in the diretory where the executable is called

            if -o is provided with just a file name without a path, and if -i is provided, then the output file will be generated in the directory where the input file is located

            --o_dir would generate a output file for each pair of restricted structure in ascending order of free energy, so output_0.txt would be the file with the minimum free energy that -o would write

            --pen are mainly used when we try to tune the optimal penalty. It is currently set to 22.96551130344778 for RNA-RNA, 34.14979695798525 for DNA-DNA and 166.0 for other types of interaction. Change this at your own risk

            -n, --hotspot_num and --r1 --r2 should not be used together
            
            --method_num forces the program to only run a specfic method. The details of each method can be found in the paper. We do not advise users to change this argument as it is mainly for internal usage. 

            --hotspot_only does not run the actual DinoKnot main program. This tells the program to generate 
            restricted structures and write to the provided file. Usually used for large sequences and need to split it 
            up and run it independently with an outside script
            example file when calling ./DinoKnot_multimodel --s1 "GCAACGAUGACAUACAUCGCUAGUCGACGC" --s2 "GCAACGAUGACAUACAUCGCUAGUCGACGCGCAACGAUGACAUACAUCGCUAGUCGACGC" --t1 RNA --t2 RNA --hotspot_only ./hotspot_file.txt --hotspot_num 7
            Seq1_hotspot_0: ____(((((_____)))))___________
            Seq1_hotspot_1: ____((((______________))))____
            Seq1_hotspot_2: ____((((_______))))___________
            Seq1_hotspot_3: _______((((___________))))____
            Seq1_hotspot_4: ____(((_________)))___________
            -----
            Seq2_hotspot_0: ______________________((((((________________________))))))__
            Seq2_hotspot_1: ______________(((((_______________)))))_____________________
            Seq2_hotspot_2: ______________________(((((__________________________)))))__
            Seq2_hotspot_3: __________________((((((________________________))))))______
            Seq2_hotspot_4: __________________________________(((((_____)))))___________
            Seq2_hotspot_5: ____(((((_____)))))_________________________________________
            Seq2_hotspot_6: _________________(((________)))_____________________________
            Seq2_hotspot_7: __________________(((((__________________________)))))______


    
    Sequence requirements:
        containing only characters GCAUT

    Structure requirements:
        -pseudoknot free
        -containing only characters ._(){}[]
        Remarks:
            Restricted structure symbols:
                () restricted base pair
                _ no restriction
    
    Type options:
        DNA
        RNA
        PMO

    Input file requirements:
            Line1: Sequence1
            Line2: Structure1
            Line3: Sequence2
            Line4: Structure2

        sample:
            GCAACGAUGACAUACAUCGCUAGUCGACGC
            (____________________________)
            GCAACGAUGACAUACAUCGCUAGUCGACGCGCAACGAUGACAUACAUCGCUAGUCGACGC
            (__________________________________________________________)

#### Example:
    assume you are in the directory where the HFold executable is located
    ./DinoKnot_multimodel -i "/home/username/Desktop/myinputfile.txt" --t1 RNA --t2 DNA
    ./DinoKnot_multimodel -i "/home/username/Desktop/myinputfile.txt" --t1 RNA --t2 DNA -o "outputfile.txt"
    ./DinoKnot_multimodel -i "/home/username/Desktop/myinputfile.txt" --t1 RNA --t2 DNA -o "/home/username/Desktop/some_folder/outputfile.txt"
    ./DinoKnot_multimodel --s1 "GCAACGAUGACAUACAUCGCUAGUCGACGC" --r1 "(____________________________)" --s2 "GCAACGAUGACAUACAUCGCUAGUCGACGCGCAACGAUGACAUACAUCGCUAGUCGACGC" --r2 "(__________________________________________________________)" --t1 RNA --t2 DNA
    ./DinoKnot_multimodel --s1 "GCAACGAUGACAUACAUCGCUAGUCGACGC" --s2 "GCAACGAUGACAUACAUCGCUAGUCGACGCGCAACGAUGACAUACAUCGCUAGUCGACGC" --t1 RNA --t2 DNA
    ./DinoKnot_multimodel --s1 "GCAACGAUGACAUACAUCGCUAGUCGACGC" --r1 "(____________________________)" --s2 "GCAACGAUGACAUACAUCGCUAGUCGACGCGCAACGAUGACAUACAUCGCUAGUCGACGC" --r2 "(__________________________________________________________)" --t1 RNA --t2 DNA -o "outputfile.txt"
    ./DinoKnot_multimodel --s1 "GCAACGAUGACAUACAUCGCUAGUCGACGC" --r1 "(____________________________)" --s2 "GCAACGAUGACAUACAUCGCUAGUCGACGCGCAACGAUGACAUACAUCGCUAGUCGACGC" --r2 "(__________________________________________________________)" --t1 RNA --t2 DNA -o "/home/username/Desktop/some_folder/outputfile.txt"
    ./DinoKnot_multimodel --s1 "GCAACGAUGACAUACAUCGCUAGUCGACGC" --s2 "GCAACGAUGACAUACAUCGCUAGUCGACGCGCAACGAUGACAUACAUCGCUAGUCGACGC" --t1 RNA --t2 RNA --hotspot_only ./hotspot_file.txt --hotspot_num 7
    ./DinoKnot_multimodel --s1 "GCAACGAUGACAUACAUCGCUAGUCGACGC" --s2 "GCAACGAUGACAUACAUCGCUAGUCGACGCGCAACGAUGACAUACAUCGCUAGUCGACGC" --t1 RNA --t2 DNA  o_dir "/home/username/Desktop/some_folder"
    ./DinoKnot_multimodel --s1 "GCAACGAUGACAUACAUCGCUAGUCGACGC" --r1 "(____________________________)" --s2 "GCAACGAUGACAUACAUCGCUAGUCGACGCGCAACGAUGACAUACAUCGCUAGUCGACGC" --r2 "(__________________________________________________________)" --t1 RNA --t2 DNA --method_num 2
    
#### Exit code:
    0       success
    1	    invalid argument error 
    3	    thread error
    4       i/o error
    5       pipe error
    6       positive energy error
    10      backtrack error
    11      error when output structure from hfold/hfold_pkonly/simfold does not match restricted
    
    error code with special meaning: http://tldp.org/LDP/abs/html/exitcodes.html
    2	    Misuse of shell builtins (according to Bash documentation)
    126	    Command invoked cannot execute
    127	    "command not found"
    128	    Invalid argument to exit	
    128+n	Fatal error signal "n"
    130	    Script terminated by Control-C
    255	    Exit status out of range (range is 0-255)
