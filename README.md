# HFold_interacting

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
1. [Download the repository](https://github.com/HosnaJabbari/HFold.git) and extract the files onto your system.
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
        HFold_interacting:
            --s1 <sequence1>
            --r1 <restricted_structure1>
            --s2 <sequence2>
            --r2 <restricted_structure2>
            --t1 <type_for_sequence1>
            --t2 <type_for_sequence2>
            -i </path/to/file>
            -o </path/to/file>

        Remarks:
            make sure the <arguments> are enclosed in "", for example --r1 "..().." instead of --r1 ..()..
            input file for -i must be .txt
            if -i is provided with just a file name without a path, it is assuming the file is in the diretory where the executable is called
            if -o is provided with just a file name without a path, the output file will be generated in the diretory where the executable is called
            if -o is provided with just a file name without a path, and if -i is provided, then the output file will be generated in the directory where the input file is located
        
        You can also include SHAPE data to be used. 
        The SHAPE data must be in a file with 1 number per line, starting with the first sequence, 5 lines of 0's, and then the data for the second sequence
        --shape ("filename") to specify a file for shape data
        --b (number) to specify an intercept for the shape data (default is -0.600000)
        --m (number) to specify a slope for the shape data (default is 1.800000)

    
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
    assume you are in the directory where the HFold executable is loacted
    ./HFold_interacting_multimodel -i "/home/username/Desktop/myinputfile.txt" --t1 RNA --t2 DNA
    ./HFold_interacting_multimodel -i "/home/username/Desktop/myinputfile.txt" --t1 RNA --t2 DNA -o "outputfile.txt"
    ./HFold_interacting_multimodel -i "/home/username/Desktop/myinputfile.txt" --t1 RNA --t2 DNA -o "/home/username/Desktop/some_folder/outputfile.txt"
    ./HFold_interacting_multimodel --s1 "GCAACGAUGACAUACAUCGCUAGUCGACGC" --r1 "(____________________________)" --s2 "GCAACGAUGACAUACAUCGCUAGUCGACGCGCAACGAUGACAUACAUCGCUAGUCGACGC" --r2 "(__________________________________________________________)" --t1 RNA --t2 DNA
    ./HFold_interacting_multimodel --s1 "GCAACGAUGACAUACAUCGCUAGUCGACGC" --r1 "(____________________________)" --s2 "GCAACGAUGACAUACAUCGCUAGUCGACGCGCAACGAUGACAUACAUCGCUAGUCGACGC" --r2 "(__________________________________________________________)" --t1 RNA --t2 DNA -o "outputfile.txt"
    ./HFold_interacting_multimodel --s1 "GCAACGAUGACAUACAUCGCUAGUCGACGC" --r1 "(____________________________)" --s2 "GCAACGAUGACAUACAUCGCUAGUCGACGCGCAACGAUGACAUACAUCGCUAGUCGACGC" --r2 "(__________________________________________________________)" --t1 RNA --t2 DNA -o "/home/username/Desktop/some_folder/outputfile.txt"
    ./HFold_interacting_multimodel --s1 "GCAACGAUGACAUACAUCGCUAGUCGACGC" --r1 "(____________________________)" --s2 "GCAACGAUGACAUACAUCGCUAGUCGACGCGCAACGAUGACAUACAUCGCUAGUCGACGC" --r2 "(__________________________________________________________)" --t1 RNA --t2 DNA --shape "shapefile" --b -0.4 --m 1.3

    
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
