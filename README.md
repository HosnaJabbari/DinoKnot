![image](https://user-images.githubusercontent.com/28453708/112767993-9160b080-8fce-11eb-98b6-906783ecb656.png)
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

Help
========================================

```
Usage: DinoKnot [sequence1] [sequence2] [options]
```

Read RNA and DNA sequences from cmdline; predict minimum\nfree energy and optimum structure

```
     --s1                   Specify the first sequence
     --r1                   Specify the pseuodoknot-free restricted structure for sequence 1 (Will not generate other hotspots)
     --s2                   Specify the second sequence that the first sequence is interacting with
     --r2                   Specify the pseuodoknot-free restricted structure for sequence 2 (Will not generate other hotspots)
     --t1                   Change the type for sequence 1 to DNA (default is RNA)
     --t2                   Change the type for sequence 2 to DNA (default is RNA)
 -p  --pen                  Specify the penalty for the interactions between the sequences
 -n  --opt                  Specify the number of suboptimal structures to output (default is hotspot-num*hotspot-num)
 -i  --input-file           Specify the input file
 -o  --output-file          Specify the path to file to output the results to
 -d  --dir                  Specify the directory for which each results will have a file
     --hotspot-num          Specify the max number of hotspots per sequence (default is 20)
     --hotspot-only         Specify the path to file to output the hotspots to
```
```
Remarks:
    Required arguments: 
    1. --s1 <sequence1>, --s2 <sequence2>

    if -i is provided with just a file name without a path, it is assuming the file is in the diretory where the executable is called

    if -o is provided with just a file name without a path, the output file will be generated in the diretory where the executable is called

    --pen are mainly used when we try to tune the optimal penalty. It is currently set to 22.96551130344778 for RNA-RNA interactions (homo-dimers), 34.14979695798525 for DNA-DNA interactions (homo-dimers) and 166.0 for hybrid interactions (hetero-dimers). Change this at your own risk.

    Sequence requirements:
        containing only characters GCAUT  

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
```

#### Example:
    assume you are in the directory where the DinoKnot executable is located
    ./DinoKnot -i "file.txt" --t2
    ./DinoKnot -i "file.txt" --t2  -o "outputfile.txt"
    ./DinoKnot -i "inputfile.txt" --t2 -o "outputfile.txt"
    ./DinoKnot --s1 "GCAACGAUGACAUACAUCGCUAGUCGACGC" --r1 "(____________________________)" --s2 "GCAACGAUGACAUACAUCGCUAGUCGACGCGCAACGAUGACAUACAUCGCUAGUCGACGC" --r2 "(__________________________________________________________)"
    ./DinoKnot --s1 "GCAACGAUGACAUACAUCGCUAGUCGACGC" --s2 "ACGATTGTGCATCAGCTGA" --t2
    ./DinoKnot --s1 "GCAACGAUGACAUACAUCGCUAGUCGACGC" --r1 "(____________________________)" --s2 "GCAACGAUGACAUACAUCGCUAGUCGACGCGCAACGAUGACAUACAUCGCUAGUCGACGC" --r2 "(__________________________________________________________)" --t1 RNA --t2 DNA -o "outputfile.txt"
    ./DinoKnot --s1 GCAACGAUGACAUACAUCGCUAGUCGACGC --r1 (____________________________) --s2 GCAACGAUGACAUACAUCGCUAGUCGACGCGCAACGAUGACAUACAUCGCUAGUCGACGC --r2 (__________________________________________________________) -o file.txt
    ./DinoKnot --s1 "GCAACGAUGACAUACAUCGCUAGUCGACGC" --s2 "GCAACGAUGACAUACAUCGCUAGUCGACGCGCAACGAUGACAUACAUCGCUAGUCGACGC" --hotspot-only file.txt --hotspot-num 7
    ./DinoKnot --s1 "GCAACGAUGACAUACAUCGCUAGUCGACGC" --s2 "GCAACGAUGACAUACAUCGCUAGUCGACGCGCAACGAUGACAUACAUCGCUAGUCGACGC" --dir /home/username/Desktop/some_folder

```
./DinoKnot --s1 "GCAACGAUGACAUACAUCGCUAGUCGACGC" --r1 "(____________________________)" --s2 "GCAACGAUGACAUACAUCGCUAGUCGACGCGCAACGAUGACAUACAUCGCUAGUCGACGC" --r2 "(__________________________________________________________)"

Seq:          GCAACGAUGACAUACAUCGCUAGUCGACGCXXXXXGCAACGAUGACAUACAUCGCUAGUCGACGCGCAACGAUGACAUACAUCGCUAGUCGACGC
Restricted_0: (____________________________).....(__________________________________________________________)
Result_0:     ....(((((.....(((((...(((((([[.....]]..[[[[[.....]]]]]...))))))......))))).....)))))........... (-16.6)



./DinoKnot --s1 "GCAACGAUGACAUACAUCGCUAGUCGACGC" --s2 "GCAACGAUGACAUACAUCGCUAGUCGACGCGCAACGAUGACAUACAUCGCUAGUCGACGC" --hotspot-only file.txt --hotspot-num 7

Seq1_hotspot_0: ____(((((_____)))))___________(-2.32)
Seq1_hotspot_1: ____((((______________))))____(-1.4)
Seq1_hotspot_2: ____((((_______))))___________(-1.11)
Seq1_hotspot_3: _______((((___________))))____(-0.56)
Seq1_hotspot_4: ____(((_________)))___________(-0.32)
---------------
Seq2_hotspot_0: ______________________((((((________________________))))))__(-5.1)
Seq2_hotspot_1: ______________(((((_______________)))))_____________________(-2.65)
Seq2_hotspot_2: ______________________(((((__________________________)))))__(-2.59)
Seq2_hotspot_3: __________________((((((________________________))))))______(-2.46)
Seq2_hotspot_4: __________________________________(((((_____)))))___________(-2.32)
Seq2_hotspot_5: ____(((((_____)))))_________________________________________(-2.32)
Seq2_hotspot_6: _________________(((________)))_____________________________(-2)

```
