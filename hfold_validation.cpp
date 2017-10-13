//created by Kevin Chang 19 June 2017
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stack>

#include "externs.h"

#include "hfold_validation.h"


//change all {} [] to () in structure
void replaceBrackets(char* structure){
	for(char* it = structure; *it; ++it) {
    char curr = *it;
    // if open bracket
    if (curr == '{' || curr == '[') {
      *it = '(';
    } else
    // if closed bracket
    if (curr == '}' || curr == ']') {
      *it = ')';
    }
  }
}

//check if structure is valid
//check if length match with sequence
//check if any characters other than ._(){}[]
//check if pk-free
bool validateStructure(char* structure, char* sequence, bool print_err){
	if(strlen(structure) > MAXSLEN){
		if(print_err){
			fprintf(stderr,"Structure length greater than max length\n");
		}
		return false;
	}

	if(strlen(structure) <= 0){
		if(print_err){
        	fprintf(stderr,"length of structure is <= 0. This shouldn't even be able to happen.\n");
		}
		return false;
	}

	//printf("strlen: %d %d\n",strlen(structure), strlen(sequence));
	if(strlen(structure) != strlen(sequence)){
		if(print_err){
			fprintf(stderr,"Length of sequence and corresponding structure must have same length\n");
		}
		return false;
	}

	//check if any characters other than ._(){}[]
	for(char* it = structure; *it; ++it) {
		char curr = *it;
		// Structure must only contain ._(){}[
		if (!(curr == '.' || curr == '_'
		      || curr == '(' || curr == '{' || curr == '['
		      || curr == ')' || curr == '}' || curr == ']'))
		{
			if(print_err){
				fprintf(stderr,"Structure must only contain ._(){}[] \n");
			}
			return false;
		}
	}

	std::string openBracketArray ("({[");
	std::string closeBracketArray (")}]");
	std::stack<char> mystack;

	for(int i=0; i<strlen(structure);i++){
		//printf("d: %c\n",structure[i]);
		if(openBracketArray.find_first_of(structure[i]) != -1){ //is open bracket
			mystack.push(structure[i]);
		}else if(closeBracketArray.find_first_of(structure[i]) != -1){ //is close bracket

        // if stack is empty that means there are more right brackets than left brackets
             if (mystack.empty() == true ) {
				 if(print_err){
                 	fprintf(stderr,"Structure is invalid: more right parentheses than left parentheses\n");
				 }
                 return false;
             }
			if(closeBracketArray.find_first_of(structure[i]) != openBracketArray.find_first_of(mystack.top())){ //if current bracket is not corresponding bracket of what we popped
				if(print_err){
					fprintf(stderr,"Structure bracket types must match and be pseudoknot free\n");
				}
				return false;
			}
			mystack.pop();
		}
	}
	if(mystack.empty() == false){
		fprintf(stderr,"Structure is invalid: more left parentheses than right parentheses\n");
		return false;
	}
	return true;


}

//check if sequence is valid with regular expression
//check length and if any characters other than GCAUT
bool validateSequence(const char* string, bool print_err){
	if(strlen(string) > MAXSLEN){
		return false;
	}

	if(strlen(string) <= 0){
		return false;
	}

  // return false if any characters other than GCAUT
  for(const char* it = string; *it; ++it) {
    const char curr = *it;
    if (!(curr == 'G' || curr == 'C' || curr == 'A' || curr == 'U' || curr == 'T')) {
		if(print_err){
        	fprintf(stderr,"Sequence contains character '%c' that is not G,C,A,U, or T.\n",curr);
		}
        return false;
    }
  }

  return true;
}


//if output_path is not a path -> take base path of input file path (if exist) and concat with output_path
//if output_path is a path -> do nothing
void addPath(char** output_path, char* input_path){
	std::string temp_out_path = *output_path;
	std::string temp_in_path = input_path;
	std::size_t out_path_found = temp_out_path.rfind("/");

	if(out_path_found == std::string::npos){ //if out path does not contain '/'
		std::size_t in_path_found = temp_in_path.rfind("/");
		if(in_path_found != std::string::npos){ //if in path contain '/'
			std::string base_path = temp_in_path.substr(0,in_path_found+1);
			std::cout << base_path << '\n';
			temp_out_path = base_path + temp_out_path;
			//std::cout << temp_out_path << '\n';
			temp_out_path.copy(*output_path, temp_out_path.length(), 0);
		}
	}

}

//assume '.' exist, if not found then do nothings
//change filename.txt to filename_timestamp.txt
void addTimestamp(char** path){
	std::string temp = *path;
	//std::cout << temp << '\n';
	std::size_t found = temp.find_last_of(".");

	time_t rawtime;
	struct tm * timeinfo;
	char buffer [80];
	time (&rawtime);
	timeinfo = localtime (&rawtime);
	strftime (buffer,80,"_%F_%X",timeinfo);
	//june 22 2017 http://www.cplusplus.com/reference/string/string/find_last_of/
	if(found != std::string::npos){
		temp.insert(found,buffer);
		temp.copy(*path,temp.length(),0);
	}
}

/*
//check if input file is in correct format and store value into corresponding variable
//return true on success, false on fail
bool validateInteractingInputFile(char* path, char* seq1, char* struc1, char* seq2, char* struc2){
    //printf("path: %s\n", path);
    FILE* fp;
    char line[MAXSLEN];

    fp = fopen(path, "r");
    if(!fp){
        fprintf("File not found\n");
        return false;
    }
    fscanf(fp,"%s\n%s\n%s\n%s",seq1,struc1,seq2,struc2);
    //printf("%s|%s|%s|%s|\n",seq1,struc1,seq2,struc2);
    fclose(fp);

    if(!(strlen(seq1) < strlen(seq2))){
        printf("Line 1 must be shorter than Line 3\n");
        return false;
    }

    if(!(validateSequence(seq1))){
        printf("Line 1 is invalid\n");
        return false;
    }

    if(!(validateSequence(seq2))){
        printf("Line 3 is invalid\n");
        return false;
    }

    if(!(validateStructure(struc1,seq1))){
        printf("Line 2 is invalid\n");
        return false;
    }

    if(!(validateStructure(struc2,seq2))){
        printf("Line 4 is invalid\n");
        return false;
    }


    return true;
}
*/

bool validateInteractingInputFile2(char* path, char* seq1, char* struc1, char* seq2, char* struc2, bool* sequence1Found, bool* structure1Found, bool* sequence2Found, bool* structure2Found){
	FILE* fp;
    char line[MAXSLEN];

    fp = fopen(path, "r");
    if(!fp){
        printf("File not found\n");
        return false;
    }

	while (fgets(line, sizeof(line), fp)) {
		if(strcmp(line, "\n") == 0){
			break;
		}

		//printf("\ncurrent line: %s\n",line);

		if(*sequence2Found && *structure2Found == false){
			//printf("#2\n");
			strncpy(struc2,line,strlen(line)-1);
			if(validateStructure(struc2,seq2, false)){
				//store struct2
				*structure2Found = true;
			}else if(validateSequence(struc2, false)){
				//store seq2
				strncpy(seq2,line,strlen(line)-1);
				*sequence2Found = true;
			}else{
				//error
				return false;
			}
		}

		if(*sequence1Found && *structure1Found == false && *structure2Found == false){
			//printf("#1\n");
		
			strncpy(struc1,line,strlen(line)-1);
			if(validateStructure(struc1,seq1, false)){
				//store struct1
				*structure1Found = true;
			}else if(validateSequence(struc1, false)){
				//store seq2
				strncpy(seq2,line,strlen(line)-1);
				*sequence2Found = true;
			}else{
				//error
				return false;
			}
		}

		if(*sequence1Found && *sequence2Found == false){
			//printf("store seq2\n");
			//store seq2
			strncpy(seq2,line,strlen(line)-1);
			*sequence2Found = true;
		}

		if(*sequence1Found == false){
			//printf("store seq1\n");
			//store seq1
			strncpy(seq1,line,strlen(line)-1);
			*sequence1Found = true;
		}
    }
	
	if(!(*sequence1Found && *sequence2Found)){
		return false;
	}
/*
	if(sequence1Found){
		printf("seq1: %s|\n",seq1);
	}
	
	if(structure1Found){
		printf("str1: %s|\n",struc1);
	}
	printf("\n");
	if(sequence2Found){
		printf("seq2: %s|\n",seq2);
	}

	if(structure2Found){
		printf("str2: %s|\n",struc2);
	}
	printf("\n");
*/
	return true;
}