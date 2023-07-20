#ifndef CMDLINE_H
#define CMDLINE_H
#include <string>

// Sequence 1
extern std::string sequence_1; 
// Sequence 2
extern std::string sequence_2;
// The restricted structure for sequence 1
extern std::string structure_1;
// The restricted structure for sequence 2
extern std::string structure_2;

// Types for sequence 1 and 2
extern int type_1;
extern int type_2;

// User set hybrid penalty
extern double hybrid_pen;

// Number of Suboptimal structures to print
extern int subopt;

// The output directory to print to
extern std::string output_dir;

// The output file to print to
extern std::string output_file;

// The input file to take from
extern std::string input_file;

// The file to print the hotspots to
extern std::string hotspot_dir;

// The number of hotspots to make for each sequence
extern int hotspot_num;


/** @brief Where the command line options are stored */
struct args_info
{
  const char *help_help; /**< @brief Print help and exit help description.  */
  const char *version_help; /**< @brief Print version and exit help description.  */
  const char *sequence1_help; /**< @brief Give the first sequence as input  */
  const char *structure1_help; /**< @brief Give a structure corresponding to the first sequence  */
  const char *sequence2_help; /**< @brief Give a second sequence as input */
  const char *structure2_help; /**< @brief Give a structure corresponding to sequence 2*/
  const char *type1_help; /**< @brief Change type for sequence 1 to DNA (base is RNA) */
  const char *type2_help; /**< @brief Change type for sequence 2 to DNA (base is RNA) */
  const char *pen_help; /**< @brief Give a hybrid-penalty to replace the current one */
  const char *subopt_help; /**< @brief Give a number of suboptimals to print  */
  const char *input_help; /**< @brief Give an input file for the results  */
  const char *output_help; /**< @brief Give an output file for the results  */
  const char *dir_help; /**< @brief Give an output directory for the results  */
  const char *h_num_help; /**< @brief Give a number for how many hotspots per sequence */
  const char *h_only_help; /**< @brief Give a file directory for outputting hotspots */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int sequence1_given ;	/**< @brief Whether sequence 1 was given.  */
  unsigned int structure1_given ;	/**< @brief Whether structure 1 was given.  */
  unsigned int sequence2_given ;	/**< @brief Whether sequence 2 was given.  */
  unsigned int structure2_given ;	/**< @brief Whether structure 2 was given.  */
  unsigned int type1_given ;	/**< @brief Whether type 1 was given.  */
  unsigned int type2_given ; /** <@brief whether type 2 was given */
  unsigned int pen_given ;	/**< @brief Whether hybrid penalty was given.  */
  unsigned int subopt_given ;	/**< @brief Whether suboptimals was given.  */
  unsigned int input_given ;	/**< @brief Whether input file was given.  */
  unsigned int output_given ;	/**< @brief Whether output file was given.  */
  unsigned int dir_given ;	/**< @brief Whether output directory was given.  */
  unsigned int h_num_given ;	/**< @brief Whether hotspot_num was given.  */
  unsigned int h_only_given ;	/**< @brief Whether hotspot_only was given.  */

  char **inputs ; /**< @brief unnamed options (options without names) */
  unsigned inputs_num ; /**< @brief unnamed options number */
} ;

/** @brief the purpose string of the program */
extern const char *gengetopt_args_info_purpose;
/** @brief the usage string of the program */
extern const char *gengetopt_args_info_usage;
/** @brief the description string of the program */
extern const char *gengetopt_args_info_description;
/** @brief all the lines making the help output */
extern const char *gengetopt_args_info_help[];


/**
 * Print the help
 */
void cmdline_parser_print_help(void);
/**
 * Print the version
 */
void cmdline_parser_print_version(void);

/**
 * The command line parser
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser (int argc, char **argv,struct args_info *args_info);


/**
 * The command line parser (version with additional parameters - deprecated)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param override whether to override possibly already present options
 * @param initialize whether to initialize the option structure my_args_info
 * @param check_required whether to check that all required options were provided
 * @return 0 if everything went fine, NON 0 if an error took place
 * @deprecated use cmdline_parser_ext() instead
 */
int cmdline_parser2 (int argc, char **argv, struct args_info *args_info);

/**
 * Initializes the passed gengetopt_args_info structure's fields
 * (also set default values for options that have a default)
 * @param args_info the structure to initialize
 */
void cmdline_parser_init (struct args_info *args_info);

/**
 * Deallocates the string fields of the gengetopt_args_info structure
 * (but does not deallocate the structure itself)
 * @param args_info the structure to deallocate
 */
void cmdline_parser_free (struct args_info *args_info);

#endif