/* © 2014-2019, ETH Zurich, Institute of Fluid Dynamics, Stephan Küchlin */ 

/*
 * Options.cpp
 *
 *  Created on: Dec 9, 2014
 *      Author: kustepha
 */

#include <vector>
#include <iostream>
//#include <mpi.h>

// ignore certain warnings in optionparser header
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#include <optionparser.h>
#pragma GCC diagnostic pop

#include "Options.h"

#include "Primitives/MyError.h"
#include "Primitives/MyString.h"

namespace SPH {

namespace Settings {

  enum class opt : unsigned int
  {
    UNKNOWN = 0,
    HELP,
    LICENSES,
    SFILE,
    IDIR,
    ODIR,
    NPC
  };

  static option::ArgStatus Required(const option::Option& option, bool msg) {

    static_cast<void>(msg);

    MYASSERT(option.arg != 0 && !(std::string(option.arg).empty()),
	     std::string("Option \"") + std::string(option.name) + std::string("\" requires an argument\n"));

    return option::ARG_OK;
  }

constexpr option::Descriptor usage[] = {
	//         const unsigned   index
	//             Index of this option's linked list in the array filled in by the parser.
	//         const int   type
	//             Used to distinguish between options with the same index. See index for details.
	//         const char *const   shortopt
	//             Each char in this string will be accepted as a short option character.
	//         const char *const   longopt
	//             The long option name (without the leading -- ).
	//         const CheckArg  check_arg
	//             For each option that matches shortopt or longopt this function will be called to check a potential argument to the option.
	//         const char *    help
	//             The usage text associated with the options in this Descriptor.
	//

	{
		static_cast<unsigned int>(opt::UNKNOWN), 0, "" ,  "",  option::Arg::None,
		"USAGE: EIGER -s \"settings-file\" [options]\n\n"  "Options:"
	},
	{
		static_cast<unsigned int>(opt::HELP),    0, "h",  "help",             option::Arg::None,
		"  --help \t-h  \tPrint usage and exit."
	},
	//          {
	//              static_cast<unsigned int>(opt::LICENSES),0, "l",  "show-licenses",    option::Arg::None,
	//              "  --show-licenses \t-l  \tPrint third-party software licenses and exit."
	//          },
	{
		static_cast<unsigned int>(opt::SFILE),   0, "s",  "settings-file",    Required,
		"  --settings-file \t-s  \tThe file containing JSON formatted settings for the program run."
	},
	{
		static_cast<unsigned int>(opt::IDIR),    0, "i",  "input-directory",  Required,
		"  --input-directory \t-i  \tThe directory relative to which program will read input data."
	},
	{
		static_cast<unsigned int>(opt::ODIR),    0, "o",  "output-directory", Required,
		"  --output-directory \t-o  \tThe directory relative to which program will write output data."
	},
	//          {
	//              static_cast<unsigned int>(opt::NPC),     0, "n",  "prealloc-per-core",Required,
	//              "  --prealloc-per-core \t-n \tThe number of particles per core for which the program should pre-allocate memory."
	//          },
	{
		static_cast<unsigned int>(opt::UNKNOWN), 0, "" ,  "",                 option::Arg::None,
		"\nExamples:\n"
		"  EIGER --settings-file \"/path-to-settings-file/settings-file\"\n"
		"  EIGER -s \"settings-file\" -i \"/data-directory\" -o \"/output-directory\"\n"
	},
	{
	0,0,0,0,0,0
	}
};

option::Option& get ( std::vector<option::Option> & opts , const opt iOpt ) {
	return opts[ static_cast<unsigned int>(iOpt) ];
}

template <typename E>
constexpr unsigned int u(E e) {
    return static_cast<unsigned int>(e);
}

Options::Options () {}


Options::Options (int argc, char* argv[]) {

	argc -= argc > 0;
	argv += argc > 0; // skip program name argv[0] if present

	option::Stats  stats(usage, argc, argv);
	std::vector<option::Option> options(stats.options_max), buffer(stats.buffer_max);
	option::Parser parse(usage, argc, argv, options.data(), buffer.data());

	MYASSERT(!parse.error(),"Fatal error parsing program command line options!");

	if ( options[ u(opt::HELP)] || argc == 0 || options[u(opt::UNKNOWN)] ) {
		option::printUsage(std::cout, usage);
		//MPI_Finalize();
		exit(0);
	}

//	if ( get( options , opt::LICENSES) ) {
//		ThirdpartyLicenses::print_thirdparty_licenses();
//		MPI_Finalize();
//		exit(0);
//	}


	MYASSERT(	options[u(opt::SFILE)] || parse.nonOptionsCount() == 1,
				"Program invoked without specifying settings file.");



	for (option::Option* opt = options[u(opt::UNKNOWN)]; opt; opt = opt->next()){
	   std::cout << "Unknown option: " << opt->name << "\n";
	}


	MYASSERT(!options[u(opt::UNKNOWN)],"Program invoked with unknown options.");


	for (int i = (options[u(opt::SFILE)] ? 0 : 1); i < parse.nonOptionsCount(); ++i) {
	  std::cout << "Non-option #" << i << ": " << parse.nonOption(i) << "\n";
	}

	MYASSERT(	parse.nonOptionsCount() == 0 ||
				(parse.nonOptionsCount() == 1 && !options[u(opt::SFILE)]),
				"Program invoked with unusable non-options.");


	// now get options
	if (options[u(opt::SFILE)]) {
		m_sFileName = std::string(options[u(opt::SFILE)].arg);
	}
	else {
		m_sFileName = std::string(parse.nonOption(0));
	}

	if (options[u(opt::IDIR)]) {
		m_iDir = MyString::makeLastCharDirSep(std::string(options[u(opt::IDIR)].arg));
	}

	if (options[u(opt::ODIR)]) {
		m_oDir = MyString::makeLastCharDirSep(std::string(options[u(opt::ODIR)].arg));
	}

//	if (options[NPC])
//	npc_ = strtoul(options[NPC].arg, 0, 0);
}

Options::~Options() = default;

std::string Options::sFilename () const {
	return m_sFileName;
}

std::string Options::iDir () const {
	return m_iDir;
}

std::string Options::oDir () const {
	return m_oDir;
}

} /* namespace Settings */
} /* namespace SPH */






























