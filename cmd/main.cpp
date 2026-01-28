/************************************************************************************
 **  Copyright (C) 2009 University of Stuttgart, Thomas Richter                    **
 **                                                                                **
 **  This software is provided 'as-is', without any express or implied             **
 **  warranty.  In no event will the authors be held liable for any damages        **
 **  arising from the use of this software.                                        **
 **                                                                                **
 **  Permission is granted to anyone to use this software for any purpose,         **
 **  including commercial applications, and to alter it and redistribute it        **
 **  freely, subject to the following restrictions:                                **
 **                                                                                **
 **  1. The origin of this software must not be misrepresented; you must not       **
 **     claim that you wrote the original software. If you use this software       **
 **     in a product, an acknowledgment in the product documentation would be      **
 **     appreciated but is not required.                                           **
 **  2. Altered source versions must be plainly marked as such, and must not be    **
 **     misrepresented as being the original software.                             **
 **  3. This notice may not be removed or altered from any source distribution.    **
 **                                                                                **
 **     Thomas Richter                                                             **
 **     thomas.richter@iis.fraunhofer.de                                           **
 **                                                                                **
 ************************************************************************************/

/*
** $Id: main.cpp,v 1.13 2025/07/09 14:26:45 thor Exp $
**
*/

/// Includes
#include "cmd/main.hpp"
#include "img/image.hpp"
#include "img/imgwriter.hpp"
#include "io/filestream.hpp"
#include "std/stdio.hpp"
#include "std/stdarg.hpp"
#include "std/stdlib.hpp"
#include "global/exceptions.hpp"
#include "std/math.hpp"
#include "measure/pooling.hpp"
///

/// StdExceptionPrinter
// A very simple class that manages the printing of exceptions.
class StdExceptionPrinter : public ExceptionPrinter {
  //
public:
  StdExceptionPrinter(void)
  {}
  ~StdExceptionPrinter(void)
  {}
  //
  // This is the real thing.
  virtual void PrintException(const char *fmt,...);
};

void StdExceptionPrinter::PrintException(const char *fmt,...)
{
  va_list args;

  va_start(args,fmt);
  vfprintf(stderr,fmt,args);
  va_end(args);
}
///

/// A structure keeping the local settings
// The following structure keeps the selected command line
// arguments for the imco project. It is filled in by command
// line parsing and only used locally here.
struct Settings {

  bool Log;
  //
  // The first input file name. This should be a ppm/pgm file. 
  char *m_pcInputName_1;
  

  // The second input file name. This should be a ppm/pgm file. 
  char *m_pcInputName_2;
  
  // Saliency map names.
  char *m_pcSaliency_1;
  char *m_pcSaliency_2;
  
  // The number of CPUs to run in parallel
  int  ncpus;
  
  //
  // Enable logging?
  int  logme;
  //
  // Output in jnd scales?
  bool jnd;
public:
  Settings(void)
    : Log(false),
      m_pcInputName_1(NULL), m_pcInputName_2(NULL),
      m_pcSaliency_1(NULL), m_pcSaliency_2(NULL),
      ncpus(1), logme(0), jnd(false)
  { 
  }
  //
  ~Settings(void)
  {
    delete[] m_pcInputName_1;
    delete[] m_pcInputName_2;
    delete[] m_pcSaliency_1;
    delete[] m_pcSaliency_2;
  }
  //
  // Print the usage rules for this program.
  void Usage(const char *progname);
  //
  // Parse off the arguments from argc,argv.
  void ParseArgs(int argc,char **argv);
};
///

/// Settings::Usage
// Print th command line synopsis for this
// fine program.
void Settings::Usage(const char *progname)
{
  printf("Usage: %s\n"
	 "\t[-CC #]\t: set the number of CPUs to use in parallel\n"
	 "\t[-v]   \t: be verbose and write dct/masking and error image (more v's, more information)\n"
	 "\t[-jnd] \t: output in JND scale rather than mDCTPSNR scale\n"
	 "\t[-s s1.pfm s2.pfm]\t: Define saliency maps of orignal and distorted image\n"
	 "\tinfile1:\t the original file name.\n"
	 "\tinfile2:\t the distorted file name.\n"
	 "%s currently understands .ppm and .pgm files.default:ssim\n",
	 progname,progname);
}
///

/// Settings::ParseArgs
// Parse off program arguments from the command line arguments.
void Settings::ParseArgs(int argc,char **argv)
{
  
  const char *arg;
  const char *progname = argv[0];
  bool failure = false;
  
  //do I need this??
  argc--;
  argv++; // drop the command name.

 while((arg = *argv++)) {
    argc--;
    // Check whether we currently hold an option.
    if (arg[0] == '-') {
      // Yes, an option. o on parsing the input.
      if (!strcmp(arg,"--")) {
	// abort parsing options, go
	// to the file name arguments.
	break;
      } else if (!strcmp(arg,"-v")) {
	logme = 1;
      } else if (!strcmp(arg,"-vv")) {
	logme = 2;
      } else if (!strcmp(arg,"-vvv")) {
	logme = 3;
      } else if (!strcmp(arg,"-jnd")) {
	jnd   = true;
      } else if (!strcmp(arg,"-s")) {
	if (argv[0] == NULL || argv[1] == NULL) {
	  failure = true;
	} else if (m_pcSaliency_1 || m_pcSaliency_2) {
	  failure = true;
	} else {
	  m_pcSaliency_1 = new char[strlen(argv[0]) + 1];
	  m_pcSaliency_2 = new char[strlen(argv[1]) + 1];
	  strcpy(m_pcSaliency_1,argv[0]);
	  strcpy(m_pcSaliency_2,argv[1]);
	  argc -= 2;
	  argv += 2;
	}
      } else if (!strcmp(arg,"-CC")) {
	if (argv[0]) {
	  ncpus = atoi(argv[0]);
	  if (ncpus <= 0)
	    failure = true;
	  argc--;
	  argv++;
	} else
	  failure = true;
      } else {
	failure = true;
	break;
      }
    } else {
      // Put this argument back.
      argc++;
      argv--;
      break;
    }
 }
    
 // Exactly two arguments must be left:
 // the 2 input files' name.
 if (failure || argc != 2) {
   Usage(progname);
   exit(10);
 }
 //
 m_pcInputName_1  = new char[strlen(argv[0]) + 1];
 m_pcInputName_2 = new char[strlen(argv[1]) + 1];
 //
 // copy the names over.
 strcpy(m_pcInputName_1,argv[0]);
 strcpy(m_pcInputName_2,argv[1]);
}
///

/// The main program loop
int main(int argc,char **argv)
{
  struct Settings settings;
  //
  try {
    class Pooling p;
    DOUBLE err;
    //
    // parse off the command line arguments here.
    settings.ParseArgs(argc,argv);
    //
    // Now check whether we encode or decode.
    class FileStream in1,in2;
    class Image img1,img2;
    //
    // Open the input stream for first image's reading.
    in1.OpenForRead(settings.m_pcInputName_1);
    //
    // Open the input stream for second image's reading.
    in2.OpenForRead(settings.m_pcInputName_2);
    //
    // Read the images from the files, and transform them
    // on the way.
    img1.OpenPNM(&in1,settings.ncpus);
    img2.OpenPNM(&in2,settings.ncpus);
    //
    if (img1.ComponentCountOf() != img2.ComponentCountOf()) {
      Throw(InvalidParameter,"main","Component counts differ, cannot compare images.\n");
    }
    //
    // ensure that the image components all have the same
    // dimensions. Can't handle other images in imco otherwise.
    if (!(img1.WidthOf()  == img2.WidthOf() &&
	  img1.HeightOf() == img2.HeightOf())) {
      Throw(InvalidParameter,"main","Image dimensions differ, cannot compare images.\n");
    } 
    //
    // Dimensions must be at least 8.
    if (img1.WidthOf() <= 8+Masking::MaskSize) {
      Throw(InvalidParameter,"main","Image not wide enough.\n");
    }
    if (img1.HeightOf() <= 8+Masking::MaskSize) {
      Throw(InvalidParameter,"main","Image not high enough.\n");
    }
    //
    if (settings.m_pcSaliency_1 && settings.m_pcSaliency_2) {
      DOUBLE s1,s2;
      {
	class Image sal1,sal2;
	class FileStream in3,in4;
	//
	in3.OpenForRead(settings.m_pcSaliency_1);
	in4.OpenForRead(settings.m_pcSaliency_2);
	//
	sal1.OpenPFM(&in3);
	sal2.OpenPFM(&in4);
	//
	s1  = sal1.FindScale();
	s2  = sal2.FindScale();
      }
      //
      //
      class Image sal1,sal2;
      class FileStream in3,in4;
      in3.OpenForRead(settings.m_pcSaliency_1);
      in4.OpenForRead(settings.m_pcSaliency_2);
      //
      sal1.OpenPFM(&in3);
      sal2.OpenPFM(&in4);
      //
      sal1.SetScaling(1.0 / s1);
      sal2.SetScaling(1.0 / s2);
      //
      err = p.Measure(&img1,&img2,&sal1,&sal2,settings.ncpus,settings.logme);
    } else {
      err = p.Measure(&img1,&img2,NULL,NULL,settings.ncpus,settings.logme);
    }

    if (settings.jnd) {
      err = 80.0 - err;
      if (err < 0.0)
	err = 0.0;
      err = 0.000475997802 * pow(err,2.3182);
    }
    printf("%g\n",err);
    //
  } catch(const CodecException &ce) {
    class StdExceptionPrinter ep;
    //
    ce.PrintException(ep);
    return 5;
  }
  return 0;
}
///
