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

#ifndef MEASURE_POOLING_HPP
#define MEASURE_POOLING_HPP

/// Includes
#include <cstdint>
#include "img/imgwriter.hpp"
#include "measure/masking.hpp"
#include "global/thread.hpp"
///

/// Forwards
class Image;
class Line;
///

/// Class Pooling
// This class performs the error pooling and final computation from the
// masking and the line data.
class Pooling : private Thread {
  // 
  //
  class Image    *m_pReference;
  class Image    *m_pDistorted;
  //
  int             m_iComponents;
  int             m_iLogging;
  //
  class ImgWriter m_DCTImgs1[8][8][3];
  class ImgWriter m_DCTImgs2[8][8][3];
  class ImgWriter m_DCTDiff[8][8][3];
  class ImgWriter m_MaskImgs1[8][8][3];
  class ImgWriter m_MaskImgs2[8][8][3];
  class ImgWriter m_SrcImgs1[3];
  class ImgWriter m_SrcImgs2[3];
  class ImgWriter m_ErrorMap;
  class ImgWriter m_LowPasses1[2][2];
  class ImgWriter m_LowPasses2[2][2];
  class Masking   m_Weights1[8][8][3];
  class Masking   m_Weights2[8][8][3];
  //
  void Run(int offset,int mod);
  //
  // Measure the visibility in a band
  // and adjust the error map.
  void MeasureInBand(class Line *refline,class Line *dstlinst,
		     const class Line *refmask,const class Line *dstmask,
		     int w,int logme,float visbase,
		     class Line *errorline);
  //
public:
  Pooling(void);
  ~Pooling(void);
  //
  // Measure the difference between the two images, return the result.
  // The saliency images are probability maps and optional (may be nullptr)
  // Optionally create output log images.
  float Measure(class Image *ref ,class Image *dist,
		 class Image *sref,class Image *dref,
		 int cores,int loglevel);
  //
};
///

///
#endif
