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
** $Id: image.hpp,v 1.9 2025/07/08 14:35:22 thor Exp $
**
*/

#ifndef IMG_IMAGE_HPP
#define IMG_IMAGE_HPP

/// Includes
#include "global/types.hpp"
#include "std/assert.hpp"
#include "ctrafo/colortransformer.hpp"
#include "dct/component.hpp"
///

/// Forwards
class Component;
class ByteStream;
class Line;
///

/// Image
// This simplistic class holds all components and all
// data of an image. This is typically RGB three component
// image data.
class Image {
  //
  // Number of components we have here.
  UWORD             m_usComponents;
  //
  // An array of components.
  class Component **m_ppDCTArray;
  //
  // An array of lines. They are temporaries.
  class Line      **m_ppLineArray;
  //
  // The input bytestream used for loading.
  class ByteStream *m_pIn;
  //
  // The dimension of the image.
  ULONG             m_ulWidth;
  ULONG             m_ulHeight;
  ULONG             m_ulY;
  UBYTE             m_ucBits;
  DOUBLE            m_dScaling;
  //
  // The color transformer.
  class ColorTransformer m_CTrafo;
  //
  // Service (not required elsewhere): Read an ascii string from the input file,
  // encoding a number. This number gets returned. Throws on error.
  LONG ReadNumber(class ByteStream *from);
  // Write an Ascii string to a bytestream.
  void WriteNumber(class ByteStream *to,LONG number);
  // Read an ASCII encoded floating point number
  DOUBLE ReadFloat(class ByteStream *from);
  // Skip blank spaces in the bytestream.
  void SkipBlanks(class ByteStream *from);
  // Skip comment lines starting with #
  void SkipComment(class ByteStream *in);
  //
public:
  // Create an image for the given number of components. It is empty
  // afterwards.
  Image(void);
  //
  // Dispose the image again.
  ~Image(void);
  //
  // Get the number of components in here.
  UWORD ComponentCountOf(void) const
  {
    return m_usComponents;
  }
  //
  // Define the scaling for the PNM probability map.
  void SetScaling(DOUBLE s)
  {
    m_dScaling = s;
  }
  //
  // Load an image from an already open (binary) PPM or PGM file
  // Throw in case the file should be invalid.
  void OpenPNM(class ByteStream *input,int cores);
  //
  // Load a PFM (extended floating point format used by the Itty/Koch saliency toolkit) image.
  void OpenPFM(class ByteStream *input);
  //
  // Read the next line from the image.
  void ReadNextLine(void);
  //
  // Read the next line of saliency data from a PFM file. This line is not DCT transformed
  // but is a probability directly.
  class Line *ReadNextPFMLine(void);
  //
  // Find the scale (maximum number) in the lines of a PFM file.
  DOUBLE FindScale(void);
  //
  // Return the next available DCT output line for the given component and DCT index.
  // Returns NULL if this line is not yet available.
  class Line *GetDCTBand(int h,int v,int component) const
  {
    assert(component < m_usComponents);
  
    return m_ppDCTArray[component]->GetDCTBand(h,v);
  }
  //
  // Return the original line of the given image, offset
  // by the half DCT block size.
  class Line *GetOriginal(int component) const
  { 
    assert(component < m_usComponents);
  
    return m_ppDCTArray[component]->GetOriginal();
  }
  //
  // Returns true if the image is handled completely.
  bool ImageDone(void) const
  {
    assert(m_pIn);

    return (m_ulY >= m_ulHeight);
  }
  //
  // Return the width of the image.
  ULONG WidthOf(void) const
  {
    assert(m_pIn);

    return m_ulWidth;
  }
  //
  // Return the height of the image.
  ULONG HeightOf(void) const
  {
    assert(m_pIn);

    return m_ulHeight;
  }
};
///

///
#endif
