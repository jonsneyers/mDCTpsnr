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
** $Id: image.cpp,v 1.11 2025/07/08 14:35:22 thor Exp $
**
*/

/// Includes
#include "options.hpp"
#include "img/image.hpp"
#include "dct/component.hpp"
#include "dct/line.hpp"
#include "ctrafo/colortransformer.hpp"
#include "global/exceptions.hpp"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "io/bytestream.hpp"
///

/// logint2
// Compute the logarithm to the base of two
// as an integer approximation.
static int logint2(unsigned int v)
{
  int c = 0;

  while(v > 1) {
    v >>= 1;
    c++;
  }
  
  return c;
}
///

/// Image::Image
// Create a new image
Image::Image(void)
  : m_usComponents(0), m_ppDCTArray(nullptr), m_ppLineArray(nullptr), m_pIn(nullptr)
{
  m_dScaling = 1.0;
}
///

/// Image::~Image
Image::~Image(void)
{
  if (m_ppDCTArray) { 
    uint16_t i;
    for(i=0;i<m_usComponents;i++) {
      delete m_ppDCTArray[i];
    }
    delete[] m_ppDCTArray;
  }
  if (m_ppLineArray) {
    uint16_t i;
    for(i=0;i<m_usComponents;i++) {
      delete m_ppLineArray[i];
    }
    delete[] m_ppLineArray;
  }
}
///

/// Image::ReadNumber
// Read an ascii string from the input file,
// encoding a number. This number gets returned. Throws on error.
int32_t Image::ReadNumber(class ByteStream *from)
{
  int32_t number   = 0;     // integer number (so far)
  bool negative = false; // sign of the number (true if negative);
  bool valid    = false; // gets true as soon as we get at least one valid digit.
  int32_t in;

  //
  // Skip any leading blanks
  SkipBlanks(from);
  //
  in = from->Get();
  if (in == '+') {
    // number is positive, do nothing
  } else if (in == '-') {
    // number is negative, invert sign.
    negative = true;
  } else {
    // no sign, put data back.
    from->LastUnDo();
  }
  // Skip another block of blanks.
  SkipBlanks(from);
  // Now read digits, one by one.
  do {
    in = from->Get();
    if (in < '0' || in > '9')
      break; // stop on first invalid input.
    // Check whether the number would overflow.
    if (number >= 214748364)
      Throw(OutOfRange,"Image::ReadNumber","input number is too large");
    // Otherwise, add the digit.
    number = number * 10 + in - '0';
    valid  = true;
  } while(true);

  // put back the last digit, we should not have read it as it is
  // not part of the number.
  from->LastUnDo();
  //
  // Now check whether we got at least a single valid digit. If
  // not, there is no valid number to deliver.
  if (!valid)
    Throw(InvalidParameter,"Image::ReadNumber","input number is invalid");
  //
  if (negative)
    number = -number;

  return number;
}
///

/// Image::ReadFloat
float Image::ReadFloat(class ByteStream *from)
{  
  char buffer[256];
  char *end = buffer+255;
  char *bp  = buffer;
  int32_t in;
  float res;
  //
  // Skip any leading blanks
  SkipBlanks(from);
  //
  do {
    in = from->Get();
    if (in == '\n' || in == ' ' || in < 0)
      break;
    
    if (bp > end)
      Throw(OutOfRange,"Image::ReadFloat","floating point number is too long");
    
    *bp++ = in;
  } while(true);

  // put back the last digit, we should not have read it as it is
  // not part of the number.
  from->LastUnDo();
  
  *bp = 0;
  res = static_cast<float>(strtod(buffer,&bp));
  if (*bp != 0)
    Throw(InvalidParameter,"Image::ReadFloat","floating point number is not well formed");

  return res;
}
///

/// Image::WriteNumber
// Write an Ascii string to a bytestream.
void Image::WriteNumber(class ByteStream *to,int32_t number)
{
  char buf[11]; // we need at most 11 digits to represent it.
  char *out = buf;
  //
  sprintf(buf,"%d",int(number));
  //
  while(*out) {
    to->Put(*out);
    out++;
  }
}
///

/// Image::SkipBlanks
// Skip blank spaces in the bytestream.
void Image::SkipBlanks(class ByteStream *from)
{
  int32_t ch;

  do {
    ch = from->Get();
  } while(ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r');
  //
  // Put the last read character back into the stream.
  from->LastUnDo();
}
///

/// Image::SkipComment
// Skip comment lines starting with #
void Image::SkipComment(class ByteStream *in)
{
  int32_t c;
  //
  //
  do {
    c = in->Get();
    // Skip blanks, white spaces, etc.
  } while(c == ' ' || c == '\t' || c == '\r');
  //
  // Check the last character. If a new line, check for comments.
  if (c == '\n') {
    // Possibly a # is following.
    while((c = in->Get()) == '#') {
      // A comment line.
      do {
	c = in->Get();
      } while(c != '\n' && c != ByteStream::Eof && c != '\r');
    }
    // put it back
    in->LastUnDo();
  } else {
    in->LastUnDo();
  }
}
///

/// Image::OpenPNM
// Load an image from an already open (binary) PPM or PGM file
// Throw in case the file should be invalid.
void Image::OpenPNM(class ByteStream *input,int cores)
{
  int32_t data;
  uint16_t i;
  int32_t precision;
  uint32_t width,height;
  uint8_t bits;
  //
  assert(m_ppDCTArray == nullptr && m_pIn == nullptr);
  //
  // Read the header of the file. This must be P6 for a
  // color image, and P5 for a grey-scale image. We currently
  // do not support ASCII images.
  data = input->Get();
  if (data != 'P')
    Throw(InvalidParameter,"Image::LoadPNM","input image stream is no valid PNM file");
  data = input->Get();
  if (data == '6') {
    // A color image. Allocate three components.
    m_usComponents     = 3;
  } else if (data == '5') {
    // A grey scale image. Allocate only one component.
    m_usComponents     = 1;
  } else {
    Throw(InvalidParameter,"Image::LoadPNM","input image is either invalid or an unsupported PNM type");
  }
  // Initialize the component array now.
  m_ppDCTArray = new class Component*[m_usComponents];
  for(i = 0;i<m_usComponents;i++)
    m_ppDCTArray[i] = nullptr;
  //
  // Initialize the line array.
  m_ppLineArray      = new class Line*[m_usComponents];
  for(i = 0;i<m_usComponents;i++)
    m_ppLineArray[i]      = nullptr;
  //
  SkipComment(input);
  // Read the width and the height off the stream. This will also
  // strip blanks we might have in the stream now.
  width     = ReadNumber(input);
  SkipComment(input);
  height    = ReadNumber(input);
  SkipComment(input);
  precision = ReadNumber(input);
  bits      = logint2(precision + 1);
  //
  // Make some consistency checks.
  if (width <= 0 || height <= 0)
    Throw(OutOfRange,"Image::LoadPNM","image dimensions are out of range");
  //
  if (precision <= 0 || precision >= 65536)
    Throw(OutOfRange,"Image::LoadPNM","image precision/bitdepth is out of range");
  //
  // Skip a single whitespace character.
  data = input->Get();
  // Check for MS-Dos line separator \r\n
  if (data == '\r') {
    data = input->Get();
    if (data != '\n') {
      // no MS-Dos separator, MacOs separator! Iek!
      input->LastUnDo();
      data = '\r';
    }
  }
  if (data != ' ' && data != '\n' && data != '\r' && data != '\t')
    Throw(InvalidParameter,"Image::LoadPNM","input image is not a valid PNM file");
  //
  // Now allocate the components.
  for(i=0;i<m_usComponents;i++) {
    m_ppDCTArray[i]       = new class Component();
    m_ppLineArray[i]      = new class Line(width);
    m_ppDCTArray[i]->SetCores(cores);
  }
  //
  // Buffer the settings so we can continue reading.
  m_pIn      = input;
  m_ulWidth  = width;
  m_ulHeight = height;
  m_ulY      = 0;
  m_ucBits   = bits;
}
///

/// Image::OpenPFM
// Load an image from an already open (binary) PFM image
void Image::OpenPFM(class ByteStream *input)
{
  int32_t data;
  uint16_t i;
  float scaling;
  uint32_t width,height;
  //
  assert(m_ppDCTArray == nullptr && m_pIn == nullptr);
  //
  // Read the header of the file. This must be P6 for a
  // color image, and P5 for a grey-scale image. We currently
  // do not support ASCII images.
  data = input->Get();
  if (data != 'P')
    Throw(InvalidParameter,"Image::LoadPFM","input image stream is no valid PFM file");
  data = input->Get();
  if (data == 'F') {
    // Only single precision PFM is understood here. And sufficient.
    m_usComponents     = 1;
  } else {
    Throw(InvalidParameter,"Image::LoadPFM","input image is either invalid or an unsupported PFM type");
  }
  //
  // Initialize the line array.
  m_ppLineArray      = new class Line*[m_usComponents];
  for(i = 0;i<m_usComponents;i++)
    m_ppLineArray[i]      = nullptr;
  //
  SkipComment(input);
  // Read the width and the height off the stream. This will also
  // strip blanks we might have in the stream now.
  width     = ReadNumber(input);
  SkipComment(input);
  height    = ReadNumber(input);
  SkipComment(input);
  scaling   = ReadFloat(input);
  //
  // Make some consistency checks.
  if (width <= 0 || height <= 0)
    Throw(OutOfRange,"Image::LoadPFM","image dimensions are out of range");
  //
  if (scaling <= 0.0f)
    Throw(OutOfRange,"Image::LoadPFM","image scaling is out of range");
  //
  // Skip a single whitespace character.
  data = input->Get();
  // Check for MS-Dos line separator \r\n
  if (data == '\r') {
    data = input->Get();
    if (data != '\n') {
      // no MS-Dos separator, MacOs separator! Iek!
      input->LastUnDo();
      data = '\r';
    }
  }
  if (data != ' ' && data != '\n' && data != '\r' && data != '\t')
    Throw(InvalidParameter,"Image::LoadPFM","input image is not a valid PFM file");
  //
  m_ppLineArray      = new class Line*[m_usComponents];
  for(i = 0;i<m_usComponents;i++)
    m_ppLineArray[i]      = new class Line(width);
  //
  // Buffer the settings so we can continue reading.
  m_pIn      = input;
  m_ulWidth  = width;
  m_ulHeight = height;
  m_ulY      = 0;
  m_ucBits   = 0;
  m_dScaling = scaling;
}
///

/// Image::AllocateAllBuffers
// Eagerly allocate all Component DCT buffers before processing starts.
// This ensures all pointers are valid (never nullptr) after warmup.
// Called after OpenPNM, before warmup loop.
void Image::AllocateAllBuffers(void)
{
  assert(m_ppDCTArray);
  assert(m_ulWidth > 0);
  
  // Allocate buffers for all components
  for (int i = 0; i < m_usComponents; i++) {
    m_ppDCTArray[i]->AllocateBuffers(m_ulWidth);
  }
}
///

/// Image::ReadNextLine
// Read the next line from a PNM image.
void Image::ReadNextLine(void)
{
  assert(m_pIn);
  assert(m_ppDCTArray);
  //
  // Now read the data, component wise interleaved.
  if (m_ulY < m_ulHeight) {
    uint32_t x;
    int i;
    int32_t data;
    int32_t precision = (1UL << m_ucBits) - 1;
    //
    for(x = 0; x < m_ulWidth;x++) {
      for(i=0;i<m_usComponents;i++) {
	data = m_pIn->Get();
	// Fall over if this is an EOF.
	if (data == ByteStream::Eof)
	  Throw(Eof,"Image::LoadPNM","unexpected EOF detected in input image");
	//
	if (m_ucBits > 8) {
	  int32_t dt = m_pIn->Get();
	  if (dt == ByteStream::Eof)
	    Throw(Eof,"Image::LoadPNM","unexpected EOF detected in input image");
	  data = (data << 8) | dt;
	}
	//
	if (data < 0 || data > precision)
	  Throw(OutOfRange,"Image::LoadPNM","the input image contains invalid pixels");
	//
	// Store as integer value (will be converted via lookup table in color transform)
	// This avoids normalization and allows direct lookup table indexing
	m_ppLineArray[i]->At(x) = float(data);
      }
    }
    //
    // Run a color tranformation? Only one or three components supported
    // here.
    if (m_usComponents == 3) {
      m_CTrafo.ForwardsTransform(m_ppLineArray[0],m_ppLineArray[1],m_ppLineArray[2],m_ucBits);
    }
    //
    // Now push into the DCT buffer.
    for(i=0;i<m_usComponents;i++) {
      m_ppDCTArray[i]->PushLine(*m_ppLineArray[i]);
    }
    m_ulY++;
  }
}
///

/// Image::ReadNextPFMLine
// Read the next line of saliency data from a PFM file. This line is not DCT transformed
// but is a probability directly.
class Line *Image::ReadNextPFMLine(void)
{ 
  assert(m_pIn);
  assert(m_ppDCTArray == nullptr);
  assert(m_usComponents == 1);
  assert(sizeof(float) == sizeof(uint32_t));

  if (m_ulY < m_ulHeight) {
    class Line *line = m_ppLineArray[0];
    assert(line);
    for(uint32_t x = 0; x < m_ulWidth;x++) {
      union {
	float f_data;
	uint32_t i_data;
      } u;
      float v;
      int32_t in1,in2,in3,in4;
      // Note that data is safed in little endianness. Yuck!
      in1 = m_pIn->Get();
      in2 = m_pIn->Get();
      in3 = m_pIn->Get();
      in4 = m_pIn->Get();
      //
      if (in4 < 0)
	Throw(Eof,"Image::ReadNextPFMLine","unexpected EOF detected in saliency map");
      //
      u.i_data = (in4 << 24) | (in3 << 16) | (in2 << 8) | (in1 << 0);
      v        = u.f_data * m_dScaling;
      if (v < 0.0f)
	v = 0.0f;
      if (v > 1.0f)
	v = 1.0f;
      line->At(x) = v;
    }
    m_ulY++;
    return line;
  }
  return nullptr;
}
///


/// Image::FindScale
float Image::FindScale(void)
{
  float max = 0.0;

  while (m_ulY < m_ulHeight) {
    class Line *line = ReadNextPFMLine();
    for(uint32_t i = 0;i < line->LengthOf();i++) {
      float v = line->At(i);
      if (v > max)
	max = v;
    }
  }

  return max;
}
///
