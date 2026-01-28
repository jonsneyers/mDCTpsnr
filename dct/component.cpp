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
** $Id: component.cpp,v 1.9 2025/07/08 14:35:22 thor Exp $
**
*/

/// Includes
#include "global/types.hpp"
#include "std/assert.hpp"
#include "std/math.hpp"
#include "dct/line.hpp"
#include "dct/component.hpp"
///

/// Component::Component
Component::Component(void)
  : m_ulY(0), m_iCores(1), m_pOriginal(NULL)
{
  ULONG x,y;
  
  memset(m_pInputLines,0,sizeof(m_pInputLines));
  memset(m_pOutputLines,0,sizeof(m_pOutputLines));

  // Initalize the window function. This is a Hamming window here.
  double xz = double(8 - 1) * 0.5;
  double yz = double(8 - 1) * 0.5;
  double xf = double(8 + 1) * 0.5;
  double yf = double(8 + 1) * 0.5;
  double norm = 0.0;

  for(y = 0;y < 8;y++) {
    for(x = 0;x < 8;x++) {
      double g  = cos((x - xz) / xf * M_PI * 0.5) * cos((y - yz) / yf * M_PI * 0.5);
      //double g  = exp((x - xz) * (x - xz) + (y - yz) * (y - yz));
      norm     += g;
      m_Window[x][y] = g;
    }
  }

  norm /= 8 * 8;

  for(y = 0;y < 8;y++) {
    for(x = 0;x < 8;x++) {
      m_Window[x][y] /= norm;
    }
  }
}
///

/// Component::~Component
Component::~Component(void)
{
  int x,y;

  for(y = 0;y < 8;y++) {
    for(x = 0;x < 8;x++) {
      delete m_pOutputLines[x][y];
      m_pOutputLines[x][y] = NULL;
    }
    delete m_pInputLines[y];
    m_pInputLines[y] = NULL;
  }
  delete m_pOriginal;
}
///

/// Component::PushLine
// Push a new line into the DCT conversion. Note that this is an overcomplete
// representation, i.e. each line triggers a new DCT transformation.
void Component::PushLine(Line &line)
{
  //
  // Must have at least eight pixels, otherwise no output here.
  assert(line.LengthOf() > 8);
  //
  // Check whether there is a target line at the target position. If not, allocate
  // one.
  if (m_pInputLines[m_ulY] == NULL) {
    m_pInputLines[m_ulY] = new Line(line.LengthOf());
  }
  //
  // Swap this into the line buffer, make the output line free again.
  m_pInputLines[m_ulY]->Swap(line);
  //
  // If all lines are full, run a DCT transformation to make the new line
  // available to the analysis.
  if (m_ulY == 7) {
    class Line *top;
    int x,y,w = line.LengthOf() - 8;
    // Ensure the output is available.
    if (m_pOutputLines[0][0] == NULL) {
      for(y = 0;y < 8;y++) {
	for(x = 0;x < 8;x++) {
	  assert(m_pOutputLines[x][y] == NULL);
	  m_pOutputLines[x][y] = new Line(w);
	}
      }
    }
    //
    // Keep the original line as an offset into the input lines, i.e. the
    // center of the windowed block.
    if (m_pOriginal == NULL) {
      m_pOriginal = new Line(w);
    }
    for(x = 0;x < w;x++) {
      m_pOriginal->At(x) = m_pInputLines[4]->Get(x+4);
    }
    //
    // Run the transformation.
    DCT();
    //
    // The topmost line is now no longer required. Rotate down so the next input
    // can go in there.
    top = m_pInputLines[0];
    for(y = 0;y < 7;y++) {
      m_pInputLines[y] = m_pInputLines[y+1];
    }
    m_pInputLines[7] = top;
  } else {
    // Otherwise wait for the next line to be filled.
    m_ulY++;
  }
}
///

/// Component::Run
void Component::Run(int offset,int mod)
{  
  int x,y,i,w; 
  double tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
  double tmp10, tmp11, tmp12, tmp13;
  double z1,z2,z3,z4,z5,z11,z13;
  double in[8][8];

  w = m_pInputLines[0]->LengthOf();

  //
  // Anchor position within the eight input lines.
  for(i = offset;i < w - 8;i+=mod) {
    DOUBLE avg = 0.0;
    // Compute the overall window average.
    for(y = 0;y < 8;y++) {
      for(x = 0;x < 8;x++) {
	double v = m_pInputLines[y]->Get(i+x);
	in[x][y] = v;
	avg += v;
      }
    }
    avg /= 64.0;
    //
    // Fixup the average, window the input.
    for(y = 0;y < 8;y++) {
      for(x = 0;x < 8;x++) {
	in[x][y] = avg * (1 - m_Window[x][y]) + in[x][y] * m_Window[x][y];
      }
    }
    //
    // Now the DCT. First process rows.
    for(y = 0;y < 8;y++) {
      tmp0  = (in[0][y] + in[7][y]);
      tmp7  = (in[0][y] - in[7][y]);
      tmp1  = (in[1][y] + in[6][y]);
      tmp6  = (in[1][y] - in[6][y]);
      tmp2  = (in[2][y] + in[5][y]);
      tmp5  = (in[2][y] - in[5][y]);
      tmp3  = (in[3][y] + in[4][y]);
      tmp4  = (in[3][y] - in[4][y]);
      // Phase 2, even part.
      tmp10 = tmp0 + tmp3;
      tmp13 = tmp0 - tmp3;
      tmp11 = tmp1 + tmp2;
      tmp12 = tmp1 - tmp2;
      // Phase 3
      m_pOutputLines[0][y]->At(i) = (tmp10 + tmp11);
      m_pOutputLines[4][y]->At(i) = (tmp10 - tmp11);
      // Phase 4,5
      z1 = (tmp12 + tmp13) * 0.707106781;
      m_pOutputLines[2][y]->At(i) = (tmp13 + z1);
      m_pOutputLines[6][y]->At(i) = (tmp13 - z1);
      // Phase 2, odd part
      tmp10 = tmp4 + tmp5;
      tmp11 = tmp5 + tmp6;
      tmp12 = tmp6 + tmp7;
      // Rotation
      z5 = (tmp10 - tmp12) * 0.382683433;
      z2 = 0.541196100 * tmp10 + z5;
      z4 = 1.306562965 * tmp12 + z5;
      z3 = tmp11 * 0.707106781;
      // Phase 5
      z11 = tmp7 + z3;
      z13 = tmp7 - z3;
      // Phase 6
      m_pOutputLines[5][y]->At(i) = (z13 + z2);	/* phase 6 */
      m_pOutputLines[3][y]->At(i) = (z13 - z2);
      m_pOutputLines[1][y]->At(i) = (z11 + z4);
      m_pOutputLines[7][y]->At(i) = (z11 - z4);
    }
    //
    // Then process columns.
    for(x = 0;x < 8;x++) {
      tmp0  = m_pOutputLines[x][0]->Get(i) + m_pOutputLines[x][7]->Get(i);
      tmp7  = m_pOutputLines[x][0]->Get(i) - m_pOutputLines[x][7]->Get(i);
      tmp1  = m_pOutputLines[x][1]->Get(i) + m_pOutputLines[x][6]->Get(i);
      tmp6  = m_pOutputLines[x][1]->Get(i) - m_pOutputLines[x][6]->Get(i);
      tmp2  = m_pOutputLines[x][2]->Get(i) + m_pOutputLines[x][5]->Get(i);
      tmp5  = m_pOutputLines[x][2]->Get(i) - m_pOutputLines[x][5]->Get(i);
      tmp3  = m_pOutputLines[x][3]->Get(i) + m_pOutputLines[x][4]->Get(i);
      tmp4  = m_pOutputLines[x][3]->Get(i) - m_pOutputLines[x][4]->Get(i);
      // Phase 2, even part.
      tmp10 = tmp0 + tmp3;
      tmp13 = tmp0 - tmp3;
      tmp11 = tmp1 + tmp2;
      tmp12 = tmp1 - tmp2;
      // Phase 3
      m_pOutputLines[x][0]->At(i) = (tmp10 + tmp11);
      m_pOutputLines[x][4]->At(i) = (tmp10 - tmp11);
      // Phase 4,5
      z1 = (tmp12 + tmp13) * 0.707106781;
      m_pOutputLines[x][2]->At(i) = (tmp13 + z1);
      m_pOutputLines[x][6]->At(i) = (tmp13 - z1);
      // Phase 2, odd part
      tmp10 = tmp4 + tmp5;
      tmp11 = tmp5 + tmp6;
      tmp12 = tmp6 + tmp7;
      // Rotation
      z5 = (tmp10 - tmp12) * 0.382683433;
      z2 = 0.541196100 * tmp10 + z5;
      z4 = 1.306562965 * tmp12 + z5;
      z3 = tmp11 * 0.707106781;
      // Phase 5
      z11 = tmp7 + z3;
      z13 = tmp7 - z3;
      // Phase 6
      m_pOutputLines[x][5]->At(i) = (z13 + z2);	/* phase 6 */
      m_pOutputLines[x][3]->At(i) = (z13 - z2);
      m_pOutputLines[x][1]->At(i) = (z11 + z4);
      m_pOutputLines[x][7]->At(i) = (z11 - z4);
    }
  }
}
///

/// Component::DCT
// Transform the buffered input lines into 64 bands using a lapped
// window DCT.
void Component::DCT(void)
{
  // Run the "real" work now.
  SplitWork(m_iCores);
}
///

