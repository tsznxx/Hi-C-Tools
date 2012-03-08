/*****************************************************************************
  common.h
  Last-modified: 08 Mar 2012 11:03:18 AM

  (c) 2012 - Yunfei Wang
  Center for Systems Biology
  Department of Molecular & Cell Biology
  University of Texas at Dallas
  tszn1984@gmail.com

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/

#ifndef COMMON_H
#define COMMON_H

#include <stdint.h>

/*************************************************
 Version definition
*************************************************/

#define VERSION "1.0"



/*************************************************
 Bed definitions
*************************************************/

#define MAXCHR 1000000000

typedef uint32_t CHRPOS;
typedef uint16_t BINLEVEL;
typedef uint32_t BIN; 
typedef uint16_t USHORT;
typedef uint32_t UINT;


/*************************************************
 Genome binning constants
*************************************************/

const BIN      _numBins   = 37450;
const BINLEVEL _binLevels = 7; 

// bins range in size from 16kb to 512Mb
// Bin  0          spans 512Mbp,   # Level 1
// Bins 1-8        span 64Mbp,     # Level 2
// Bins 9-72       span 8Mbp,      # Level 3
// Bins 73-584     span 1Mbp       # Level 4
// Bins 585-4680   span 128Kbp     # Level 5
// Bins 4681-37449 span 16Kbp      # Level 6
const BIN _binOffsetsExtended[] = {32678+4096+512+64+8+1, 4096+512+64+8+1, 512+64+8+1, 64+8+1, 8+1, 1, 0};
//const BIN _binOffsetsExtended[] = {4096+512+64+8+1, 4096+512+64+8+1, 512+64+8+1, 64+8+1, 8+1, 1, 0};

const USHORT _binFirstShift = 14;       /* How much to shift to get to finest bin. */
const USHORT _binNextShift  = 3;        /* How much to shift to get to next larger bin. */


/*************************************************
 Parameter check macro
*************************************************/

#define PARAMETER_CHECK(param, paramLen, actualLen) (actualLen == paramLen) && (!strcmp(argv[i], param))


/*************************************************

*************************************************/










#endif //COMMON_H

