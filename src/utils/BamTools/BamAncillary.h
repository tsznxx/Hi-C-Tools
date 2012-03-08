/*****************************************************************************
  bamAncillary.h

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licensed under the GNU General Public License 2.0+ license.
******************************************************************************/
#include "bedUtils.h"
#include "stringUtils.h"
#include "BamAux.h"

namespace BamTools {
    void getBamBlocks(const BamAlignment &bam, const RefVector &refs, 
                  BedVec &blocks, bool includeDeletions = true);
}
