/*****************************************************************************
  bamTowig.cpp
  Last-modified: 26 Jul 2012 11:48:00 AM

  (c) 2012 - Yunfei Wang
  Center for Systems Biology
  Department of Molecular & Cell Biology
  University of Texas at Dallas
  tszn1984@gmail.com

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/

#include <iostream>
#include <string>
using namespace std;

#include "version.h"
#include "BamReader.h"
#include "BamAncillary.h"
#include "BamAux.h"
#include "bedFile.h"
using namespace BamTools;

// is Overlap. Since the bed is sorted, comparison will be easier.
bool isOverlap(refstart,refend,start,end){
	if ( start < refend )
		return true;
	return false;
}

// Convert BAM to Wig
void convertBamToBed( const string &bamfile, const string &bamTag){
	// Open BAM file
	BamReader reader;
	reader.Open(bamFile);

	// get header & reference information
	string header = reader.GetHeaderText();
	RefVector refs = reader.GetReferenceData();

	// rip through the BAM file and convert each mapped entry to BED
	BamAlignment bam;

	// Variables for region list
	unsigned int refstart=refend=0;
	unsigned int start,end;
	unsigned int librarySize=0;
	while (reader.GetNextAlignment(bam)) 
	{
		if (bam.IsMapped() == true)
		{
			librarySize++;
			// strand
			string strand = "+";
			if (bam.IsReverseStrand() == true) strand = "-";

			// start, end
			start= bam.Position;
			end = bam.GetEndPosition(false);

			// extend
			if ( extend ){
				if (strand == "+") 
				{
					end += extend;
					end = end > chromlst[chrom]? chromlst[chrom] : end;
				}
				else
				{
					start -=extend;
					start = start>=0 ? start : 0;
				}
			}

			// Merge overlapped regions
			if ( chrom==lastchrom  && start < refend ){ // overlap with the last bed
				refend = end > refend ? end : refend;
				positions.push_back(start);
				positions.push_back(end);
			}
			else {
				coverages[chrom][refstart] = calculateCoverage(positions);
				positions.clear();
				refstart = start;
				refend = end;
			}
			positions.push_back(start);
			positions.push_back(end);
		}
	}





int main(int argc, char* argv[]){
    return 0;
}

