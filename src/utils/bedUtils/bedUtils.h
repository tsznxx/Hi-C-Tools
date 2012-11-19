/*****************************************************************************
  bedFile.h
  Last-modified: 14 Mar 2012 04:01:16 PM

  (c) 2012 - Yunfei Wang
  System biology center
  Department of Cellular and Molecular Biology
  University of Texas at Dallas
  tszn1984@gmail.com

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/

#ifndef BEDUTILS_H
#define BEDUTILS_H

// standard includes
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <iomanip>


// "local" includes
#include "common.h"
#include "gzstream.h"
#include "fileUtils.h"
#include "seqUtils.h"

using namespace std;


//*************************************************
// Common data structures
//*************************************************

struct DEPTH {
    CHRPOS starts;
    CHRPOS ends;
};


/*
    Structure for regular Bed records
*/
class Bed 
{
public:
    // Regular Bed fields
    string chrom;
    CHRPOS start;
    CHRPOS end;
    string name;
    float  score;
    string strand;

    // Add'l fields for Bed12 and/or custom Bed annotations
    vector<string> otherFields;
    string desc;

public:
    // constructors

    // Bed6
    Bed(const string &chrom="", const CHRPOS &start=0, const CHRPOS &end=0, const string &name="", const float &score=0.0, const string &strand=".", const string &desc="");

    // Bed more
    Bed( const string &chrom, const CHRPOS &start, const CHRPOS &end, const string &name, const float &score, const string &strand, const vector<string> &elems, const string &desc="");

    // Bed from Vector
    Bed(const vector<string> &elems, const string &desc="");
	
    // Destructor
    ~Bed(void);

/*
	// Set Bed6
    void update(const string &chrom, const CHRPOS &start, const CHRPOS &end, const string &name="", const float &score=0.0, const string &strand=".", const string &desc="");

	// Set Bed more
	void update(const string &chrom, const CHRPOS &start, const CHRPOS &end, const string &name, const float &score, const string &strand, const vector<string> &elems, const string &desc="");

	// Set Bed from Vector
	void update(const vector<string> &elems, const string &desc="");
*/
	
	// Size
    CHRPOS size(void) const { return end-start;};

    // IsOverlap
    bool isOverlap(const Bed &B, bool forcestrand=false) const
    {
		if (forcestrand && ( (this->strand=="+" && B.strand=="-")|| (this->strand=="-" && B.strand=="+")) )
			return false;
        if (B.chrom!=chrom)
            return false;
        if (end<B.start)
            return false;
        if (start>B.end)
            return false;
		return true;
    }

    // Overlap length
    CHRPOS overlapLength(const Bed &B, bool forcestrand=false) const
    {
        if (isOverlap(B, forcestrand))
            return size() + B.size() -  max(end,B.end) + min(start,B.start);
        else
            return 0;
    }

    // Distance between two Beds
    CHRPOS distance(const Bed &B) const
    {
        if (B.chrom!=chrom)
            return MAXCHR;
        if (end < B.start)
            return B.start-end;
        if (start >B.end)
            return start-B.end;
        return 0;
    }

    // Operators overloaded. "<" is enough for sort.
    bool operator < (const Bed &B)  const
    {
        if (chrom!=B.chrom)
            return (chrom<B.chrom ? true: false);
        if (start!=B.start)
            return (start<B.start ? true: false);
        if (end!=B.end)
            return (end<B.end ? true:false);
        if (score<B.score)
            return true;
        return false;
    }

    //bool operator > (const Bed &B)  const  { return chrom>B.chrom || start > B.start || end >B.end || score >B.score;};
    //bool operator == (const Bed &B) const  { return chrom==B.chrom || start == B.start || end ==B.end; };
    
    // IOstream  <<
    friend inline ostream & operator << (ostream & os, Bed &B )
    {
        os << B.chrom << "\t" << B.start << "\t" << B.end << "\t" << B.name << "\t" << setprecision(4) << B.score << "\t" << B.strand;
        return os;
    }

    // toString
    string toString(bool bedonly=true)
    {
        stringstream ss;
        ss << *this;
        if(!bedonly)
        {
            vector<string>::iterator it;
            for(it = otherFields.begin();it != otherFields.end(); it++)
                ss << "\t" << *it;
        }
        return ss.str();
    }


    // get BIN
    BIN getBIN() const
    {
        BINLEVEL i;
		CHRPOS tstart,tend;
        tstart= start    >> _binFirstShift;
        tend=  (end-1) >> _binFirstShift;
        for( i=0;i<_binLevels;i++)
        {
            if(tstart==tend)
				break;
             tstart>>=_binNextShift;
             tend>>=_binNextShift;
        }
		return _binOffsetsExtended[i] + tstart;
    }

    // get Seq
    string getSeq(string twobitFile) const
    {
        string seq;
        seq.reserve(end-start);
        
        // get the sequence from TwoBit file.
		//
		// Add code here!!!!
		//
        return seq;
    }

    //More methods here ... 
}; // Bed
 
// enum to flag the state of a given line in a Bed file.
enum BedLineStatus
{
    Bed_INVALID = -1,
    Bed_HEADER  = 0,
    Bed_BLANK   = 1,
    Bed_VALID   = 2
};

// enum to indicate the type of file we are dealing with
enum FileType
{
    BED_FILETYPE,
    TAB_FILETYPE,
    GFF_FILETYPE,
    VCF_FILETYPE,
	BOWTIE_FILETYPE
};

//*************************************************
// Data structure typedefs
//*************************************************

typedef vector<Bed>    BedVec;
typedef vector<Bed*>   pBedVec;
//typedef vector<MATE> mateVector;

typedef map<BIN, BedVec,    std::less<BIN> > binBeds;
//typedef map<BIN, mateVector, std::less<BIN> > binsToMates;

typedef map<string, binBeds, std::less<string> >    BedMap;
//typedef map<string, binsToMates, std::less<string> > masterMateMap;
typedef map<string, BedVec, std::less<string> >     BedMapNoBin;

// Hits of inersectBed < CHRPOS overlap, pBedVec overlappedBeds>
typedef map<CHRPOS,pBedVec, std::greater<CHRPOS> > OverlapMap;

/************************************************
 * BedUtils
************************************************/

namespace BedUtils
{
	// Transform column line to Bed format.
	void colToBed(Bed &tbed, const vector<string> &elems);
	void bowtieToBed(Bed &tbed, const vector<string> &elems);
	void tabToBed(Bed &tbed, const vector<string> &elems);

	// BedVec
	void loadBedFileToVec(BedVec &bedvec, const string &bedfile);

	// BedMap related
	void loadBedFileToMap   (BedMap &bedmap, const string &bedfile);
	void loadBedVecToMap (BedMap &bedmap, const BedVec & bedvector);

	// IntersectBed
	void intersectBed (const Bed &tbed, BedMap &bedmap, OverlapMap &overlaps, const bool &forcestrand = false);
	
};

//public:
//
//    // Constructor
//    BedFile(string &);
//
//    // Destructor
//    ~BedFile(void);
//
//    /********* File management ********/
//    // Open a Bed file for reading (creates an istream pointer)
//    void Open(void);
//
//    // Close an opened Bed file.
//    void Close(void);
//
//    // are the any intervals left in the file?
//    bool Empty(void);
//
//    // Rewind the pointer back to the beginning of the file
//    void Rewind(void);
//
//    // Jump to a specific byte in the file
//    void Seek(unsigned long offset);
//
//    // dump the header, which is collected as part of Open()
//    void PrintHeader(void);
//
//    // get the next line in the file. splits a line in _bedFields
//    void GetLine(void);
//
//    // Get the next Bed entry in an opened Bed file.
//    bool getNextBed (Bed &bed, bool forceSorted = false);
//    
//    // Returns the next MERGED (i.e., non-overlapping) interval in an opened Bed file
//    // NOTE: assumes input file is sorted by chrom then start
//    bool getNextMergedBed(Bed &merged_bed);
//
//    // load a Bed file into a map keyed by chrom, then bin. value is vector of Beds
//    void loadBedFileIntoMap();
//
//    // load a Bed file into a map keyed by chrom, then bin. value is vector of BedCOVs
//    void loadBedCovFileIntoMap();
//
//    // load a Bed file into a map keyed by chrom, then bin. value is vector of BedCOVLISTs
//    void loadBedCovListFileIntoMap();
//
//    // load a Bed file into a map keyed by chrom. value is vector of Beds
//    void loadBedFileIntoMapNoBin();
//
//    // Given a chrom, start, end and strand for a single feature,
//    // search for all overlapping features in another Bed file.
//    // Searches through each relevant genome bin on the same chromosome
//    // as the single feature. Note: Adapted from kent source "binKeeperFind"
//    void FindOverlapsPerBin(string chrom, CHRPOS start, CHRPOS end, string strand, vector<Bed> &hits, bool sameStrand, bool diffStrand);
//
//    // return true if at least one overlap was found.  otherwise, return false.
//    bool FindOneOrMoreOverlapsPerBin(string chrom, CHRPOS start, CHRPOS end, string strand,
//                                        bool sameStrand, bool diffStrand, float overlapFraction = 0.0);
//
//    // return true if at least one __reciprocal__ overlap was found.  otherwise, return false.
//    bool FindOneOrMoreReciprocalOverlapsPerBin(string chrom, CHRPOS start, CHRPOS end, string strand,
//                                                    bool sameStrand, bool diffStrand, float overlapFraction = 0.0);
//
//    // Given a chrom, start, end and strand for a single feature,
//    // increment a the number of hits for each feature in B file
//    // that the feature overlaps
//    void countHits(const Bed &a, bool sameStrand = false, bool diffStrand = false, bool countsOnly = false);
//
//    // same as above, but has special logic that processes a set of
//    // Bed "blocks" from a single entry so as to avoid over-counting
//    // each "block" of a single BAM/Bed12 as distinct coverage.  That is,
//    // if one read has four block, we only want to count the coverage as
//    // coming from one read, not four.
//    void countSplitHits(const vector<Bed> &bedBlock, bool sameStrand = false, bool diffStrand = false, bool countsOnly = false);
//
//    // Given a chrom, start, end and strand for a single feature,
//    // increment a the number of hits for each feature in B file
//    // that the feature overlaps
//    void countListHits(const Bed &a, int index, bool sameStrand, bool diffStrand);
//
//    // the bedfile with which this instance is associated
//    string bedFile;
//    unsigned int bedType;  // 3-6, 12 for Bed
//                           // 9 for GFF
//    bool isZeroBased;
//
//    // Main data structires used by BedTools
//    masterBedCovMap      bedCovMap;
//    masterBedCovListMap  bedCovListMap;
//    masterBedMap         bedMap;
//    masterBedMapNoBin    bedMapNoBin;
//    
//    BedLineStatus _status;
//    int _lineNum;
//private:
//
//    // data
//    bool _isGff;
//    bool _isVcf;
//    bool _typeIsKnown;        // do we know the type?   (i.e., Bed, GFF, VCF)
//    FileType   _fileType;     // what is the file type? (Bed? GFF? VCF?)
//    istream   *_bedStream;
//    string _bedLine;
//
//    Bed _nullBed;
//    string _header;
//    bool _firstLine;
//    vector<string> _bedFields;
//    int _merged_start;
//    int _merged_end;
//    string _merged_chrom;
//    int _prev_start;
//    string _prev_chrom;
//
//    void setZeroBased(bool zeroBased);
//    void setGff (bool isGff);
//    void setVcf (bool isVcf);
//    void setFileType (FileType type);
//    void updateType (int colNums);
//
//    /************ Private utilities ***********************/
//    void GetHeader(void);
//
//    /******************************************************
//    Private definitions to circumvent linker issues with
//    templated member functions.
//    *******************************************************/
//
//    /*
//        parseLine: converts a lineVector into either Bed or BedCOV (templated, hence in header to avoid linker issues.)
//    */
//    template <typename T>
//    inline BedLineStatus parseLine (T &bed, const vector<string> &lineVector) {
//        
//        // clear out the data from the last line.
//        bed = _nullBed;
//        unsigned int numFields = lineVector.size();
//
//        // bail out if we have a blank line
//        if (numFields == 0) {
//            return Bed_BLANK;
//        }
//        // bail out if we have a comment line
//        if ( (lineVector[0].find("#")       == 0) ||
//             (lineVector[0].find("browser") == 0) ||
//             (lineVector[0].find("track")   == 0) 
//           )
//        {
//            return Bed_HEADER;
//        }
//
//        if (numFields >= 3) {
//            // line parsing for all lines after the first non-header line
//            if (_typeIsKnown == true) {
//                switch(_fileType) {
//                    case Bed_FILETYPE:
//                        if (parseBedLine(bed, lineVector, _lineNum, numFields) == true) return Bed_VALID;
//                    case VCF_FILETYPE:
//                        if (parseVcfLine(bed, lineVector, _lineNum, numFields) == true) return Bed_VALID;
//                    case GFF_FILETYPE:
//                        if (parseGffLine(bed, lineVector, _lineNum, numFields) == true) return Bed_VALID;
//                    default:
//                        printf("ERROR: file type encountered. Exiting\n");
//                        exit(1);
//                }
//            }
//            // line parsing for first non-header line: figure out file contents
//            else {
//                // it's Bed format if columns 2 and 3 are integers
//                if (isInteger(lineVector[1]) && isInteger(lineVector[2])) {
//                    setGff(false);
//                    setZeroBased(true);
//                    setFileType(Bed_FILETYPE);
//                    updateType(numFields);       // we now expect numFields columns in each line
//                    if (parseBedLine(bed, lineVector, _lineNum, numFields) == true) return Bed_VALID;
//                }
//                // it's VCF, assuming the second column is numeric and there are at least 8 fields.
//                else if (isInteger(lineVector[1]) && numFields >= 8) {
//                    setGff(false);
//                    setVcf(true);
//                    setZeroBased(false);
//                    setFileType(VCF_FILETYPE);
//                    updateType(numFields);       // we now expect numFields columns in each line
//                    if (parseVcfLine(bed, lineVector, _lineNum, numFields) == true) return Bed_VALID;
//                }
//                // it's GFF, assuming columns columns 4 and 5 are numeric and we have 9 fields total.
//                else if ((numFields >= 8) && isInteger(lineVector[3]) && isInteger(lineVector[4])) {
//                    setGff(true);
//                    setZeroBased(false);
//                    setFileType(GFF_FILETYPE);
//                    updateType(numFields);       // we now expect numFields columns in each line
//                    if (parseGffLine(bed, lineVector, _lineNum, numFields) == true) return Bed_VALID;
//                }
//                else {
//                    cerr << "Unexpected file format.  Please use tab-delimited Bed, GFF, or VCF. " <<
//                            "Perhaps you have non-integer starts or ends at line " << _lineNum << "?" << endl;
//                    exit(1);
//                }
//            }
//        }
//        else {
//            cerr << "It looks as though you have less than 3 columns at line: " << _lineNum << ".  Are you sure your files are tab-delimited?" << endl;
//            exit(1);
//        }
//        // default
//        return Bed_INVALID;
//    }
//
//
//    /*
//        parseBedLine: converts a lineVector into either Bed or BedCOV (templated, hence in header to avoid linker issues.)
//    */
//    template <typename T>
//    inline bool parseBedLine (T &bed, const vector<string> &lineVector, int _lineNum, unsigned int numFields) {
//
//        // process as long as the number of fields in this
//        // line matches what we expect for this file.
//        if (numFields == this->bedType) {
//            bed.chrom = lineVector[0];
//            int i;
//            i = atoi(lineVector[1].c_str());
//            if (i<0) {
//                 cerr << "Error: malformed Bed entry at line " << _lineNum << ". Start Coordinate detected that is < 0. Exiting." << endl;
//                 exit(1);
//            }
//            bed.start = (CHRPOS)i;
//            i = atoi(lineVector[2].c_str());
//            if (i<0) {
//                cerr << "Error: malformed Bed entry at line " << _lineNum << ". End Coordinate detected that is < 0. Exiting." << endl;
//                exit(1);
//            }
//            bed.end = (CHRPOS)i;
//            
//            // handle starts == end (e.g., insertions in reference genome)
//            if (bed.start == bed.end) {
//                bed.start--;
//                bed.end++;
//                bed.zeroLength = true;
//            }
//            
//            if (this->bedType == 4) {
//                bed.name = lineVector[3];
//            }
//            else if (this->bedType == 5) {
//                bed.name = lineVector[3];
//                bed.score = lineVector[4];
//            }
//            else if (this->bedType == 6) {
//                bed.name = lineVector[3];
//                bed.score = lineVector[4];
//                bed.strand = lineVector[5];
//            }
//            else if (this->bedType > 6) {
//                bed.name = lineVector[3];
//                bed.score = lineVector[4];
//                bed.strand = lineVector[5];
//                for (unsigned int i = 6; i < lineVector.size(); ++i) {
//                    bed.otherFields.push_back(lineVector[i]);
//                }
//            }
//            else if (this->bedType != 3) {
//                cerr << "Error: unexpected number of fields at line: " << _lineNum
//                     << ".  Verify that your files are TAB-delimited.  Exiting..." << endl;
//                exit(1);
//            }
//
//            // sanity checks.
//            if (bed.start <= bed.end) {
//                return true;
//            }
//            else {
//                cerr << "Error: malformed Bed entry at line " << _lineNum << ". Start was greater than end. Exiting." << endl;
//                exit(1);
//            }
//        }
//        else if (numFields == 1) {
//            cerr << "Only one Bed field detected: " << _lineNum << ".  Verify that your files are TAB-delimited.  Exiting..." << endl;
//            exit(1);
//        }
//        else if ((numFields != this->bedType) && (numFields != 0)) {
//            cerr << "Differing number of Bed fields encountered at line: " << _lineNum << ".  Exiting..." << endl;
//            exit(1);
//        }
//        else if ((numFields < 3) && (numFields != 0)) {
//            cerr << "TAB delimited Bed file with at least 3 fields (chrom, start, end) is required at line: "<< _lineNum << ".  Exiting..." << endl;
//            exit(1);
//        }
//        return false;
//    }
//
//
//    /*
//        parseVcfLine: converts a lineVector into either Bed or BedCOV (templated, hence in header to avoid linker issues.)
//    */
//    template <typename T>
//    inline bool parseVcfLine (T &bed, const vector<string> &lineVector, int _lineNum, unsigned int numFields) {
//        if (numFields == this->bedType) {
//            bed.chrom  = lineVector[0];
//            bed.start  = atoi(lineVector[1].c_str()) - 1;  // VCF is one-based
//            bed.end    = bed.start + lineVector[3].size(); // VCF 4.0 stores the size of the affected REF allele.
//            bed.strand = "+";
//            // construct the name from the ref and alt alleles.
//            // if it's an annotated variant, add the rsId as well.
//            bed.name   = lineVector[3] + "/" + lineVector[4];
//            if (lineVector[2] != ".") {
//                bed.name += "_" + lineVector[2];
//            }
//
//            if (this->bedType > 2) {
//                for (unsigned int i = 2; i < numFields; ++i)
//                    bed.otherFields.push_back(lineVector[i]);
//            }
//
//            if ((bed.start <= bed.end) && (bed.start >= 0) && (bed.end >= 0)) {
//                return true;
//            }
//            else if (bed.start > bed.end) {
//                cerr << "Error: malformed VCF entry at line " << _lineNum << ". Start was greater than end. Exiting." << endl;
//                exit(1);
//            }
//            else if ( (bed.start < 0) || (bed.end < 0) ) {
//                cerr << "Error: malformed VCF entry at line " << _lineNum << ". Coordinate detected that is < 0. Exiting." << endl;
//                exit(1);
//            }
//        }
//        else if (numFields == 1) {
//            cerr << "Only one VCF field detected: " << _lineNum << ".  Verify that your files are TAB-delimited.  Exiting..." << endl;
//            exit(1);
//        }
//        else if ((numFields != this->bedType) && (numFields != 0)) {
//            cerr << "Differing number of VCF fields encountered at line: " << _lineNum << ".  Exiting..." << endl;
//            exit(1);
//        }
//        else if ((numFields < 2) && (numFields != 0)) {
//            cerr << "TAB delimited VCF file with at least 2 fields (chrom, pos) is required at line: "<< _lineNum << ".  Exiting..." << endl;
//            exit(1);
//        }
//        return false;
//    }
//
//
//
//    /*
//        parseGffLine: converts a lineVector into either Bed or BedCOV (templated, hence in header to avoid linker issues.)
//    */
//    template <typename T>
//    inline bool parseGffLine (T &bed, const vector<string> &lineVector, int lineNum, unsigned int numFields) {
//        if (numFields == this->bedType) {
//            if (this->bedType >= 8 && _isGff) {
//                bed.chrom = lineVector[0];
//                if (isInteger(lineVector[3]))
//                    bed.start  = atoi(lineVector[3].c_str());
//                if (isInteger(lineVector[4]))
//                    bed.end  = atoi(lineVector[4].c_str());
//                bed.name   = lineVector[2];
//                bed.score  = lineVector[5];
//                bed.strand = lineVector[6].c_str();
//                bed.otherFields.push_back(lineVector[1]);  // add GFF "source". unused in Bed
//                bed.otherFields.push_back(lineVector[7]);  // add GFF "fname". unused in Bed
//                // handle the optional 9th field.
//                if (this->bedType == 9)
//                    bed.otherFields.push_back(lineVector[8]);  // add GFF "group". unused in Bed
//                bed.start--;
//            }
//            else {
//                cerr << "Error: unexpected number of fields at line: " << lineNum <<
//                        ".  Verify that your files are TAB-delimited and that your GFF file has 8 or 9 fields.  Exiting..." << endl;
//                exit(1);
//            }
//            if (bed.start > bed.end) {
//                cerr << "Error: malformed GFF entry at line " << lineNum << ". Start was greater than end. Exiting." << endl;
//                exit(1);
//            }
//            if ( (bed.start < 0) || (bed.end < 0) ) {
//                cerr << "Error: malformed GFF entry at line " << lineNum << ". Coordinate detected that is < 1. Exiting." << endl;
//                exit(1);
//            }
//            return true;
//        }
//        else if (numFields == 1) {
//            cerr << "Only one GFF field detected: " << lineNum << ".  Verify that your files are TAB-delimited.  Exiting..." << endl;
//            exit(1);
//        }
//        else if ((numFields != this->bedType) && (numFields != 0)) {
//            cerr << "Differing number of GFF fields encountered at line: " << lineNum << ".  Exiting..." << endl;
//            exit(1);
//        }
//        else if ((numFields < 8) && (numFields != 0)) {
//            cerr << "TAB delimited GFF file with 8 or 9 fields is required at line: "<< lineNum << ".  Exiting..." << endl;
//            exit(1);
//        }
//        return false;
//    }
//
//
//public:
//
//    /*
//        reportBedTab
//
//        Writes the _original_ Bed entry with a TAB
//        at the end of the line.
//        Works for Bed3 - Bed6.
//    */
//    template <typename T>
//    inline void reportBedTab(const T &bed) {
//        
//        // if it is azeroLength feature, we need to
//        // correct the start and end coords to what they were
//        // in the original file
//        CHRPOS start = bed.start;
//        CHRPOS end   = bed.end;
//        if (bed.zeroLength) {
//            if (_isGff == false)
//                start++;
//            end--;
//        }
//        
//        // Bed
//        if (_isGff == false && _isVcf == false) {
//            if (this->bedType == 3) {
//                printf ("%s\t%d\t%d\t", bed.chrom.c_str(), start, end);
//            }
//            else if (this->bedType == 4) {
//                printf ("%s\t%d\t%d\t%s\t", bed.chrom.c_str(), start, end, bed.name.c_str());
//            }
//            else if (this->bedType == 5) {
//                printf ("%s\t%d\t%d\t%s\t%s\t", bed.chrom.c_str(), start, end, bed.name.c_str(),
//                                                bed.score.c_str());
//            }
//            else if (this->bedType == 6) {
//                printf ("%s\t%d\t%d\t%s\t%s\t%s\t", bed.chrom.c_str(), start, end, bed.name.c_str(),
//                                                    bed.score.c_str(), bed.strand.c_str());
//            }
//            else if (this->bedType > 6) {
//                printf ("%s\t%d\t%d\t%s\t%s\t%s\t", bed.chrom.c_str(), start, end, bed.name.c_str(),
//                                                    bed.score.c_str(), bed.strand.c_str());
//
//                vector<string>::const_iterator othIt = bed.otherFields.begin();
//                vector<string>::const_iterator othEnd = bed.otherFields.end();
//                for ( ; othIt != othEnd; ++othIt) {
//                    printf("%s\t", othIt->c_str());
//                }
//            }
//        }
//        // VCF
//        else if (_isGff == false && _isVcf == true) {
//            printf ("%s\t%d\t", bed.chrom.c_str(), start+1);
//
//            vector<string>::const_iterator othIt = bed.otherFields.begin();
//            vector<string>::const_iterator othEnd = bed.otherFields.end();
//            for ( ; othIt != othEnd; ++othIt) {
//                printf("%s\t", othIt->c_str());
//            }
//        }
//        // GFF
//        else if (_isGff == true) {
//            // "GFF-8"
//            if (this->bedType == 8) {
//                printf ("%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t", bed.chrom.c_str(), bed.otherFields[0].c_str(),
//                                                                 bed.name.c_str(), start+1, end,
//                                                                 bed.score.c_str(), bed.strand.c_str(),
//                                                                 bed.otherFields[1].c_str());
//            }
//            // "GFF-9"
//            else if (this->bedType == 9) {
//                printf ("%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t", bed.chrom.c_str(), bed.otherFields[0].c_str(),
//                                                                 bed.name.c_str(), start+1, end,
//                                                                 bed.score.c_str(), bed.strand.c_str(),
//                                                                 bed.otherFields[1].c_str(), bed.otherFields[2].c_str());
//            }
//        }
//    }
//
//
//
//    /*
//        reportBedNewLine
//
//        Writes the _original_ Bed entry with a NEWLINE
//        at the end of the line.
//        Works for Bed3 - Bed6.
//    */
//    template <typename T>
//    inline void reportBedNewLine(const T &bed) {
//        
//        // if it is azeroLength feature, we need to
//        // correct the start and end coords to what they were
//        // in the original file
//        CHRPOS start = bed.start;
//        CHRPOS end   = bed.end;
//        if (bed.zeroLength) {
//            if (_isGff == false)
//                start++;
//            end--;
//        }
//        
//        //Bed
//        if (_isGff == false && _isVcf == false) {
//            if (this->bedType == 3) {
//                printf ("%s\t%d\t%d\n", bed.chrom.c_str(), start, end);
//            }
//            else if (this->bedType == 4) {
//                printf ("%s\t%d\t%d\t%s\n", bed.chrom.c_str(), start, end, bed.name.c_str());
//            }
//            else if (this->bedType == 5) {
//                printf ("%s\t%d\t%d\t%s\t%s\n", bed.chrom.c_str(), start, end, bed.name.c_str(),
//                                                bed.score.c_str());
//            }
//            else if (this->bedType == 6) {
//                printf ("%s\t%d\t%d\t%s\t%s\t%s\n", bed.chrom.c_str(), start, end, bed.name.c_str(),
//                                                    bed.score.c_str(), bed.strand.c_str());
//            }
//            else if (this->bedType > 6) {
//                printf ("%s\t%d\t%d\t%s\t%s\t%s", bed.chrom.c_str(), start, end, bed.name.c_str(),
//                                                    bed.score.c_str(), bed.strand.c_str());
//
//                vector<string>::const_iterator othIt = bed.otherFields.begin();
//                vector<string>::const_iterator othEnd = bed.otherFields.end();
//                for ( ; othIt != othEnd; ++othIt) {
//                    printf("\t%s", othIt->c_str());
//                }
//                printf("\n");
//            }
//        }
//        // VCF
//        else if (_isGff == false && _isVcf == true) {
//            printf ("%s\t%d\t", bed.chrom.c_str(), start+1);
//
//            vector<string>::const_iterator othIt = bed.otherFields.begin();
//            vector<string>::const_iterator othEnd = bed.otherFields.end();
//            for ( ; othIt != othEnd; ++othIt) {
//                printf("%s\t", othIt->c_str());
//            }
//            printf("\n");
//        }
//        // GFF
//        else if (_isGff == true) {
//            // "GFF-8"
//            if (this->bedType == 8) {
//                printf ("%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\n", bed.chrom.c_str(), bed.otherFields[0].c_str(),
//                                                                 bed.name.c_str(), start+1, end,
//                                                                 bed.score.c_str(), bed.strand.c_str(),
//                                                                 bed.otherFields[1].c_str());
//            }
//            // "GFF-9"
//            else if (this->bedType == 9) {
//                printf ("%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\n", bed.chrom.c_str(), bed.otherFields[0].c_str(),
//                                                                 bed.name.c_str(), start+1, end,
//                                                                 bed.score.c_str(), bed.strand.c_str(),
//                                                                 bed.otherFields[1].c_str(), bed.otherFields[2].c_str());
//            }
//        }
//    }
//
//
//
//    /*
//        reportBedRangeNewLine
//
//        Writes a custom start->end for a Bed entry
//        with a NEWLINE at the end of the line.
//
//        Works for Bed3 - Bed6.
//    */
//    template <typename T>
//    inline void reportBedRangeTab(const T &bed, CHRPOS start, CHRPOS end) {
//        
//        // if it is azeroLength feature, we need to
//        // correct the start and end coords to what they were
//        // in the original file
//        if (bed.zeroLength) {
//            start = bed.start + 1;
//            end   = bed.end - 1;
//        }
//        // Bed
//        if (_isGff == false && _isVcf == false) {
//            if (this->bedType == 3) {
//                printf ("%s\t%d\t%d\t", bed.chrom.c_str(), start, end);
//            }
//            else if (this->bedType == 4) {
//                printf ("%s\t%d\t%d\t%s\t", bed.chrom.c_str(), start, end, bed.name.c_str());
//            }
//            else if (this->bedType == 5) {
//                printf ("%s\t%d\t%d\t%s\t%s\t", bed.chrom.c_str(), start, end, bed.name.c_str(),
//                                                bed.score.c_str());
//            }
//            else if (this->bedType == 6) {
//                printf ("%s\t%d\t%d\t%s\t%s\t%s\t", bed.chrom.c_str(), start, end, bed.name.c_str(),
//                                                    bed.score.c_str(), bed.strand.c_str());
//            }
//            else if (this->bedType > 6) {
//                printf ("%s\t%d\t%d\t%s\t%s\t%s\t", bed.chrom.c_str(), start, end, bed.name.c_str(),
//                                                    bed.score.c_str(), bed.strand.c_str());
//
//                vector<string>::const_iterator othIt = bed.otherFields.begin();
//                vector<string>::const_iterator othEnd = bed.otherFields.end();
//                for ( ; othIt != othEnd; ++othIt) {
//                    printf("%s\t", othIt->c_str());
//                }
//            }
//        }
//        // VCF
//        else if (_isGff == false && _isVcf == true) {
//            printf ("%s\t%d\t", bed.chrom.c_str(), bed.start+1);
//
//            vector<string>::const_iterator othIt = bed.otherFields.begin();
//            vector<string>::const_iterator othEnd = bed.otherFields.end();
//            for ( ; othIt != othEnd; ++othIt) {
//                printf("%s\t", othIt->c_str());
//            }
//        }
//        // GFF
//        else if (_isGff == true) {
//            // "GFF-8"
//            if (this->bedType == 8) {
//                printf ("%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t", bed.chrom.c_str(), bed.otherFields[0].c_str(),
//                                                             bed.name.c_str(), start+1, end,
//                                                             bed.score.c_str(), bed.strand.c_str(),
//                                                             bed.otherFields[1].c_str());
//            }
//            // "GFF-9"
//            else if (this->bedType == 9) {
//                printf ("%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t", bed.chrom.c_str(), bed.otherFields[0].c_str(),
//                                                             bed.name.c_str(), start+1, end,
//                                                             bed.score.c_str(), bed.strand.c_str(),
//                                                             bed.otherFields[1].c_str(), bed.otherFields[2].c_str());
//            }
//        }
//    }
//
//
//
//    /*
//        reportBedRangeTab
//
//        Writes a custom start->end for a Bed entry
//        with a TAB at the end of the line.
//
//        Works for Bed3 - Bed6.
//    */
//    template <typename T>
//    inline void reportBedRangeNewLine(const T &bed, CHRPOS start, CHRPOS end) {
//        
//        // if it is azeroLength feature, we need to
//        // correct the start and end coords to what they were
//        // in the original file
//        if (bed.zeroLength) {
//            start = bed.start + 1;
//            end   = bed.end - 1;
//        }
//        // Bed
//        if (_isGff == false && _isVcf == false) {
//            if (this->bedType == 3) {
//                printf ("%s\t%d\t%d\n", bed.chrom.c_str(), start, end);
//            }
//            else if (this->bedType == 4) {
//                printf ("%s\t%d\t%d\t%s\n", bed.chrom.c_str(), start, end, bed.name.c_str());
//            }
//            else if (this->bedType == 5) {
//                printf ("%s\t%d\t%d\t%s\t%s\n", bed.chrom.c_str(), start, end, bed.name.c_str(),
//                                                bed.score.c_str());
//            }
//            else if (this->bedType == 6) {
//                printf ("%s\t%d\t%d\t%s\t%s\t%s\n", bed.chrom.c_str(), start, end, bed.name.c_str(),
//                                                    bed.score.c_str(), bed.strand.c_str());
//            }
//            else if (this->bedType > 6) {
//                printf ("%s\t%d\t%d\t%s\t%s\t%s", bed.chrom.c_str(), start, end, bed.name.c_str(),
//                                                    bed.score.c_str(), bed.strand.c_str());
//
//                vector<string>::const_iterator othIt = bed.otherFields.begin();
//                vector<string>::const_iterator othEnd = bed.otherFields.end();
//                for ( ; othIt != othEnd; ++othIt) {
//                    printf("\t%s", othIt->c_str());
//                }
//                printf("\n");
//            }
//        }
//        // VCF
//        else if (_isGff == false && _isVcf == true) {
//            printf ("%s\t%d\t", bed.chrom.c_str(), bed.start+1);
//
//            vector<string>::const_iterator othIt = bed.otherFields.begin();
//            vector<string>::const_iterator othEnd = bed.otherFields.end();
//            for ( ; othIt != othEnd; ++othIt) {
//                printf("%s\t", othIt->c_str());
//            }
//            printf("\n");
//        }
//        // GFF
//        else if (_isGff == true) {
//            // "GFF-9"
//            if (this->bedType == 8) {
//                printf ("%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\n", bed.chrom.c_str(), bed.otherFields[0].c_str(),
//                                                             bed.name.c_str(), start+1, end,
//                                                             bed.score.c_str(), bed.strand.c_str(),
//                                                             bed.otherFields[1].c_str());
//            }
//            // "GFF-8"
//            else if (this->bedType == 9) {
//                printf ("%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\n", bed.chrom.c_str(), bed.otherFields[0].c_str(),
//                                                             bed.name.c_str(), start+1, end,
//                                                             bed.score.c_str(), bed.strand.c_str(),
//                                                             bed.otherFields[1].c_str(), bed.otherFields[2].c_str());
//            }
//        }
//    }
//
//
//    /*
//        reportNullBedTab
//    */
//    void reportNullBedTab() {
//
//        if (_isGff == false && _isVcf == false) {
//            if (this->bedType == 3) {
//                printf (".\t-1\t-1\t");
//            }
//            else if (this->bedType == 4) {
//                printf (".\t-1\t-1\t.\t");
//            }
//            else if (this->bedType == 5) {
//                printf (".\t-1\t-1\t.\t-1\t");
//            }
//            else if (this->bedType == 6) {
//                printf (".\t-1\t-1\t.\t-1\t.\t");
//            }
//            else if (this->bedType > 6) {
//                printf (".\t-1\t-1\t.\t-1\t.\t");
//                for (unsigned int i = 6; i < this->bedType; ++i) {
//                    printf(".\t");
//                }
//            }
//        }
//        else if (_isGff == true && _isVcf == false) {
//            if (this->bedType == 8) {
//                printf (".\t.\t.\t-1\t-1\t-1\t.\t.\t");
//            }
//            else if (this->bedType == 9) {
//                printf (".\t.\t.\t-1\t-1\t-1\t.\t.\t.\t");
//            }
//        }
//    }
//
//
//    /*
//        reportNullBedTab
//    */
//    void reportNullBedNewLine() {
//
//        if (_isGff == false && _isVcf == false) {
//            if (this->bedType == 3) {
//                printf (".\t-1\t-1\n");
//            }
//            else if (this->bedType == 4) {
//                printf (".\t-1\t-1\t.\n");
//            }
//            else if (this->bedType == 5) {
//                printf (".\t-1\t-1\t.\t-1\n");
//            }
//            else if (this->bedType == 6) {
//                printf (".\t-1\t-1\t.\t-1\t.\n");
//            }
//            else if (this->bedType > 6) {
//                printf (".\t-1\t-1\t.\t-1\t.");
//                for (unsigned int i = 6; i < this->bedType; ++i) {
//                    printf("\t.");
//                }
//                printf("\n");
//            }
//        }
//        else if (_isGff == true && _isVcf == false) {
//            if (this->bedType == 8) {
//                printf (".\t.\t.\t-1\t-1\t-1\t.\t.\n");
//            }
//            else if (this->bedType == 9) {
//                printf (".\t.\t.\t-1\t-1\t-1\t.\t.\t.\n");
//            }
//        }
//    }
//
//
//};

#endif /* BEDUTILS_H */
