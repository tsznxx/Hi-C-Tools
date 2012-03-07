#include "bedUtils.h"

using namespace std;

//////////////////////////////////////
// Bed
//////////////////////////////////////

// constructors
// Bed6
Bed::Bed (const string & chrom, const CHRPOS & start, const CHRPOS & end, const string & name, const float &score, const string & strand, const string &desc):
        chrom (chrom),
        start (start),
        end (end),
        name (name),
        score (score),
        strand (strand),
		desc(desc)
{
}

Bed::Bed (const string & chrom, const CHRPOS & start, const CHRPOS & end, const string & name, const float &score, const string & strand, const vector < string > &elems, const string &desc):
        chrom (chrom),
        start (start),
        end (end),
        name (name),
        score (score),
        strand (strand),
		desc(desc)
{
	otherFields.assign (elems.begin (), elems.end ());
}


// Bed from Vector
Bed::Bed (const vector < string > &elems, const string &desc)
{
  int l;
  l = elems.size ();
  if (l < 3)
  {
      cerr << "Vector doesn't contain enough elements!" << endl;
      exit (1);
  }
  chrom = elems[0];
  start = atoi (elems[1].c_str ());
  end = atoi (elems[2].c_str ());
  name = l > 3 ? elems[3] : "";
  score = l > 4 ? atof (elems[4].c_str ()) : 0.0;
  strand = l > 5 ? elems[5] : ".";

  if (l > 6)
    otherFields.assign (elems.begin () + 6, elems.end ());
  this->desc=desc;
}

/*
// Set Bed6
void Bed::setBed (const string & chrom, const CHRPOS & start, const CHRPOS & end, const string & name, const float &score, const string & strand, const string &desc)
{

	this->chrom  = chrom;
	this->start  = start;
	this-> end   = end;
	this->name   = name;
	this->score  = score;
	this->strand = strand;
	this->desc   = desc;
}

// Set Bed more
void Bed::setBed (const string & chrom, const CHRPOS & start, const CHRPOS & end,
        const string & name, const float &score, const string & strand,
        const vector < string > &elems, const string &desc)
{
	this->chrom  = chrom;
	this->start  = start;
	this-> end   = end;
	this->name   = name;
	this->score  = score;
	this->strand = strand;
	this->desc   = desc;

	otherFields.assign (elems.begin (), elems.end ());
}

// Set Bed from vector
void Bed::setBed (const vector < string > &elems, const string &desc)
{
  int l;
  l = elems.size ();
  if (l < 3)
  {
      cerr << "Vector doesn't contain enough elements!" << endl;
      exit (1);
  }
  chrom = elems[0];
  start = atoi (elems[1].c_str ());
  end = atoi (elems[2].c_str ());
  name = l > 3 ? elems[3] : "";
  score = l > 4 ? atof (elems[4].c_str ()) : 0.0;
  strand = l > 5 ? elems[5] : ".";

  if (l > 6)
    otherFields.assign (elems.begin () + 6, elems.end ());
  this->desc=desc;
}
*/

//Desctructor
Bed::~Bed (void)
{
}


/***************************************
 BedUtils
***************************************/

void BedUtils::bowtieToBed(Bed &tbed, const LINE_ELEMS &elems)
{
	//Bed([x[2],x[3],x[3]+len(x[4]),x[0],1,x[1]])	
	tbed.chrom  = elems[2];
	tbed.start  = ToValue<CHRPOS>(elems[3]);
	tbed.end    = tbed.start+elems[4].size();
	tbed.name   = elems[0];
	tbed.strand = elems[1];
}

// Load Bed file into BedVector
void loadBedFileToVec( BedVec
// BedMap related
// Load Bed file into BedMap
void BedUtils::loadBedFileToMap   (BedMap &bedmap, const string &bedfile)
{
	LINE_ELEMS elems;
	ColumnReader bedreader(bedfile);

	bedreader.Open();
	while(bedreader.getNext(elems)!=LINE_INVALID)
	{
		Bed tbed(elems);
		bedmap[tbed.chrom][tbed.getBIN()].push_back(tbed);
	}
}

// Load BedVector into BedMap
void BedUtils::loadBedVectorToMap   (BedMap &bedmap, const BedVector & bedvector)
{
	for(BedVector::const_iterator it=bedvector.begin(); it!=bedvector.end();it++)
		bedmap[it->chrom][it->getBIN()].push_back(*it);
}

// Intersect Bed
void intersectBed (const Bed &tbed, BedMap &bedmap,BedVector &Hits, const bool &forcestrand)
{
	CHRPOS overlap,maxoverlap=0;
	BIN startBin,endBin,offset;
	BedVector::const_iterator bedItr;
	BedVector::const_iterator bedEnd;

	startBin = tbed.start   >> _binFirstShift;
	endBin   = (tbed.end-1) >> _binFirstShift;

	for (BINLEVEL i = 0; i < _binLevels; ++i) // while(startBin!=endBin)
	{
		offset = _binOffsetsExtended[i];
		for (BIN j = (startBin+offset); j <= (endBin+offset); ++j)
		{
			bedItr = bedmap[tbed.chrom][j].begin();
			bedEnd = bedmap[tbed.chrom][j].end();
			for(;bedItr != bedEnd; bedItr++)
			{
				overlap = tbed.overlapLength(*bedItr, forcestrand);
				if(overlap)
				{
					if(overlap>maxoverlap) // Make sure the best hit is the first element in Hits.
					{
						maxoverlap=overlap;
						Hits.insert(Hits.begin(),*bedItr);
					}
					else
						Hits.push_back(*bedItr);
				}
			}
		}
		startBin >>= _binNextShift;
		endBin >>= _binNextShift;
	}
}
















