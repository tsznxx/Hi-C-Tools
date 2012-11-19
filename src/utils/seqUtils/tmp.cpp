#include <iostream>
#include <list>
#include <algorithm>
#include <cctype>

#include "seqUtils.h"

using namespace std;


#ifdef _MSC_VER
	#define ftello _ftelli64
	#define fseeko _fseeki64
#endif


/*###Class Fasta###
*************************************/

Fasta::Fasta(const string &tid, const string &tseq):id(tid),seq(tseq){
}

// Format output (lineSize=0: print the whole sequence in one line
void Fasta::print(CHRPOS lineSize){
	int i,curlen;

	curlen=seq.size()-lineSize;
	cout << ">" << id << endl;
	for(i=0; i< curlen;i+=lineSize)
	{
		cout << seq.substr(i,lineSize) << endl;
	}
	cout << seq.substr(i) << endl;
}

/*###Class Fastq###
*************************************/

// Constructor
Fastq::Fastq(const string &tid, const string &tseq, const string &tphred):Fasta(tid,tseq){
	phred=tphred;
}

/*###Class Genome###
*********************************/

struct Genome::FastaPrivate {
  
    struct GenomeIndexData {
        string  Name;
        int32_t Length;
        int64_t Offset;
        int32_t LineLength;
        int32_t ByteLength; // LineLength + newline character(s) - varies on OS where file was generated
    };
  
    // data members
    FILE* Stream;
    bool IsOpen;
    
    FILE* IndexStream;
    bool HasIndex;
    bool IsIndexOpen;
  
    vector<GenomeIndexData> Index;
    
    // ctor
    FastaPrivate(void);
    ~FastaPrivate(void);
    
    // 'public' API methods
    bool Close(void);
    bool CreateIndex(const string& indexFilename);
    bool GetBase(const int& refId, const int& position, char& base);
    bool GetSequence(const int& refId, const int& start, const int& stop, string& sequence);
    bool Open(const string& filename, const string& indexFilename);
    
    // internal methods
    private:
        void Chomp(char* sequence);
        bool GetNameFromHeader(const string& header, string& name);
        bool GetNextHeader(string& header);
        bool GetNextSequence(string& sequence);
        bool LoadIndexData(void);
        bool Rewind(void);
        bool WriteIndexData(void);
};

Genome::FastaPrivate::FastaPrivate(void) : IsOpen(false), HasIndex(false), IsIndexOpen(false){
}

Genome::FastaPrivate::~FastaPrivate(void) {
    Close();
}

// remove any trailing newlines
void Genome::FastaPrivate::Chomp(char* sequence) {
    static const int CHAR_LF = 10;
    static const int CHAR_CR = 13;
  
    size_t seqLength = strlen(sequence);
    if ( seqLength == 0 ) return;
    --seqLength; // ignore null terminator
  
    while ( sequence[seqLength] == CHAR_LF || 
            sequence[seqLength] == CHAR_CR 
          ) 
    {
        sequence[seqLength] = 0;
        --seqLength;
        if (seqLength < 0) 
            break;
    }
}

bool Genome::FastaPrivate::Close(void) { 
    // close fasta file
    if ( IsOpen ) {
        fclose(Stream);
        IsOpen = false;
    }

    // close index file
    if ( HasIndex && IsIndexOpen ) {
        fclose(IndexStream);
        HasIndex = false;
        IsIndexOpen = false;
    }
  
    // return success
    return true;
}

bool Genome::FastaPrivate::CreateIndex(const string& indexFilename) {  
    // check that file is open
    if ( !IsOpen ) {
        cerr << "FASTA error : cannot create index, FASTA file not open" << endl;
        return false;
    }
  
    // rewind FASTA file
    if ( !Rewind() ) {
        cerr << "FASTA error : could not rewind FASTA file" << endl;
        return false;
    }
    
    // clear out prior index data
    Index.clear();
    
    // -------------------------------------------
    // calculate lineLength & byteLength
    
    int lineLength = 0;
    int byteLength = 0;
    
    // skip over header
    char buffer[1024];
    if ( fgets(buffer, 1024, Stream) == 0 ) {
        cerr << "FASTA error : could not read from file" << endl;
        return false;
    }
    if ( feof(Stream) ) return false;
    if ( buffer[0] != '>' ) { 
        cerr << "FASTA error : expected header ('>'), instead : " << buffer[0] << endl;
        return false;
    }
  
    // read in first line of sequence  
    char c = fgetc(Stream);
    while ( (c >= 0) && (c != '\n') ) {
        ++byteLength;
        if (isgraph(c)) ++lineLength;
        c = fgetc(Stream);
    }
    ++byteLength; // store newline
    
    // rewind FASTA file
    if ( !Rewind() ) {
        cerr << "FASTA error : could not rewind FASTA file" << endl;
        return false;
    }
    
    // iterate through fasta entries
    int currentId   = 0;
    string header   = "";
    string sequence = "";
    while ( GetNextHeader(header) ) {
        
        // ---------------------------
        // build index entry data
        GenomeIndexData data;
        
        // store file offset of beginning of DNA sequence (after header)
        data.Offset = ftello(Stream);
        
        // parse header, store sequence name in data.Name
        if ( !GetNameFromHeader(header, data.Name) ) {
            cerr << "FASTA error : could not parse read name from FASTA header" << endl;
            return false;
        }
        
        // retrieve FASTA sequence
        if ( !GetNextSequence(sequence) ) {
            cerr << "FASTA error : could not read in next sequence from FASTA file" << endl;
            return false;
        }
        
        // store sequence length & line/byte lengths
        data.Length = sequence.length();
        data.LineLength = lineLength;
        data.ByteLength = byteLength;
        
        // store index entry
        Index.push_back(data);
        
        // update ref Id
        ++currentId;
    }
    
    // open index file
    if ( !indexFilename.empty() ) {
        IndexStream = fopen(indexFilename.c_str(), "wb");
        if ( !IndexStream ) {
            cerr << "FASTA error : Could not open " << indexFilename << " for writing." << endl;
            return false;
        }
        IsIndexOpen = true;
    }
    
    // write index data
    if ( !WriteIndexData() ) return false;
    HasIndex = true;
    
    // close index file
    fclose(IndexStream);
    IsIndexOpen = false;
    
    // return succes status
    return true;
}

bool Genome::FastaPrivate::GetBase(const int& refId, const int& position, char& base) {  
    // make sure FASTA file is open
    if ( !IsOpen ) {
        cerr << "FASTA error : file not open for reading" << endl;
        return false;
    }
  
    // use index if available
    if ( HasIndex && !Index.empty() ) {
        
        // validate reference id 
        if ( (refId < 0) || (refId >= (int)Index.size()) ) {
            cerr << "FASTA error: invalid refId specified: " << refId << endl;
            return false;
        }
        
        // retrieve reference index data
        const GenomeIndexData& referenceData = Index.at(refId);
        
        // validate position 
        if ( (position < 0) || (position > referenceData.Length) ) {
            cerr << "FASTA error: invalid position specified: " << position << endl;
            return false;
        }
        
        // seek to beginning of sequence data
        if ( fseeko(Stream, referenceData.Offset, SEEK_SET) != 0 ) {
            cerr << "FASTA error : could not sek in file" << endl;
            return false;
        }
      
        // retrieve sequence
        string sequence = "";
        if ( !GetNextSequence(sequence) ) {
            cerr << "FASTA error : could not retrieve base from FASTA file" << endl;
            return false;
        }
        
        // set base & return success
        base = sequence.at(position);
        return true;
    }
    
    // else plow through sequentially
    else {
      
        // rewind FASTA file
        if ( !Rewind() ) {
            cerr << "FASTA error : could not rewind FASTA file" << endl;
            return false;
        }
        
        // iterate through fasta entries
        int currentId = 0;
        string header = "";
        string sequence = "";
        
        // get first entry
        GetNextHeader(header);
        GetNextSequence(sequence);
        
        while ( currentId != refId ) {
            GetNextHeader(header);
            GetNextSequence(sequence);
            ++currentId;
        }
        
        // get desired base from sequence 
        // TODO: error reporting on invalid position
        if ( currentId == refId && (sequence.length() >= (size_t)position) ) {          
            base = sequence.at(position);
            return true;
        }
      
        // could not get sequence
        return false;
    }
 
    // return success
    return true;
}

bool Genome::FastaPrivate::GetNameFromHeader(const string& header, string& name) {
    // get rid of the leading greater than sign
    string s = header.substr(1);

    // extract the first non-whitespace segment
    char* pName = (char*)s.data();
    unsigned int nameLen = (unsigned int)s.size();

    unsigned int start = 0;
    while ( (pName[start] == 32) || (pName[start] == 9) || (pName[start] == 10) || (pName[start] == 13) ) {
        start++;
        if ( start == nameLen ) 
            break;
    }

    unsigned int stop  = start;
    if ( stop < nameLen ) {
        while( (pName[stop] != 32) && (pName[stop] != 9) && (pName[stop] != 10) && (pName[stop] != 13) ) {
            stop++;
            if ( stop == nameLen ) 
                break;
        }
    }

    if ( start == stop ) {
        cerr << "FASTA error : could not parse read name from FASTA header" << endl;
        return false;
    }

    name = s.substr(start, stop - start).c_str();
    return true;
}

bool Genome::FastaPrivate::GetNextHeader(string& header) {
  
    // validate input stream
    if ( !IsOpen || feof(Stream) ) 
        return false;
    
    // read in header line
    char buffer[1024];
    if ( fgets(buffer, 1024, Stream) == 0 ) {
        cerr << "FASTA error : could not read from file" << endl;
        return false;
    }
    
    // make sure it's a FASTA header
    if ( buffer[0] != '>' ) { 
        cerr << "FASTA error : expected header ('>'), instead : " << buffer[0] << endl;
        return false;
    }
  
    // import buffer contents to header string
    stringstream headerBuffer("");
    headerBuffer << buffer;
    header = headerBuffer.str();
  
    // return success
    return true;
}

bool Genome::FastaPrivate::GetNextSequence(string& sequence) {
  
    // validate input stream
    if ( !IsOpen || feof(Stream) ) 
        return false;
    
    // read in sequence  
    char buffer[1024];
    ostringstream seqBuffer("");
    while(true) {
        
        char ch = fgetc(Stream);
        ungetc(ch, Stream);
        if( (ch == '>') || feof(Stream) ) 
              break;       
        
        if ( fgets(buffer, 1024, Stream) == 0 ) {
            cerr << "FASTA error : could not read from file" << endl;
            return false;
        }
        
        Chomp(buffer);
        seqBuffer << buffer;
    }
    
    // import buffer contents to sequence string
    sequence = seqBuffer.str();
  
    // return success
    return true;
}

bool Genome::FastaPrivate::GetSequence(const int& refId, const int& start, const int& stop, string& sequence) {
 
    // make sure FASTA file is open
    if ( !IsOpen ) {
        cerr << "FASTA error : file not open for reading" << endl;
        return false;
    }
  
    // use index if available
    if ( HasIndex && !Index.empty() ) {
      
        // validate reference id 
        if ( (refId < 0) || (refId >= (int)Index.size()) ) {
            cerr << "FASTA error: invalid refId specified: " << refId << endl;
            return false;
        }
        
        // retrieve reference index data
        const GenomeIndexData& referenceData = Index.at(refId);
        
        // validate stop position 
        if ( (start < 0) || (start > stop) || (stop > referenceData.Length) ) {
            cerr << "FASTA error: invalid start/stop positions specified: " << start << ", " << stop << endl;
            return false;
        }
        
        // seek to beginning of sequence data
        if ( fseeko(Stream, referenceData.Offset, SEEK_SET) != 0 ) {
            cerr << "FASTA error : could not sek in file" << endl;
            return false;
        }
      
        // retrieve full sequence
        string fullSequence = "";
        if ( !GetNextSequence(fullSequence) ) {
            cerr << "FASTA error : could not retrieve sequence from FASTA file" << endl;
            return false;
        }
        
        // set sub-sequence & return success
        const int seqLength = (stop - start) + 1;
        sequence = fullSequence.substr(start, seqLength);
        return true;
    }
    
    // else plow through sequentially
    else {
      
        // rewind FASTA file
        if ( !Rewind() ) {
            cerr << "FASTA error : could not rewind FASTA file" << endl;
            return false;
        }
     
        // iterate through fasta entries
        int currentId = 0;
        string header = "";
        string fullSequence = "";
        
        // get first entry
        GetNextHeader(header);
        GetNextSequence(fullSequence);
        
        while ( currentId != refId ) {
            GetNextHeader(header);
            GetNextSequence(fullSequence);
            ++currentId;
        }
        
        // get desired substring from sequence
        // TODO: error reporting on invalid start/stop positions
        if ( currentId == refId && (fullSequence.length() >= (size_t)stop) ) {          
            const int seqLength = (stop - start) + 1;
            sequence = fullSequence.substr(start, seqLength);
            return true;
        }
      
        // could not get sequence
        return false;
    }
  
    // return success
    return true;
}

bool Genome::FastaPrivate::LoadIndexData(void) {
  
    // skip if no index file available
    if ( !IsIndexOpen ) return false; 
  
    // clear any prior index data
    Index.clear();
  
    char buffer[1024];
    stringstream indexBuffer;
    while ( true ) {
      
        char c = fgetc(IndexStream);
        if ( (c == '\n') || feof(IndexStream) ) break;
        ungetc(c, IndexStream);
      
        // clear index buffer
        indexBuffer.str("");
        
        // read line from index file
        if ( fgets(buffer, 1024, IndexStream) == 0 ) {
            cerr << "FASTA LoadIndexData() error : could not read from index file" << endl;
            HasIndex = false;
            return false;
        }
      
        // store line in indexBuffer
        indexBuffer << buffer;
        
        // retrieve fasta index data from line
        GenomeIndexData data;
        indexBuffer >> data.Name;
        indexBuffer >> data.Length;
        indexBuffer >> data.Offset;
        indexBuffer >> data.LineLength;
        indexBuffer >> data.ByteLength;
        
        // store index entry
        Index.push_back(data);
    }
    
    return true;
}

bool Genome::FastaPrivate::Open(const string& filename, const string& indexFilename) {
 
    bool success = true;
  
    // open FASTA filename
    Stream = fopen(filename.c_str(), "rb");
    if ( !Stream ) {
        cerr << "FASTA error: Could not open " << filename << " for reading" << endl;
        return false;
    }
    IsOpen = true;
    success &= IsOpen;
    
    // open index file if it exists
    if ( !indexFilename.empty() ) {
        IndexStream = fopen(indexFilename.c_str(), "rb");
        if ( !IndexStream ) {
            cerr << "FASTA error : Could not open " << indexFilename << " for reading." << endl;
            return false;
        }
        IsIndexOpen = true;
        success &= IsIndexOpen;
        
        // attempt to load index data
        HasIndex = LoadIndexData();
        success &= HasIndex;
    }
    
    // return success status
    return success;
}

bool Genome::FastaPrivate::Rewind(void) {
    if ( !IsOpen ) return false;
    return ( fseeko(Stream, 0, SEEK_SET) == 0 );
}

bool Genome::FastaPrivate::WriteIndexData(void) {
 
    // skip if no index file available
    if ( !IsIndexOpen ) return false; 
  
    // iterate over index entries
    bool success = true;
    stringstream indexBuffer;
    vector<GenomeIndexData>::const_iterator indexIter = Index.begin();
    vector<GenomeIndexData>::const_iterator indexEnd  = Index.end();
    for ( ; indexIter != indexEnd; ++indexIter ) {
      
        // clear stream
        indexBuffer.str("");
      
        // write data to stream
        const GenomeIndexData& data = (*indexIter);
        indexBuffer << data.Name << "\t"
                    << data.Length << "\t"
                    << data.Offset << "\t"
                    << data.LineLength << "\t"
                    << data.ByteLength << endl;
                    
        // write stream to file
        success &= ( fputs(indexBuffer.str().c_str(), IndexStream) >= 0 );
    }
  
    // return success status
    return success;
}

// --------------------------------
// Genome implementation

Genome::Genome(void) {
    d = new FastaPrivate;
}

Genome::~Genome(void) {
    delete d;
    d = 0;
}

bool Genome::Close(void) { 
    return d->Close();
}

bool Genome::CreateIndex(const string& indexFilename) {
    return d->CreateIndex(indexFilename);
}

bool Genome::GetBase(const int& refId, const int& position, char& base) {
    return d->GetBase(refId, position, base);
}

bool Genome::GetSequence(const int& refId, const int& start, const int& stop, string& sequence) {
    return d->GetSequence(refId, start, stop, sequence);
}

bool Genome::Open(const string& filename, const string& indexFilename) {
    return d->Open(filename, indexFilename);
}

/*###Class SeqReader###
************************************/

// Constructor
SeqReader::SeqReader(const string &fname):Reader(fname){
}

// Destructor
SeqReader::~SeqReader(void){
	close();
}

// Get next fasta record
bool SeqReader::getNext(Fasta &fa){
	// Initiation
	string curline,tseq;
	static string tid;
	list<string> seqlst;
	bool flag=true;
	
	// Get next Fasta
	while (flag)
	{
		if(_inStream->good())
			getline(*_inStream,curline);
		else // The last run
		{
			curline=">";
			flag=false;
		}

		if (*(curline.begin())=='>')
		{
			if (tid!="")
			{
				fa.id=tid;
				tseq.reserve(seqlst.size()* seqlst.begin()->size());
				for(list<string>::iterator it=seqlst.begin(); it!=seqlst.end();it++)
					tseq.append(*it);
				fa.seq=tseq;
				tid=curline.substr(1);
				return true;
			}
			else // The first run
				tid=curline.substr(1);
		}
		else
		{
			transform(curline.begin(), curline.end(), curline.begin(), (int (*)(int))toupper);
			seqlst.push_back(curline);
		}
	}
	if (tid!="" && tseq!="") // The last run
	{
		fa.id=tid;
		tseq.reserve(seqlst.size()* seqlst.begin()->size());
		for(list<string>::iterator it=seqlst.begin(); it!=seqlst.end();it++)
			tseq.append(*it);
		fa.seq=tseq;
		return true;
	}
	return false;
}

// Get next fastq record
bool SeqReader::getNext(Fastq &fq){
	// Initiations
	string curline;

	// Read new record
	if (_inStream->good())
	{
		getline(*_inStream,curline);
		if (curline.size())
		{
			fq.id=curline.substr(1);
		}
		else
			return false;
		getline(*_inStream,fq.seq);
		getline(*_inStream,curline);
		getline(*_inStream,fq.phred);
		return true;
	}
	else
		return false;
}

/*Namespace SeqUtils
******************************************/

// DNA/RNA to uppercase
string SeqUtils::toUpper(const string &seq){
	string upper;
	transform(seq.begin(),seq.end(),upper.begin(),(int (*)(int)) toupper );
	return upper;
}

// DNA/RNA to lowercase
string SeqUtils::toLower(const string &seq){
	string lower;
	transform(seq.begin(),seq.end(),lower.begin(),(int (*)(int)) tolower );
	return lower;
}

// DNA to RNA
string SeqUtils::toRNA  (const string &seq){
	string rna=seq;
	size_t found;

	found=rna.find_first_of("tT");
	while(found!=string::npos)
	{
		rna[found++]++; // 't' to 'u' or 'T' to 'U'
		found=rna.find_first_of("tT",found);
	}
	return rna;
}

// RNA to DNA
string SeqUtils::toDNA  (const string &seq){
	string dna=seq;
	size_t found;

	found=dna.find_first_of("uU");
	while(found!=string::npos)
	{
		dna[found++]--; // 'u' to 't' or 'U' to 'T'
		found=dna.find_first_of("uU",found);
	}
	return dna;
}

// reverse complermentary table
char transtable(char i){
	static string trans="ACGNTacgntTGCNAtgcna";
	static size_t found;

	found=trans.find_first_of(i);
	if(found!=string::npos)
		return trans[found+10];
	else
		return 'N';
}

// reverse complementary sequence
string SeqUtils::rc (const string &seq){
	string rcseq=seq;
	reverse(rcseq.begin(),rcseq.end());
	transform(rcseq.begin(),rcseq.end(),rcseq.begin(),transtable);
	return rcseq;
}

// GC content
float SeqUtils::GC (const string &seq){
	size_t gc=0;
	size_t found;
	found=seq.find_first_of("CGcg");
	while(found!=string::npos)
	{
		gc++;
		found=seq.find_first_of("CGcg",found+1);
	}
	return 100.0*gc/seq.size();
}

// TM value
float SeqUtils::TM(const string &seq){
	// Add code here.
	return 60;
}

// Molecular weight
float SeqUtils::MW(const string &seq){
	// Add code here.
	return 3000;
}

