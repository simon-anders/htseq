/* File: DynamicStepVector.i */
%module DynamicStepVector
%include "std_string.i"
%include "std_vector.i"
%include "cpointer.i"
%{
#define SWIG_FILE_WITH_INIT
#include "bamtools/BamReader.h"
#include "DynamicStepVector.hpp"
//HTSeq::DSV< long int, int >;
//HTSeq::DSV< long int, bool >;
//HTSeq::DSV< long int, double >;
%}

namespace std {
   %template(strvector) vector<string>;
   %template(intvector) vector<int>;
   %template(dblvector) vector<double>;
};


template< typename TValue >
std::string display_vector( std::vector< TValue > const & );
   
template< typename TKey, typename TValue >
class DSV{
public:
    typedef std::map< TKey, Value< TValue >* > Map;
    typedef SingleValue< TValue > TSV;
    typedef MultipleValue< TValue > TMV;

    DSV(void);
    ~DSV(void);

    DSV( size_t t );

    DSV( DSV< TKey, TValue > const & other );
    
    void add( TKey const & from, TKey const & to, TValue const & offset );

    void clear();

    void invert( TKey const & from, TKey const & to, TValue const & offset );
    
    TValue get( TKey const & key );
    
    std::vector< TValue > get( TKey const & from, TKey const & to );

    template< typename TFkt >
    void apply( TKey const & from, TKey const & to, TFkt & fkt );
    
    void set( TKey const & key, TValue const & val );
    
    void set( TKey const & from, TKey const & to, TValue const & val );
    
    void refurbish( TKey const & key );
    
    typename Map::iterator get_iter( TKey const & key );
            
    std::string info();
    
    std::vector< std::string > rinfo();
    
    size_t size();
    
    Map const & get_steps();
    
private:
    size_t threshold;
    std::map< TKey, Value< TValue >* > steps;
};

%template(intDSV) DSV< long int, int >;
//%template(boolDSV) DSV< long int, bool >;
%template(floatDSV) DSV< long int, double >;

namespace BamTools {

class BamReader {

    // constructor / destructor
    public:
        BamReader(void);
        ~BamReader(void);

    // public interface
    public:

        // ----------------------
        // BAM file operations
        // ----------------------

        // close BAM file
        void Close(void);
        // performs random-access jump to reference, position
        bool Jump(int refID, int position = 0);
        // opens BAM file (and optional BAM index file, if provided)
        void Open(const std::string& filename, const std::string& indexFilename = "");
        // returns file pointer to beginning of alignments
        bool Rewind(void);
        // sets a region of interest (with left & right bound reference/position)
        // attempts a Jump() to left bound as well
        // returns success/failure of Jump()
        bool SetRegion(const BamTools::BamRegion& region);
        bool SetRegion(const int& leftRefID, const int& leftBound, const int& rightRefID, const int& rightBound);

        // ----------------------
        // access alignment data
        // ----------------------

        // retrieves next available alignment (returns success/fail)
        bool GetNextAlignment(BamTools::BamAlignment& bAlignment);
        
        // retrieves next available alignment core data (returns success/fail)
        // ** DOES NOT parse any character data (bases, qualities, tag data)
        //    these can be accessed, if necessary, from the supportData 
        // useful for operations requiring ONLY positional or other alignment-related information
        bool GetNextAlignmentCore(BamTools::BamAlignment& bAlignment);

        // ----------------------
        // access auxiliary data
        // ----------------------

        // returns SAM header text
        const std::string GetHeaderText(void) const;
        // returns number of reference sequences
        int GetReferenceCount(void) const;
        // returns vector of reference objects
        const BamTools::RefVector GetReferenceData(void) const;
        // returns reference id (used for BamReader::Jump()) for the given reference name
        int GetReferenceID(const std::string& refName) const;
        // returns the name of the file associated with this BamReader
        const std::string GetFilename(void) const;

        // ----------------------
        // BAM index operations
        // ----------------------

        // creates index for BAM file, saves to file (default = bamFilename + ".bai")
        bool CreateIndex(void);

    // private implementation
    private:
        struct BamReaderPrivate;
        BamReaderPrivate* d;
};

struct BamAlignment {

    // constructors & destructor
    public:
        BamAlignment(void);
        BamAlignment(const BamAlignment& other);
        ~BamAlignment(void);

    // Queries against alignment flags
    public:        
        bool IsDuplicate(void) const;		// Returns true if this read is a PCR duplicate       
        bool IsFailedQC(void) const;		// Returns true if this read failed quality control      
        bool IsFirstMate(void) const;		// Returns true if alignment is first mate on read        
        bool IsMapped(void) const;		// Returns true if alignment is mapped        
        bool IsMateMapped(void) const;		// Returns true if alignment's mate is mapped        
        bool IsMateReverseStrand(void) const;	// Returns true if alignment's mate mapped to reverse strand        
        bool IsPaired(void) const;		// Returns true if alignment part of paired-end read        
        bool IsPrimaryAlignment(void) const;	// Returns true if reported position is primary alignment       
        bool IsProperPair(void) const;		// Returns true if alignment is part of read that satisfied paired-end resolution     
        bool IsReverseStrand(void) const;	// Returns true if alignment mapped to reverse strand
        bool IsSecondMate(void) const;		// Returns true if alignment is second mate on read

    // Manipulate alignment flags
    public:        
        void SetIsDuplicate(bool ok);		// Sets "PCR duplicate" flag        
        void SetIsFailedQC(bool ok);		// Sets "failed quality control" flag        
        void SetIsFirstMate(bool ok);		// Sets "alignment is first mate" flag        
        void SetIsMateUnmapped(bool ok);	// Sets "alignment's mate is mapped" flag        
        void SetIsMateReverseStrand(bool ok);	// Sets "alignment's mate mapped to reverse strand" flag        
        void SetIsPaired(bool ok);		// Sets "alignment part of paired-end read" flag        
	void SetIsProperPair(bool ok);		// Sets "alignment is part of read that satisfied paired-end resolution" flag        
        void SetIsReverseStrand(bool ok);	// Sets "alignment mapped to reverse strand" flag        
        void SetIsSecondaryAlignment(bool ok);	// Sets "position is primary alignment" flag        
        void SetIsSecondMate(bool ok);		// Sets "alignment is second mate on read" flag        
        void SetIsUnmapped(bool ok);		// Sets "alignment is mapped" flag

    // Tag data access methods
    public:
        bool GetEditDistance(uint8_t& editDistance) const;	// get "NM" tag data - contributed by Aaron Quinlan
        bool GetReadGroup(std::string& readGroup) const;	// get "RG" tag data
        
        bool GetTag(const std::string& tag, std::string& destination);
        template<typename T> bool GetTag(const std::string& tag, T& destination);

    // Additional data access methods
    public:
	int GetEndPosition(bool usePadded = false) const;	// calculates alignment end position, based on starting position and CIGAR operations

    // 'internal' utility methods 
    private:
        static void SkipToNextTag(const char storageType, char* &pTagData, unsigned int& numBytesParsed);

    // Data members
    public:
        std::string     Name;              // Read name
        int             Length;            // Query length
        std::string     QueryBases;        // 'Original' sequence (as reported from sequencing machine)
        std::string     AlignedBases;      // 'Aligned' sequence (includes any indels, padding, clipping)
        std::string     Qualities;         // FASTQ qualities (ASCII characters, not numeric values)
        std::string     TagData;           // Tag data (accessor methods will pull the requested information out)
        int             RefID;             // ID number for reference sequence
        int             Position;          // Position (0-based) where alignment starts
        unsigned short  Bin;               // Bin in BAM file where this alignment resides
        unsigned short  MapQuality;        // Mapping quality score
        unsigned int    AlignmentFlag;     // Alignment bit-flag - see Is<something>() methods to query this value, SetIs<something>() methods to manipulate 
        std::vector<BamTools::CigarOp> CigarData; // CIGAR operations for this alignment
        int             MateRefID;         // ID number for reference sequence where alignment's mate was aligned
        int             MatePosition;      // Position (0-based) where alignment's mate starts
        int             InsertSize;        // Mate-pair insert size
        
        struct BamAlignmentSupportData {
      
            // data members
            std::string AllCharData;
            uint32_t    BlockLength;
            uint32_t    NumCigarOperations;
            uint32_t    QueryNameLength;
            uint32_t    QuerySequenceLength;
            bool        IsParsed;
            
            // constructor
            BamAlignmentSupportData(void)
                : BlockLength(0)
                , NumCigarOperations(0)
                , QueryNameLength(0)
                , QuerySequenceLength(0)
                , IsParsed(false)
            { }
        };
        
        BamTools::BamAlignment::BamAlignmentSupportData SupportData;  // Contains raw character data & lengths 

    // Alignment flag query constants
    // Use the get/set methods above instead
    private:
        enum { PAIRED        = 1
             , PROPER_PAIR   = 2
             , UNMAPPED      = 4
             , MATE_UNMAPPED = 8
             , REVERSE       = 16
             , MATE_REVERSE  = 32
             , READ_1        = 64
             , READ_2        = 128
             , SECONDARY     = 256
             , QC_FAILED     = 512
             , DUPLICATE     = 1024 
             };
};
}

%{

%}

