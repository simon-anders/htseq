// ***************************************************************************
// bamtools_header.h (c) 2010 Derek Barnett, Erik Garrison
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 1 June 2010
// ---------------------------------------------------------------------------
// Prints the SAM-style header from a single BAM file ( or merged header from
// multiple BAM files) to stdout
// ***************************************************************************

#ifndef BAMTOOLS_HEADER_H
#define BAMTOOLS_HEADER_H

#include "bamtools_tool.h"

namespace BamTools {
  
class HeaderTool : public AbstractTool {
  
    public:
        HeaderTool(void);
        ~HeaderTool(void);
  
    public:
        int Help(void);
        int Run(int argc, char* argv[]); 
        
    private:
        struct HeaderSettings;
        HeaderSettings* m_settings;
};
  
} // namespace BamTools

#endif // BAMTOOLS_HEADER_H
