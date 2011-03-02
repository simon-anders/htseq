import re
import HTSeq

_re_vcf_meta_comment = re.compile( "^##([a-zA-Z]+)\=(.*)$" )

_re_vcf_meta_descr = re.compile('ID=[^,]+,?|Number=[^,]+,?|Type=[^,]+,?|Description="[^"]+",?')

_re_vcf_meta_types = re.compile( "[INFO|FILTER|FORMAT]" )

_vcf_typemap = {
    "Integer":int,
    "Float":float,
    "String":str,
    "Flag":bool
}

class VariantCall( object ):
    
    def __init__( self, line, nsamples = 0, sampleids=[] ):
        if nsamples == 0:
            self.format = None
            self.chrom, self.pos, self.id, self.ref, self.alt, self.qual, self.filter, self.info = line.rstrip("\n").split("\t", 7)
        else:
            lsplit = line.rstrip("\n").split("\t")
            self.chrom, self.pos, self.id, self.ref, self.alt, self.qual, self.filter, self.info = lsplit[:8]
            self.format = lsplit[8].split(":")
            self.samples = {}
            spos=9
            for sid in sampleids:
                self.samples[ sid ] = dict( ( name, value ) for (name, value) in itertools.izip( self.format, lsplit[spos].split(":") ) )
                spos += 1
        self.pos = GenomicPosition( self.chrom, int(self.pos) )
        self.alt = self.alt.split(",")
    
    def __descr__( self ):
        return "<VariantCall at %s, ref '%s', alt %s >" % (str(self.pos).rstrip("/."), self.ref, str(self.alt).strip("[]"))
    
    def __str__( self ):
        return "%s:'%s'->%s" % (str(self.pos).rstrip("/."), self.ref, str(self.alt).strip("[]"))
    
    def unpack_info( self, infodict ):
        tmp = {}
        for token in self.info.split(";"):
            if re.compile("=").search(token):
                token = token.split("=")
                if infodict.has_key( token[0] ):
                    tmp[token[0]] = map( infodict[token[0]], token[1].split(",") )
                else:
                    tmp[token[0]] = token[1].split(",")
                if len( tmp[ token[0] ] ) == 1:
                    tmp[token[0]] = tmp[token[0]][0]
            else: #Flag attribute found
                tmp[token] = True
        diff = set( infodict.keys() ).difference( set( tmp.keys() ) )
        for key in diff:
            if infodict[key] == bool:
                tmp[key] = False
        self.info = tmp

class VCF_Reader( FileOrSequence ):

    def __init__( self, filename_or_sequence ):
        FileOrSequence.__init__( self, filename_or_sequence )
        self.metadata = {}
        self.info = {}
        self.filters = {}
        self.formats = {}
        self.nsamples = 0
        self.sampleids = []
        
    def make_info_dict( self ):
        self.infodict = dict( ( key, _vcf_typemap[self.info[key]["Type"]] ) for key in self.info.keys() )
    
    def parse_meta( self ):
        for line in FileOrSequence.__iter__( self ):
            if line.startswith( '#' ):
                if line.startswith( "##" ):
                    mo = _re_vcf_meta_comment.match( line )
                    if mo:
                        value = mo.group(2)
                        if mo.group(1) == "INFO":
                            value = dict( e.rstrip(",").split("=",1) for e in _re_vcf_meta_descr.findall(value) )
                            key = value["ID"]
                            del value["ID"]
                            self.info[ key ] = value
                        elif mo.group(1) == "FILTER":
                            value = dict( e.rstrip(",").split("=",1) for e in _re_vcf_meta_descr.findall(value) )
                            key = value["ID"]
                            del value["ID"]
                            self.filters[ key ] = value
                        elif mo.group(1) == "FORMAT":
                            value = dict( e.rstrip(",").split("=",1) for e in _re_vcf_meta_descr.findall(value) )
                            key = value["ID"]
                            del value["ID"]
                            self.formats[ key ] = value
                        else:
                            self.metadata[ mo.group(1) ] = mo.group(2)
                else:
                    self.sampleids = line.rstrip("\t\n").split("\t")[9:]
                    self.nsamples = len( self.sampleids )
                continue
            else:
                break
    
    def __iter__( self ):
        for line in FileOrSequence.__iter__( self ):
            if line == "\n" or line.startswith( '#' ):
                continue
            vc = VariantCall( line, self.nsamples, self.sampleids )
            yield vc


