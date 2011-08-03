window.BioJS = function(){
    if (typeof(Array.map) !== 'function') {

        /**
         * Mapping functions onto arrays
         *
         * @param   {Function} callback the callback function
         *                              given array item as first argument
         * @param   {Object}   ctxt    [optional] object to use as
         *                             this pointer for callback call
         * @returns {Array}             a new array with the results of
         *                              callback mapped to every element of this
         * @ignore
         */
        Array.prototype.map = function (callback, ctxt) {
            var ix;
            var context;
            var mapped = new Array(this.length);

            for (ix = 0; ix < this.length; ix++) {
                // If second agrument was provided, use the right context
                if (arguments.length > 1) {
                    context = ctxt;
                } else {
                    context = this[ix];
                }

                mapped[ix] = callback.call(context, this[ix]);

            }

            return mapped;
        };
    }
    // Translation code was from
    // http://plindenbaum.blogspot.com/2010/12/server-side-javascript-translating-dna.html

    function TranslationTable(name,ncbiString) {
        this.name=name;
        this.ncbiString=ncbiString;
    }

    TranslationTable.STANDARD=new TranslationTable(
        "standard",
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    );
    TranslationTable.prototype.base2index=function(c) {
        switch(c) {
            case 'T':case 't': return 0;
            case 'C':case 'c': return 1;
            case 'A':case 'a': return 2;
            case 'G':case 'g': return 3;
            default: return -1;
        }
    };
    TranslationTable.prototype.translate=function(b1, b2, b3) {
        if (arguments.length==1) {
            var prot="";
            for(var i=0; i+2<b1.length; i+=3) {
                prot+= this.translate(b1[i],b1[i+1],b1[i+2]);
            }
            return prot;
        }
        var base1= this.base2index(b1);
        var base2= this.base2index(b2);
        var base3= this.base2index(b3);
        if (base1==-1 || base2==-1 || base3==-1) { return '?'; }
        else { return this.ncbiString[base1*16+base2*4+base3]; }
    };

    var PROTEIN_TAGS = [
        ["FLAG","DYKDDDDK"],
        ["FLAG","DYKDHDI"],
        ["FLAG","DYKDHDG"],
        ["VSVG Tag","YTDIEMNRLGK"],
        ["SV40 NLS","PKKKRKV"],
        ["SV40 NLS","PKKKRKVG"],
        ["T7 Tag", "MASMTGGQQMG"],
        ["NLS", "KKRKV"],
        ["HA","YPYDVPDYA"],
        ["6xHIS","HHHHHH"],
        ["Myc","EQKLISEEDL"],
        ["TEV","ENLYFQG"],
        ["TEV","ENPYFQG"],
        ["Myr","MGSNKSKPKDASQRR"],
        ["Myr","MGSSKSKPKDPSQRA"],
        ["V5","GKPIPNPLLGLDST"],
        ["S15","KETAAAKFERQHMDS"],
        ["Strep Tag","WSHPQFEK"]
    ];

    function __repeat(s,n) {
        var r="";
        for (var a=0;a<n;a++) {
            r+=s;
        }
        return r;
    }

    function __normalize(c) {
        if (c == '*' ||
            c == 'D' ||
            c == 'H' ||
            c == 'M' ||
            c == 'R' ||
            c == 'V' ||
            c == 'W') { return 'A'; }
        if (c == 'B' ||
            c == 'Y' ||
            c == 'S') { return 'C'; }
        if (c == 'K') { return 'G'; }
        if (c == 'U') { return 'T'; }
        if (c == 'd' ||
            c == 'h' ||
            c == 'm' ||
            c == 'r' ||
            c == 'v' ||
            c == 'w') { return 'a'; }
        if (c == 'b' ||
            c == 'y' ||
            c == 's') { return 'c'; }
        if (c == 'k') { return 'g'; }
        if (c == 'u') { return 't'; }
        return c;
    }

    // Returns array of sequence, split from input sequence
    function wrap(seq_string,width) {
        if (seq_string) {
            var res = [];
            var line = '';
            for (var i=0; i<seq_string.length; i++) {
                line += seq_string.charAt(i);
                if (line.length >= width){
                    res.push(line);
                    line = '';
                }
            }
            if (line.length > 0) { res.push(line); }
            return res;
        }
        return [];
    }

    function DNASequence(seq_string) { this.__sequence = seq_string; }

    DNASequence.prototype.sequence=function() { return this.__sequence; };
    DNASequence.prototype.length=function() { return this.__sequence.length; };

    DNASequence.prototype.reverse_complement=function(){
        if (this.__reverse_complement) { return this.__reverse_complement; }
        var r = this.__sequence.split("").reverse();
        var rc = r.map(function(c) {
            c = __normalize(c);
            if (c == 'a') { return 't'; }
            if (c == 'A') { return 'T'; }
            if (c == 'g') { return 'c'; }
            if (c == 'G') { return 'C'; }
            if (c == 't') { return 'a'; }
            if (c == 'T') { return 'A'; }
            if (c == 'c') { return 'g'; }
            if (c == 'C') { return 'G'; }
            if (c == 'n') { return 'n'; }
            if (c == 'N') { return 'N'; }
            return c;
        });
        rc = rc.join("");
        this.__reverse_complement = new DNASequence(rc);
        return this.__reverse_complement;
    };

    // Returns DNA to protein translation of sequence
    DNASequence.prototype.translate=function(){
        if (this.__translation) { return this.__translation; }
        var p = TranslationTable.STANDARD.translate(this.__sequence);
        this.__translation = new ProteinSequence(p);
        return this.__translation;
    };

    // Returns 1-indexed bp position of the first occurance of the
    // query sequence, or -1. start, also 1-indexed, can be the start
    // position to search, or -1.
    DNASequence.prototype.find=function(query,start){
        if (query === undefined || query.length === 0) { return -1; }
        if (this.__sequence_lc === undefined) {
            this.__sequence_lc = this.__sequence.toLowerCase();
        }
        var n;
        if (start > 0) {
            n = this.__sequence_lc.indexOf(query.toLowerCase(),start-1);
        }
        else {
            n = this.__sequence_lc.indexOf(query.toLowerCase());
        }
        if (n >= 0) { return n+1; }
        return n;
    };

    // Returns substring
    DNASequence.prototype.substring=function(i,j){
        var a;
        if (j === undefined) { a = this.__sequence.substring(i); }
        else { a = this.__sequence.substring(i,j); }
        return new DNASequence(a);
    };

    // Format DNA sequence to HTML
    DNASequence.prototype.format_html=function() {
        var line_width = 80;
        if (this.__html) { return this.__html; }
        var lines_vec = wrap(this.__sequence,line_width);
        s = '';
        for (var i=0; i<lines_vec.length; i++) {
            s += lines_vec[i]+'<br/>';
        }
        this.__html = s;
        return this.__html;
    };

    // Format DNA sequence, with amino acid sequence overlay.
    // seq_start and seq_end should be 1-indexed base pair numbers in
    // 5' to 3' direction of *this* sequence, so if are displaying
    // reverse complement, then end should be smaller than start.
    DNASequence.prototype.format_html_with_aa=function(seq_start,seq_end) {
        if (this.__aa_html) { return this.__aa_html; }

        if (seq_start === undefined) {
            seq_start = 1;
        }
        if (seq_end === undefined) {
            seq_end = this.__sequence.length+seq_start-1;
        }
        var aa = this.translate();

        // line_width MUST BE multiple of 3
        var line_width = 60;
        var dna_vec = wrap(this.__sequence,line_width);
        var aa_vec = wrap(aa.sequence(),Math.floor(line_width/3));

        var left_markers = [];
        var right_markers = [];
        s = '';
        for (var i=0; i<dna_vec.length; i++) {
            if (i<aa_vec.length) {
                if (i+1 == dna_vec.length) {
                    left_markers.push(i*(line_width/3)+1);
                    right_markers.push(aa.length());
                    if (seq_start < seq_end) {
                        left_markers.push(seq_start+i*line_width);
                        right_markers.push(seq_end);
                    }
                    else {
                        left_markers.push(seq_start-i*line_width);
                        right_markers.push(seq_end);
                    }
                }
                else {
                    left_markers.push(i*(line_width/3)+1);
                    right_markers.push((i+1)*(line_width/3));
                    if (seq_start < seq_end) {
                        left_markers.push(seq_start+i*line_width);
                        right_markers.push(seq_start+(i+1)*line_width-1);
                    }
                    else {
                        left_markers.push(seq_start-i*line_width);
                        right_markers.push(seq_start-(i+1)*line_width+1);
                    }
                }
                var p = aa_vec[i];
                var l = p.split('').join('&nbsp;&nbsp;');
                s += l+'<br/>';
            }
            s += '<span class="giraffe-seq-overlay">'+dna_vec[i]+'</span><br/>';
        }

        var table = '<table><tr>'+
            '<td class="giraffe-bp-marker giraffe-bp-marker-left">'+
            left_markers.join('<br/>')+
            '</td><td>'+s+'</td>'+
            '<td class="giraffe-bp-marker giraffe-bp-marker-right">'+
            right_markers.join('<br/>')+
            '</td></tr></table>';

        this.__aa_html = table;
        return this.__aa_html;
    };

    function ProteinSequence(seq_string) { this.__sequence = seq_string; }

    ProteinSequence.prototype.sequence=function() { return this.__sequence; };
    ProteinSequence.prototype.length=function() { return this.__sequence.length; };

    // Format protein to HTML, with bp markers. Also highlight tags
    // with giraffe-tag CSS class.
    ProteinSequence.prototype.format_html_with_bp=function() {
        // We have to construct three columns all in ONE SINGLE row,
        // left and right most columns for bp markers, and middle for
        // sequence. This allows 1) user selection (for copy/paste) of
        // just sequence in the middle w/o bp markers, and 2) using
        // span to highlight sequence fragments across lines.

        var line_width = 50;
        var seg_width = 10;
        if (this.__html) { return this.__html; }
        if (this.__sequence === undefined) { return ""; }

        var left_markers = [];
        var right_markers = [];

        var res = '';
        var in_tag = false;
        var tag_length = 0;
        var line = 0;
        var seg = 0;
        for (var i=0; i<this.__sequence.length; i++) {
            if (line === 0) {
                left_markers.push((i+1));
            } else if (seg === 0) {
                res += '&nbsp;';
            }

            if (!in_tag) {
                for (var j = 0; j < PROTEIN_TAGS.length; j++) {
                    var l = PROTEIN_TAGS[j][1].length;
                    if (this.__sequence.substr(i,l).toLowerCase() ==
                        PROTEIN_TAGS[j][1].toLowerCase()) {
                        in_tag = true;
                        tag_length = l;
                        res += '<span class="giraffe-tag" title="'+PROTEIN_TAGS[j][0]+'">';
                    }
                }
            }
            res += this.__sequence.charAt(i);
            line++;
            seg++;

            if (in_tag) { // found a new tag or was already in new tag
                tag_length--;
                if (tag_length === 0) {
                    in_tag = false;
                    res += '</span>';
                }
            }
            if (line == line_width) {
                res += '<br/>';
                right_markers.push ((i+1));
                line = 0;
                seg = 0;
            }
            else if (seg == seg_width) { seg = 0; }
        }

        // just to be sure, but we really should never be here...
        if (in_tag) { res += '</span>'; }

        if (line > 0) { res += '<br/>'; }

        var table = '<table><tr>'+
            '<td class="giraffe-bp-marker giraffe-bp-marker-left">'+
            left_markers.join('<br/>')+
            '</td><td>'+res+'</td>'+
            '<td class="giraffe-bp-marker giraffe-bp-marker-right">'+
            right_markers.join('<br/>')+
            '</td></tr></table>';

        this.__html = table;
        return this.__html;
    };

    // Returns FASTA format of sequence
    function fasta(seq_object,name,html) {
        var line_width = 80;
        if (html === undefined) { html = true; }
        var lines_vec = wrap(seq_object.sequence(),line_width);
        var s = '';
        if (html) { s = '<p>&gt;'+name+'<br/>'; }
        else { s = '>'+name+'\n'; }
        for (var i=0; i<lines_vec.length; i++) {
            s += lines_vec[i];
            if (html) { s += '<br/>'; } else { s += '\n'; }
        }
        if (html) { s += '</p>'; }
        return s;
    }

    // Returns GenBank format of sequence.
    //
    // features array should be an array of objects. Each object
    // should have these keys:
    //    label
    //    type -- must be a valid GenBank feature type
    //    start,end -- always in 5' to 3', we will adjust for clockwise
    //    clockwise -- true or false
    //    clockwise_sequence -- optional, but if type is CDS, then you
    //                          need to supply this if you want to see
    //                          a translation shown.
    //    gene -- optional
    //    notes -- optional
    //
    function genbank(seq_object,name,html,features) {
        if (html === undefined) { html = true; }
        if (features === undefined) { features = []; }

        var i, j;

        var sp = "&nbsp;";
        var delim = "<br/>";
        var tab = sp+sp+sp+sp;
        if (!html) { sp = " "; delim = "\n"; tab = "\t"; }

        s = 'LOCUS'+__repeat(sp,7)+name+tab+seq_object.length()+' bp '+
                    tab+'DNA'+tab+'SYN'+delim+
            'DEFINITION'+__repeat(sp,2)+name+delim+
            'ACCESSION'+__repeat(sp,3)+delim+
            'KEYWORDS'+__repeat(sp,4)+delim+
            'SOURCE'+__repeat(sp,6)+delim+
            __repeat(sp,2)+'ORGANISM'+__repeat(sp,2)+
            'other sequences; artificial sequences; vectors.'+delim;

        function __gb(p) {
            // need special formatting: lw=58, first=44
            var res = [];
            var line = '';
            var ctr = 0;
            for (var i=0; i<p.length; i++) {
                line += p.charAt(i);
                ctr++;
                if ((res.length === 0 && ctr >= 44) || ctr >= 58) {
                    res.push(line);
                    line = __repeat(sp,21);
                    ctr = 0;
                }
            }
            if (ctr > 0) { res.push(line); }
            return res.join(delim);
        }

        if (features.length) {
            s += "FEATURES"+__repeat(sp,13) +
                "Location/Qualifiers"+delim+__repeat(sp,5)+"source";
            s += __repeat(sp, 16-"source".length);
            s += "1.."+seq_object.length()+delim+__repeat(sp,21) +
                "/organism=\""+name+"\""+delim +
                __repeat(sp,21)+"/mol_type=\"other DNA\""+delim;
        }

        for (i = 0; i < features.length; i++) {
            var type = features[i].type;
            var tran;
            if (type == 'CDS' && features[i].clockwise_sequence &&
                features[i].clockwise_sequence !== '') {
                if (features[i].clockwise) {
                    tran = new BioJS.DNASequence(features[i].clockwise_sequence).translate();
                }
                else {
                    tran = new BioJS.DNASequence(features[i].clockwise_sequence)
                        .reverse_complement().translate();
                }
            }
            s += __repeat(sp,5)+type;
            var nsp = 16 - type.length;
            if (nsp < 0) {
                nsp = 0;
            }
            s += __repeat(sp, nsp);
            if (features[i].clockwise) {
                s += features[i].start+".."+features[i].end+delim;
            }
            else {
                s += "complement("+features[i].start+".."+features[i].end+")"+delim;
            }
            s += __repeat(sp,21)+"/label=\""+features[i].label+"\""+delim;
            if (features[i].gene && features[i].gene !== '') {
                s += __repeat(sp,21)+"/gene=\""+features[i].gene+"\""+delim;
            }
            if (features[i].notes && features[i].notes !== '') {
                s += __repeat(sp,21)+"/note=\""+features[i].notes+"\""+delim;
            }
            // User annotated CDS may not be in-frame...
            if (tran) {
                s += __repeat(sp,21)+"/translation=\""+__gb(tran.sequence())+"\""+delim;
            }
        }

        s += 'ORIGIN'+delim;

        var lines_10 = wrap(seq_object.sequence(),10);
        for (i=0,j=0; i<lines_10.length; i++) {
            if (j === 0) {
                var start = i*10+1;
                if (start < 10) { s += __repeat(sp,4); }
                else if (start < 100) { s += __repeat(sp,3); }
                else if (start < 1000) { s += __repeat(sp,2); }
                else if (start < 10000) { s += __repeat(sp,1); }
                s += start;
            }
            s += sp+lines_10[i];
            j++;
            if (j == 6) { s += delim; j = 0; }
        }
        s += delim+'//'+delim;
        return s;
    }

    function NCBI_blastn_form(seq_object) {
        var form =
            '<form action="http://www.ncbi.nlm.nih.gov/BLAST/Blast.cgi"'+
            ' method="POST" target="_blank">'+
            '<input type="hidden" name="DATABASE" value="nr" />'+
            '<input type="hidden" name="PAGE" value="Nucleotides" />'+
            '<input type="hidden" name="PROGRAM" value="blastn" />'+
            '<input type="hidden" name="SERVICE" value="plain" />'+
            '<input type="hidden" name="GET_SEQUENCE" value="yes" />'+
            '<input type="hidden" name="QUERY" value="'+seq_object.sequence()+'" />'+
            '<input type="submit" value="BLASTN" />'+
            '</form>';
        return form;
    }

    function NCBI_blastx_form(seq_object) {
        var form =
            '<form action="http://www.ncbi.nlm.nih.gov/BLAST/Blast.cgi"'+
            ' method="POST" target="_blank">'+
            '<input type="hidden" name="DATABASE" value="nr" />'+
            '<input type="hidden" name="PAGE" value="Nucleotides" />'+
            '<input type="hidden" name="PROGRAM" value="blastx" />'+
            '<input type="hidden" name="SERVICE" value="plain" />'+
            '<input type="hidden" name="GET_SEQUENCE" value="yes" />'+
            '<input type="hidden" name="QUERY" value="'+seq_object.sequence()+'" />'+
            '<input type="submit" value="BLASTX" />'+
            '</form>';
        return form;
    }

    function NCBI_blast2_form(seq_object) {
        var form =
            '<form action="http://blast.ncbi.nlm.nih.gov/Blast.cgi"'+
            ' method="post" target="_blank">'+
            '<input type="hidden" name="PAGE" value=blastn />'+
            '<input type="hidden" name="PROGRAM" value=blastn />'+
            '<input type="hidden" name="BLAST_PROGRAMS" value=blastn />'+
            '<input type="hidden" name="PAGE_TYPE" value=BlastSearch />'+
            '<input type="hidden" name="BLAST_SPEC" value=blast2seq />'+
            '<input type="hidden" name="QUERY" value="'+seq_object.sequence()+'" />'+
            '<textarea name="SUBJECTS"></textarea>'+
            '<input type="submit" value="BLAST2" />'+
            '</form>';
        return form;
    }

    function NCBI_blastp_form(seq_object) {
        var form =
            '<form action="http://www.ncbi.nlm.nih.gov/BLAST/Blast.cgi"'+
            ' method="post" target="_blank">'+
            '<input type="hidden" name="DATABASE" value="nr" />'+
            '<input type="hidden" name="PAGE" value="Protein" />'+
            '<input type="hidden" name="PROGRAM" value="blastp" />'+
            '<input type="hidden" name="SERVICE" value="plain" />'+
            '<input type="hidden" name="GET_SEQUENCE" value="yes" />'+
            '<input type="hidden" name="QUERY" value="'+seq_object.sequence()+'" />'+
            '<input type="submit" value="BLASTP" />'+
            '</form>';
        return form;
    }

    function NCBI_recent_results_link() {
        return '<a href="http://www.ncbi.nlm.nih.gov/BLAST/Blast.cgi?CMD=GetSaved&RECENT_RESULTS=on" target="_blank"></a>';
    }

    return {
        'fasta' : fasta,
        'genbank' : genbank,
        'wrap' : wrap,
        'DNASequence' : DNASequence,
        'ProteinSequence' : ProteinSequence,
        'NCBI_blastn_form' : NCBI_blastn_form,
        'NCBI_blastx_form' : NCBI_blastx_form,
        'NCBI_blast2_form' : NCBI_blast2_form,
        'NCBI_blastp_form' : NCBI_blastp_form,
        'NCBI_recent_results_link' : NCBI_recent_results_link
    };
}();

// vi: set expandtab:ts=4:sw=4:sts=4
