
window.BioJS = function(){
    // Translation code was from 
    // http://plindenbaum.blogspot.com/2010/12/server-side-javascript-translating-dna.html

    function TranslationTable(name,ncbiString) {
   	    this.name=name;
   	    this.ncbiString=ncbiString;
    };
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
    }

    var PROTEIN_TAGS = [
        ["FLAG","DYKDDDDK"],
        ["FLAG","DYKDHDI"],
        ["FLAG","DYKDHDG"],
        ["HA","YPYDVPDYA"],
        ["6xHIS","HHHHHH"],
        ["Myc","EQKLISEEDL"],
        ["TEV","ENLYFQG"],
        ["Myr","MGSNKSKPKDASQRR"],
        ["Myr","MGSSKSKPKDPSQRA"],
        ["V5","GKPIPNPLLGLDST"],
        ["S15","KETAAAKFERQHMDS"],
        ["Strep Tag","WSHPQFEK"]
    ];

    function __repeat(s,n) {
        var r=""; for (var a=0;a<n;a++) r+=s; return r;
    }

    function __normalize(c) {
        if (c == '*' ||
            c == 'D' ||
            c == 'H' ||
            c == 'M' ||
            c == 'N' ||
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
            c == 'n' ||
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
                line += seq_string[i];
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

    DNASequence.prototype.sequence=function() { return this.__sequence; }
    DNASequence.prototype.length=function() { return this.__sequence.length; }

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
    }

    // Returns substring
    DNASequence.prototype.substring=function(i,j){
        return new DNASequence(this.__sequence.substring(i,j));
    }

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
    }

    function ProteinSequence(seq_string) { this.__sequence = seq_string; }

    ProteinSequence.prototype.sequence=function() { return this.__sequence; }
    ProteinSequence.prototype.length=function() { return this.__sequence.length; }

    // Format protein to HTML, with bp markers. Also highlight tags
    // with giraffe-tag CSS class.
    ProteinSequence.prototype.format_html_with_bp=function() {
        // We have to construct three columns all in ONE SINGLE row,
        // left and right most columns for bp markers, and middle for
        // sequence. This allows 1) user selection (for copy/paste) of
        // just sequence in the middle w/o bp markers, and 2) using
        // span to highlight sequence fragments across lines.

        var line_width = 60;
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
            if (line == 0) { left_markers.push((i+1)); }
            else if (seg == 0) { res += '&nbsp;'; }

            if (!in_tag) {
                for (var j in PROTEIN_TAGS) {
                    var l = PROTEIN_TAGS[j][1].length;
                    if (this.__sequence.substr(i,l).toLowerCase() ==
                        PROTEIN_TAGS[j][1].toLowerCase()) {
                        in_tag = true;
                        tag_length = l;
                        res += '<span class="giraffe-tag" title="'+PROTEIN_TAGS[j][0]+'">';
                    }
                }
            }
            res += this.__sequence[i];
            line++;
            seg++;

            if (in_tag) { // found a new tag or was already in new tag
                tag_length--; 
                if (tag_length == 0) {
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
    }

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

    // Returns GenBank format of sequence. Currently does not handle
    // feature list, we will work on that.
    function genbank(seq_object,name,html) {
        if (html === undefined) { html = true; }

        var sp = "&nbsp;";
        var delim = "<br/>";
        var tab = sp+sp+sp+sp;
        if (!html) { sp = " "; delim = "\n"; tab = "\t"; }

        s =
            'LOCUS'+__repeat(sp,7)+name+tab+seq_object.length()+' bp '+
                    tab+'DNA'+tab+'SYN'+delim+
            'DEFINITION'+__repeat(sp,2)+name+delim+
            'ACCESSION'+__repeat(sp,3)+delim+
            'KEYWORDS'+__repeat(sp,4)+delim+
            'SOURCE'+__repeat(sp,6)+delim+
            __repeat(sp,2)+'ORGANISM'+__repeat(sp,2)+
            'other sequences; artificial sequences; vectors.'+delim+
            'ORIGIN'+delim;

        var lines_10 = wrap(seq_object.sequence(),10);
        for (var i=0,j=0; i<lines_10.length; i++) {
            if (j == 0) {
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
            ' method="post" target="_blank">'+
            '<input type="hidden" name="DATABASE" value="nr" />'+
            '<input type="hidden" name="PAGE" value="Nucleotides" />'+
            '<input type="hidden" name="PROGRAM" value="blastn" />'+
            '<input type="hidden" name="SERVICE" value="plain" />'+
            '<input type="hidden" name="GET_SEQUENCE" value="yes" />'+
            '<input type="hidden" name="QUERY" value="'+seq_object.sequence()+'">'+
            '<input type="submit" value="BLASTN">'+
            '</form>';
        return form;
    }

    function NCBI_blastx_form(seq_object) {
        var form =
            '<form action="http://www.ncbi.nlm.nih.gov/BLAST/Blast.cgi"'+
            ' method="post" target="_blank">'+
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
    }
}();


