
window.BioJS = function(){

    // Returns reverse complement of sequence
    function reverse_complement(seq) {
    }

    // Returns FASTA format of sequence
    function fasta_html(seq,name) {
        var lines_80 = wrap(seq,80);
        var s = '<p>&gt;'+name+'<br/>';
        for (var i=0; i<lines_80.length; i++) {
            s += lines_80[i]+'<br/>';
        }
        s += '</p>';
        return s;
    }

    // Returns GenBank format of sequence
    function genbank_html(seq) {
    }

    // Returns DNA to protein translation of sequence
    function translate(seq) {
    }

    // Returns array of sequence, split from input sequence
    function wrap(seq,width) {
        var res = [];
        var line = '';
        for (var i=0; i<seq.length; i++) {
            line += seq[i];
            if (line.length >= width){
                res.push(line);
                line = '';
            }
        }
        if (line.length > 0) { res.push(line); }
        return res;
    }

    return {
        'reverse_complement' : reverse_complement,
        'fasta_html' : fasta_html,
        'genbank_html' : genbank_html,
        'translate' : translate,
        'wrap' : wrap
    }
}();


