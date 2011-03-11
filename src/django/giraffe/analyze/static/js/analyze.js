// Call:
//    GiraffeAnalyze(jQuery,gd,{...});
//  or
//    jQuery(document).ready(function($){
//        GiraffeAnalyze($,gd,{...});
//    })
//
// Options:
//    dom_id: dom ID to attach the analyzer window
// 

(function(){window.GiraffeAnalyze = function ($,gd,options) {
    var dom_id = 'giraffe-analyze';
    if ('dom_id' in options) { dom_id = options['dom_id']; }
    var name = 'Sequence';
    if ('name' in options) { name = options['name']; }

    function sequence_tab(dom) {
        var viewer = $('<div></div>');
        var fasta = $('<div></div>');
        var genbank = $('<div></div>');
        var rc = $('<div></div>');

        var link_viewer = $('<a href="#">Viewer</a>').click(function(){
            hide_all();
            $(this).addClass('giraffe-link-on');
            $(viewer).show();
        });
        var link_fasta = $('<a href="#">FASTA</a>').click(function(){
            hide_all();
            $(this).addClass('giraffe-link-on');
            $(fasta).show();
        });
        var link_genbank = $('<a href="#">GenBank</a>').click(function(){
            hide_all();
            $(this).addClass('giraffe-link-on');
            $(genbank).show();
        });
        var link_rc = $('<a href="#">Reverse Complement</a>').click(function(){
            hide_all();
            $(this).addClass('giraffe-link-on');
            $(rc).show();
        });

        var all_divs = [viewer,fasta,genbank,rc];
        var all_links = [link_viewer,link_fasta,link_genbank,link_rc];

        function hide_all() {
            for (var i in all_divs) { $(all_divs[i]).hide(); } 
            for (var i in all_links) { $(all_links[i]).removeClass('giraffe-link-on'); } 
        }

        var links = $('<p></p>')
            .append(link_viewer).append(' | ')
            .append(link_fasta).append(' | ')
            .append(link_genbank).append(' | ')
            .append(link_rc);

        $(dom)
            .append('<p>Current sequence: '+gd.sequence.length+' base pairs</p>')
            .append(links)
            .append(viewer)
            .append(fasta)
            .append(genbank)
            .append(rc);

        $(fasta).addClass('giraffe-seq').append(BioJS.fasta_html(gd.sequence,name));

        hide_all();
        $(link_fasta).addClass('giraffe-link-on');
        $(fasta).show();
    }

    function map_tab(dom) {
        var dom_id_c = 'giraffe-'+Math.floor(Math.random()*100000000);
        var dom_map = $('<div id="'+dom_id_c+'" class="giraffe-analyze-circular-map"></div>');
        $(dom).append(dom_map)
        gd.draw_circular_map({
            'map_dom_id' : dom_id_c,
            'plasmid_name' : name,
            'fade_time' : 300,
            'cutters': [1]
        });
    }

    function full_widget() {
        // Create each tab
        var dom_id_sequence = 'giraffe-'+Math.floor(Math.random()*100000000);
        var dom_tab_sequence = $('<div id="'+dom_id_sequence+'"></div>');
        var dom_id_map = 'giraffe-'+Math.floor(Math.random()*100000000);
        var dom_tab_map = $('<div id="'+dom_id_map+'"></div>');
        var dom_id_blast = 'giraffe-'+Math.floor(Math.random()*100000000);
        var dom_tab_blast = $('<div id="'+dom_id_blast+'"></div>');
        var dom_id_align = 'giraffe-'+Math.floor(Math.random()*100000000);
        var dom_tab_align = $('<div id="'+dom_id_align+'"></div>');
        var dom_id_digest = 'giraffe-'+Math.floor(Math.random()*100000000);
        var dom_tab_digest = $('<div id="'+dom_id_digest+'"></div>');
        var dom_id_translate = 'giraffe-'+Math.floor(Math.random()*100000000);
        var dom_tab_translate = $('<div id="'+dom_id_translate+'"></div>');
        // Main tab bar
        var dom_tabs = $('<div></div>').append(
            '<ul>'+
            '<li><a href="#'+dom_id_sequence+'">Sequence</a></li>'+
            '<li><a href="#'+dom_id_map+'">Map and Features</a></li>'+
            '<li><a href="#'+dom_id_blast+'">Blast</a></li>'+
            '<li><a href="#'+dom_id_align+'">Align</a></li>'+
            '<li><a href="#'+dom_id_digest+'">Digest</a></li>'+
            '<li><a href="#'+dom_id_translate+'">Translate</a></li>'+
            '</ul>'
        ).append(dom_tab_sequence)
         .append(dom_tab_map)
         .append(dom_tab_blast)
         .append(dom_tab_align)
         .append(dom_tab_digest)
         .append(dom_tab_translate);

        $('#'+dom_id).append(dom_tabs);
        $(dom_tabs).tabs();

        map_tab(dom_tab_map);
        sequence_tab(dom_tab_sequence);
    }

    full_widget();

}})();

