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
    var sequence = new BioJS.DNASequence(gd.sequence);

    function Switch_Panes(panes) {
        var divs = [];
        var links = [];

        function hide_all() {
            for (var i in divs) { $(divs[i]).hide(); } 
            for (var i in links) { $(links[i]).removeClass('giraffe-link-on'); } 
        }
        function show(i) {
            hide_all();
            $(links[i]).addClass('giraffe-link-on');
            $(divs[i]).show();
        }
        function link(i) { return links[i]; }
        function pane(i) { return divs[i]; }

        var links_dom = $('<p></p>');
        var divs_dom = $('<div></div>');

        for (var i in panes) {
            divs[i] = $('<div></div>');
            links[i] = $('<a pane="'+i+'" href="#">'+panes[i]+'</a>').click(function() {
                var i = $(this).attr('pane'); show(i);
            });
            $(links_dom).append(links[i]);
            if (i < panes.length-1) { $(links_dom).append(' | '); }
            $(divs_dom).append(divs[i]);
        }

        return {
            'links' : links_dom,
            'panes' : divs_dom,
            'pane' : pane,
            'link' : link,
            'show' : show,
            'hide_all' : hide_all,
        }
    }

    function sequence_tab(dom) {
        panes = Switch_Panes(['Viewer', 'Fasta', 'GenBank', 'Reverse Complement']);

        $(dom).append('<p>Current sequence: '+sequence.length()+' base pairs</p>')
              .append(panes.links)
              .append(panes.panes);

        $(panes.pane(1))
            .addClass('giraffe-seq').append(BioJS.fasta(sequence,name));
        $(panes.pane(2))
            .addClass('giraffe-seq').append(BioJS.genbank(sequence,name));
        $(panes.pane(3))
            .addClass('giraffe-seq').append(
                BioJS.fasta(sequence.reverse_complement(),name)
            );
        panes.hide_all();
        panes.show(1);
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

    function translate_tab(dom) {
        panes = Switch_Panes(
            ['ORFs',
             "5'3' Frame 1",
             "5'3' Frame 2",
             "5'3' Frame 3",
             "3'5' Frame 1",
             "3'5' Frame 2",
             "3'5' Frame 3"]
        );

        $(dom).append(panes.links)
              .append(panes.panes);

        $(panes.pane(1))
            .addClass('giraffe-seq').append(sequence.translate().format_html());
        $(panes.pane(2))
            .addClass('giraffe-seq').append(sequence.substring(1).translate().format_html());
        $(panes.pane(3))
            .addClass('giraffe-seq').append(sequence.substring(2).translate().format_html());
        $(panes.pane(4))
            .addClass('giraffe-seq').append(
                sequence.reverse_complement().translate().format_html());
        $(panes.pane(5))
            .addClass('giraffe-seq').append(
                sequence.reverse_complement().substring(1).translate().format_html());
        $(panes.pane(6))
            .addClass('giraffe-seq').append(
                sequence.reverse_complement().substring(2).translate().format_html());

        panes.hide_all();
        panes.show(1);
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
        translate_tab(dom_tab_translate);

    }

    full_widget();

}})();

