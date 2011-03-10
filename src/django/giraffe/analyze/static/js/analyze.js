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

(function(){window.GiraffeAnalyze = function ($,gd_object,options) {
    var dom_id = 'giraffe-analyze';
    if ('dom_id' in options) { dom_id = options['dom_id']; }

    // Layout the big widget

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

    var dom_id_c = 'giraffe-'+Math.floor(Math.random()*100000000);
    var dom_map = $('<div id="'+dom_id_c+'" class="giraffe-analyze-circular-map"></div>');
    $(dom_tab_map).append(dom_map)
    gd.draw_circular_map({
        'map_dom_id' : dom_id_c,
        'plasmid_name' : 'Test',
        'fade_time' : 300,
        'cutters': [1]
    });

    $(dom_tab_sequence).append(gd.sequence);

}})();

