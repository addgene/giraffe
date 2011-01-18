function blat_api_closure(features) {
    var $ = jQuery
	
    function draw() {
        alert(features.length)
	}

    function panel() {
        $('#blat_panel').addClass('ui-tabs')
          .append(
            '<ul class="ui-tabs-nav">'+
            '<li><a href="#tabs-1">A</a></li>'+
            '<li><a href="#tabs-2">B</a></li>'+
            '<li><a href="#tabs-3">C</a></li>'+
            '</ul>'
          )
          .append(
            '<div id="tabs-1" class="ui-tabs-panel">Tabs A'+
                '<div id="accordian">'+
                    '<h3><a href="#">Section 1</a></h3>'+
                    '<div>Section 1</div>'+
                    '<h3><a href="#">Section 2</a></h3>'+
                    '<div>Section 2</div>'+
                    '<h3><a href="#">Section 3</a></h3>'+
                    '<div>Section 3</div>'+
                '</div>'+
                '</div>'+
            '<div id="tabs-2" class="ui-tabs-panel ui-tabs-hide">Tabs B</div>'+
            '<div id="tabs-3" class="ui-tabs-panel ui-tabs-hide">Tabs C</div>'
          ).tabs()
       $('#accordian').accordion()
    }

    function about() { alert('BLAT API') }

    return {
        'draw' : draw,
        'panel' : panel,
        'about' : about
    }
}

function blat_api(features) {
    // global assignment
    blat_closure = blat_api_closure(features)
    return blat_closure
}

