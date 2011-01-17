function blat_api_closure(features) {
	
    function draw() {
        alert(features.length)
	}
    function about() { alert('BLAT API') }

    return {
        'draw' : draw,
        'about' : about
    }
}

function blat_api(features) {
    // global assignment
    blat_closure = blat_api_closure(features)
    return blat_closure
}

