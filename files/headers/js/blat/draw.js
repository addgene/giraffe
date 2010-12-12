jQuery(document).ready(function($){
	// Creates canvas 320x200 at 10, 50
	var paper = Raphael(0, 0, 1000, 1000);
	// Creates circle at x = 50, y = 40, with radius 10
	var circle = paper.circle(400, 400, 200);
	// Sets the stroke attribute of the circle to white
	circle.attr("stroke", "#f00");
})
