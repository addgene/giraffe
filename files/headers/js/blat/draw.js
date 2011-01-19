jQuery(document).ready(function ($) {
	// Global parameters
	// Paper setup
	var map_width = 640;
	var map_height = 480;
	var cx = map_width/2;
	var cy = map_height/2;

	// Global plasmid info
	var seq_length = parseInt($("#sequence-data-length").html());

	// Loop radii
	var radius_spacing = 20; // spacing
	var plasmid_radius = 200;
	var inner_radius = plasmid_radius - radius_spacing; 
	var outer_radius = plasmid_radius + radius_spacing;

	// Feature visual properties
	var feature_width = 15;
	var enzyme_width = 20;
	var feature_opacity = 1.0/3.0;
	var head_width = 25;
	var head_length = 7;

	// Overlaps
	var min_overlap_pct = 0.05;
	var min_overlap_feature_size = 10;
	
	// Tic marks
	var tic_mark_length = 15;
	var tic_mark_radius = inner_radius - tic_mark_length/2;
	var tic_label_radius = tic_mark_radius - 1.5*tic_mark_length;

	// Colors
	var color_bg_text = "#aaa";
	var color_plasmid = "#000";
	var color_feature = "#f00";
	var color_primer  = "#090";
	var color_origin  = "#333";
	var color_enzyme  = "#009";

	// Feature Types
	var ft_gene = "Gene";
	var ft_regulatory = "Regulatory";
	var ft_enzyme = "Enzyme";
	var ft_primer = "Primer";
	var ft_promoter = "Promoter";
	var ft_terminator = "Terminator";
	var ft_origin = "Origin";
	var ft_feature = "Feature";
	var ft_exact_feature = "Exact Feature";

	// SVG Object
	// For dealing with the SVG syntax
	var svg = {
		move: function (x, y) {
			return this.to_path(['M', x, y]);
		},
		arc: function (r, x, y) {
			return this.to_path(['A', r, r, 0, 0, 0, x, y]);
		},
		line: function (x, y) {
			return this.to_path(['L', x, y]);
		},
		close: function () {
			return 'z';
		},
		to_path: function (plist) {
			return plist[0] + plist.slice(1).join(' ');
		}
	};

	// Conversion Object
	// Groups conversions into handy section
	var convert = {
		pos_to_angle: function (p) {
			//     start at the top of the circle
			return 90 - (p/seq_length) * 360;
		},
		angle_to_pos: function (a) {
			//     start at the top of the circle
			return Math.round((seq_length/360) * ((360 + 90 - a) % 360));
		},
		/* Angle a is in degrees from the horizontal, counterclockwise, and
		 * r is relative to the center of the paper */
		polar_to_rect: function (r, a) {
			var rect = {};
			rect.x = cx + r * Math.cos(Raphael.rad(a));
			// Coordinates increase as you go down, so y is flipped.
			rect.y = cy - r * Math.sin(Raphael.rad(a));
			return rect;
		}
	};


	// Feature class
	function Feature(feature_list) {
		var name = feature_list[0];
		var start = parseInt(feature_list[1]);
		var end = parseInt(feature_list[2]);
		var type = feature_list[3];
		var clockwise = (start <= end);

		// Ensure that start and end go clockwise, regardless
		// of what the feature does
		if (!clockwise) {
			var tmp = start;
			start = end;
			end = tmp;
		}

		// Radius is public, unlike other properties, which are permanent
		this.radius = plasmid_radius; // Default to plasmid radius, can be changed
		                              // later on by other methods

		// Accessors for private properties set at creation
		this.name = function() { return name; };
		this.start = function() { return start; };
		this.end = function() { return end; };
		this.type = function() { return type; };
		this.clockwise = function() { return clockwise; };

		// Feature drawing
		var f = this; // JavaScript bug workaround
		this.draw = function () {

			// Property selection
			var color = color_feature;
			var width = feature_width;
			var draw_head = false;
			switch(type) {
				case ft_promoter:
					draw_head = true; // Promotors are the only primer-colored
				case ft_primer:       // features with heads
				case ft_terminator:
					color = color_primer;
					break;
				case ft_regulatory:
				case ft_origin:
					color = color_origin;
					break;
				case ft_enzyme:
					color = color_enzyme;
					width = enzyme_width;
					break;
				case ft_gene:
					draw_head = true;
					break;
			}

			// Convert from sequence positions to angles
			var a0 = convert.pos_to_angle(start);
			var a1 = convert.pos_to_angle(end);

			// Create the draw feature, a set which will have the head 
			// and arc pushed onto it as necessary.
			var draw_feature = paper.set();
			
			// Arrowhead drawing, if needed
			if (draw_head) {

				// Arrow tip point lines up with a0 or a1,
				// We need to figure out how many radians the arrow takes up
				// in order to adjust a0 or a1 by that amount, and to set the 
				// base of the triangle even with that angle
				var r_p = Math.sqrt(f.radius*f.radius + 
						head_length*head_length);
				// "height" of the arrowhead, in degrees
				var a_b;
				var a_p = Raphael.deg(Math.asin(head_length/r_p));
				// Adjust the appropriate edge to compensate for the arrowhead
				if (clockwise) {
					a_b = (a1 + a_p) % 360 ; // base angle
					a_p = a1;       // point angle
					a1  = a_b;      // adjust arc edge
				} else {
					a_b = (a0 - a_p) % 360 ; // base angle
					a_p = a0;       // point angle
					a0  = a_b;      // adjust arc edge
				}
				var xy_p = convert.polar_to_rect(f.radius, a_p);
				
				// bottom and top points, rectangular
				var xy_b = convert.polar_to_rect(f.radius - head_width/2.0, a_b);
				var xy_t = convert.polar_to_rect(f.radius + head_width/2.0, a_b);

				// Unlike the arc, the head is traced with a line, and
				// then created entirely with the fill color
				var head = paper.path(svg.move(xy_p.x, xy_p.y) +
									  svg.line(xy_b.x, xy_b.y) +
									  svg.line(xy_t.x, xy_t.y) + 
									  svg.close());
				head.attr({"stroke-width": 0,
						   "fill":         color});
				draw_feature.push(head);
			}

			// Arc drawing
			if (a1 < a0) { // Compensating for the head may have "taken up" all
				           // the room on the plasmid, in which case no arc needs
						   // to be drawn

				// Rectangular coordinates of the edges of the arc: 
				// arcs are drawn counterclockwise, even though the plasmid
				// sequence increases clockwise, so we flip the
				// indices
				var xy0 = convert.polar_to_rect(f.radius, a1);
				var xy1 = convert.polar_to_rect(f.radius, a0);

				// The arc has no fill-color: it's just a thick line
				var arc = paper.path(svg.move(xy0.x, xy0.y) +
									 svg.arc(f.radius, xy1.x, xy1.y));
				arc.attr({"stroke-width": width});

				draw_feature.push(arc);
			}

			// Apply the feature-wise properties to the whole feature
			draw_feature.attr({"stroke":         color,
			                   "stroke-linecap": "butt",
			                   "opacity":        feature_opacity,
			                   "title":          name});

			// Toggle solid/light upon click
			draw_feature.click(function (event) {
				if (this.attr("opacity") < 1.0) {
					draw_feature.animate({"opacity": 1}, 300);
				} else {
					draw_feature.animate({"opacity": feature_opacity}, 300);
				}
			});

		}
	}

	// Circle setup
	function draw_plasmid() {
		function draw_tic_mark(a) {
			var r0 = tic_mark_radius - tic_mark_length/2;
			var r1 = tic_mark_radius + tic_mark_length/2;
			var xy0 = convert.polar_to_rect(r0,a);
			var xy1 = convert.polar_to_rect(r1,a);
			var tic = paper.path(svg.move(xy0.x, xy0.y) +
								 svg.line(xy1.x, xy1.y));
			tic.attr({"stroke": color_bg_text});

			var xyl = convert.polar_to_rect(tic_label_radius, a);
			var label = paper.text(xyl.x, xyl.y, String(convert.angle_to_pos(a)));
			label.attr({"fill": color_bg_text});
			if (a < 90 || a > 270) { // Right half of wheel: align right
				label.attr({"text-anchor": "end"});
			} else if (a > 90 && a < 270) { // Left half of wheel: align left
				label.attr({"text-anchor": "start"});
			} // Top and bottom default to middle, which is correct
		}

		var plasmid = paper.circle(cx, cy, plasmid_radius);
		plasmid.attr("stroke", color_plasmid);
		var plasmid_label = paper.text(cx, cy, seq_length + " bp");
		plasmid_label.attr({"fill":      color_plasmid,
							"font-size": "12pt"});

		for (ang = 0; ang < 360; ang += 30) {
			draw_tic_mark(ang);
		}
	}

	// Load info for features
	function parse_features() {
		features = [];

		// Parse through the plasmid table and draw the features
		$("table.table-list>tbody>tr").each(function (rx) {
			row = [];
			$(this).find("td").each(function (fx) {
				foo = $(this).html();
				row.push($(this).html());
			})
			var feat = new Feature(row);
			features.push(feat);
		})

		return features;
	}
	
	// Move features that overlap to other radii.
	function resolve_conflicts(features) {
		var conflicts;
		var rad = plasmid_radius; // current radius
		var rx = 1;               // radius counter
		do {
			var new_rad = rad + Math.pow(-1, rx) * rx * radius_spacing;

			conflicts = 0; // Assume you have no conflicts until you find some
			var biggest_size = 0;
			var biggest_feature;
			var furthest_point = 0;
			for (fx in features) {
				var f = features[fx];
				if (f.radius == rad && f.type != ft_enzyme) { 
					var new_size = f.end() - f.start() + 1;
					var overlap = furthest_point - f.start();
					if (overlap <= 0) { 
						// We've cleared all potential conflicts: reset
						// the indicators
						biggest_size = new_size;
						biggest_feature = f;
						furthest_point = f.end();
					} else if (biggest_size > min_overlap_feature_size &&
					           new_size > min_overlap_feature_size &&
							   overlap/biggest_size > min_overlap_pct &&
							   overlap/new_size > min_overlap_pct) {
						// Overlap: conflict!
						if (new_size > biggest_size) { // This feature is top dog,
						                               // move the original to the
						                               // new radius
							biggest_feature.radius = new_rad;
							biggest_size = new_size;
							biggest_feature = f;
							furthest_point = f.end();
							conflicts++;

						} else { // The original feature is top dog. move the new
						         // feature to the new radius
							f.radius = new_rad;
							conflicts++;
						}

					}
				}
			}

			// Keep alternating between inside and outside the plasmid.
			rad = new_rad;
			rx++;
		} while (conflicts > 0); // Keep adding levels of resolution

	}

	///////////////////////////////////////////////////////////////////
	// MAIN ENTRY POINT
	var paper = Raphael("plasmid-map", map_width, map_height);
	draw_plasmid();

	features = parse_features();
	resolve_conflicts(features);
	
	for (fx in features) {
		features[fx].draw();
	}

})
