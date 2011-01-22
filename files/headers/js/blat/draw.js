jQuery(document).ready(function ($) {
	var _debug = false;
	///////////////////////////////////////////////////////////////////
	// Global parameters
	// Paper setup
	var map_width = 640;
	var map_height = 640;
	var cx = map_width/2;
	var cy = map_height/2;

	// Global plasmid info
	var seq_length = parseInt($("#sequence-data-length").html());

	// Loop radii
	var radius_spacing = 20; // spacing
	var plasmid_radius = 200;
	var inner_radius = plasmid_radius - radius_spacing; 
	var outer_radius = plasmid_radius + radius_spacing;
	var label_radius_offset = 0;

	// Feature visual properties
	var feature_width = 15;
	var enzyme_width = 25;
	var enzyme_weight = 1; // Restriction enzymes are drawn differently
	                       // This controls their "thickness" on the map
	var feature_opacity = 1.0/3.0;
	var enzyme_opacity = 1.0;
	var head_width = 25;
	var head_length = 7;

	// Animation properties
	var fade_time = 300;

	// Overlaps
	var min_overlap_cutoff = -0.1;// in degrees
	var min_overlap_pct = 0.01;
	var min_overlap_feature_size = 0.5; // in degrees
	
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
	var color_enzyme  = "#00c";

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

	///////////////////////////////////////////////////////////////////
	// Internals start here
	
	// Global label height list: keeps track of the current heights of each of
	// the 8 label lists, so that we know at what height to add the next label.
	var label_heights = [0, 0, 0, 0, 0, 0, 0, 0];
	
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
		seq_length_to_angle: function (l) {
			//     just like pos_to_angle, but without caring about the start
			return (l/seq_length) * 360;
		},
		angle_to_pos: function (a) {
			//     start at the top of the circle
			return Math.round(1 + ((seq_length - 1)/360) * ((360 + 90 - a) % 360));
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


	///////////////////////////////////////////////////////////////////
	// Feature class
	function Feature(feature_list) {

		// Private data members, from the argument list
		var _name = feature_list[0];
		var _start = parseInt(feature_list[1]);
		var _end = parseInt(feature_list[2]);
		var _type = feature_list[3];
		var _clockwise = (_start <= _end);

		// Because we store the clockwise information in a separate variable,
		// ensure that start and end go clockwise (end > start), regardless
		// of what the original feature data was
		if (!_clockwise) {
			var _tmp = _start;
			_start = _end;
			_end = _tmp;
		}

		// Type-based property selection
		var _color = color_feature;
		var _width = feature_width;
		var _draw_head = false;
		var _opacity = feature_opacity;
		var _opaque = false; // holds opacity for clicks
		switch(_type) {
			case ft_promoter:
			case ft_primer:
				_draw_head = true; // Promotors and primers are the only primer-colored
			case ft_terminator:   // features with heads
				_color = color_primer;
				break;
			case ft_regulatory:
			case ft_origin:
				_color = color_origin;
				break;
			case ft_enzyme:
				_color = color_enzyme;
				_width = enzyme_width;
				_opacity = enzyme_opacity;
				break;
			case ft_gene:
				_draw_head = true;
				break;
		}

		var _this = this; // another copy of the this pointer, to 
		              // work around a JavaScript bug that reassigns (this)
		              // improperly

		// Radius is public, unlike other properties, which are permanent
		this.radius = plasmid_radius; // Default to plasmid radius, can be changed
		                              // later on by other methods

		// Accessors for private properties set at creation
		this.name = function() { return _name; };
		this.start = function() { return _start; };
		this.end = function() { return _end; };
		this.type = function() { return _type; };
		this.clockwise = function() { return _clockwise; };

		// Calculated properties
		
		// Degree conversion, for overlap calculation:
		// for these functions, the sequence starts at 90 degrees and goes down.
		this.start_degrees = function() {
			var sd;
			// Take the minimum head size into account. Only need to do this 
			// when the head is drawn and pointing clockwise, to
			// "push the start back."
			if (_draw_head && _clockwise) { 
				sd = convert.pos_to_angle(_end) + _this.size_degrees();
			} else { // Headless feature, or head is pointing the wrong way.
				     // Just give its typical start position
				sd = convert.pos_to_angle(_start);
			}
			return sd;
		}
		this.end_degrees = function() {
			var ed;
			// Take the minimum head size into account. Only need to do this 
			// when the head is drawn and pointing counterclockwise, to 
			// "push the end forward."
			if (_draw_head && !_clockwise) { // Take the minimum head size into account
				ed = convert.pos_to_angle(_start) - _this.size_degrees();
			} else { // Headless feature, or head is pointing the wrong way.
				     // Just give its typical end position
				ed = convert.pos_to_angle(_end);
			}
			return ed;
		}
		this.size_degrees = function() {
			var szd; // size in degrees
			// Normal definition of size
			szd = convert.seq_length_to_angle(_this.end() - _this.start() + 1);

			// Head size: return this if it's bigger
			if (_draw_head) {
				// Convert the head length into degrees, just as you do
				// in the draw() method. Must recalcualte every time, as
				// radius may have changed
				var r_p = Math.sqrt(_this.radius*_this.radius + 
						head_length*head_length);
				var hszd = Raphael.deg(Math.asin(head_length/r_p));
				if (hszd > szd)
					szd = hszd;
			}

			return szd;
		}


		// The visual object to modify when accessing the feature.
		var draw_feature = paper.set();

		// Feature drawing
		this.draw = function () {


			// Convert from sequence positions to angles
			var a0 = convert.pos_to_angle(_start);
			var a1 = convert.pos_to_angle(_end);

			// Create the draw feature, a set which will have the head 
			// and arc pushed onto it as necessary.
			
			// Arrowhead drawing, if needed
			if (_draw_head) {

				// Arrow tip point lines up with a0 or a1 and it points
				// tangent to the circle.
				// We need to figure out how many radians the arrow takes up
				// in order to adjust a0 or a1 by that amount, and to set the 
				// base of the triangle even with that angle
				var r_p = Math.sqrt(_this.radius*_this.radius + 
						head_length*head_length);
				// "height" of the arrowhead, in degrees
				var a_b;
				var a_p = Raphael.deg(Math.asin(head_length/r_p));
				// Adjust the appropriate edge to compensate for the arrowhead
				if (_clockwise) {
					a_b = (a1 + a_p) % 360 ; // base angle
					a_p = a1;       // point angle
					a1  = a_b;      // adjust arc edge
				} else {
					a_b = (a0 - a_p) % 360 ; // base angle
					a_p = a0;       // point angle
					a0  = a_b;      // adjust arc edge
				}
				var xy_p = convert.polar_to_rect(_this.radius, a_p);
				
				// bottom and top points, rectangular
				var xy_b = convert.polar_to_rect(_this.radius - head_width/2.0, a_b);
				var xy_t = convert.polar_to_rect(_this.radius + head_width/2.0, a_b);

				// Unlike the arc, the head is traced with a line, and
				// then created entirely with the fill color
				var head = paper.path(svg.move(xy_p.x, xy_p.y) +
									  svg.line(xy_b.x, xy_b.y) +
									  svg.line(xy_t.x, xy_t.y) + 
									  svg.close());
				head.attr({"stroke-width": 0,
						   "fill":         _color});
				draw_feature.push(head);
			}

			// Arc drawing
			if (a1 < a0 && _type != ft_enzyme) { 
				// Compensating for the head may have "taken up" all
				// the room on the plasmid, in which case no arc needs
				// to be drawn

				// Rectangular coordinates of the edges of the arc: 
				// arcs are drawn counterclockwise, even though the plasmid
				// sequence increases clockwise, so we flip the
				// indices
				var xy0 = convert.polar_to_rect(_this.radius, a1);
				var xy1 = convert.polar_to_rect(_this.radius, a0);

				// The arc has no fill-color: it's just a thick line
				var arc = paper.path(svg.move(xy0.x, xy0.y) +
									 svg.arc(_this.radius, xy1.x, xy1.y));
				arc.attr({"stroke-width": _width});

				draw_feature.push(arc);
			} else if (_type == ft_enzyme) { 
				// Restriction enzymes get drawn on their own
				var xy0 = convert.polar_to_rect(_this.radius - enzyme_width/2.0, 
						(a0+a1)/2.0);
				var xy1 = convert.polar_to_rect(_this.radius + enzyme_width/2.0, 
						(a0+a1)/2.0);
				// Not really an arc, just a line, but left this way
				// for consistency
				var arc = paper.path(svg.move(xy0.x, xy0.y) +
									 svg.line(xy1.x, xy1.y));
				arc.attr({"stroke-width": enzyme_weight});
				arc.toBack();

				draw_feature.push(arc);
				_opaque = true;
			}

			// Apply the feature-wise properties to the whole feature
			draw_feature.attr({"stroke":         _color,
			                   "stroke-linecap": "butt",
			                   "opacity":        _opacity,
			                   "title":          _name});


		} // END Feature::draw()

		// Draw the label associated with that feature
		this.draw_label = function (r_l) {
			// Figure out the center of the feature
			a_c = (_this.start_degrees() + _this.end_degrees()) / 2.0;
			xy0 = convert.polar_to_rect(_this.radius, a_c);
			
			// Figure out the label position: divide the grid up into eight
			// sections
			section = Math.floor((90 - a_c) / 45);
			section_angle = 90 - 45/2.0 - section*45;
			feature_section_extent = (a_c - (section_angle + 45/2.0))/45;

			xy1 = convert.polar_to_rect(r_l, section_angle);

			y_shift = label_heights[section];
			// TODO: FIX THIS TO ACCOUNT FOR OTHER LABELS!
			if (xy1.y > cy) { // Lower half: add below
				xy1.y += y_shift;
			} else { // Upper half: add above
				xy1.y -= y_shift;
			}

			// Draw the line to the label position
			var label_line = paper.path(svg.move(xy0.x, xy0.y) +
										svg.line(xy1.x, xy1.y));
			label_line.attr({"stroke": color_bg_text,
			                 "opacity": feature_opacity});


			var label = paper.text(xy1.x, xy1.y, _name);
			if (a_c < -90 && a_c > -270) { // Left half of wheel: align right
				label.attr({"text-anchor": "end"});
			} else if (a_c < 90 && a_c > -90) { // Right half of wheel: align left
				label.attr({"text-anchor": "start"});
			} // Top and bottom default to middle, which is correct

			// Update the label heights
			label_heights[section] += label.getBBox().height;

			label.attr({"fill": _color,
			            "opacity": feature_opacity});

			draw_feature.push(label);
			draw_feature.push(label_line);

		} // END Feature::draw_label(r_l)

		// Register click and hover handlers to make the feature "active"
		this.register_handlers = function () {
			// Toggle solid/light upon click
			draw_feature.click(function (event) {
				if (_opaque) {
					draw_feature.animate({"opacity": feature_opacity}, fade_time);
					_opaque = false;
				} else {
					draw_feature.animate({"opacity": 1}, fade_time);
					_opaque = true;
				}
			});
			draw_feature.hover(
				function (event) {
					if (!_opaque)
						draw_feature.animate({"opacity": 1}, fade_time);
				}, 
				function (event) {
					if (!_opaque)
						draw_feature.animate({"opacity": feature_opacity}, fade_time);
				} 
			);
		} // END Feature::register_handlers()
	} // END Feature Class

	
	///////////////////////////////////////////////////////////////////
	// Global functions

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

		for (var ang = 0; ang < 360; ang += 30) {
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
			var feat = new Feature(row.slice(2)); // row[0:2] are the checkboxes
			features.push(feat);
		})

		return features;
	}
	
	// Move features that overlap to other radii.
	function resolve_conflicts(features) {
		var conflicts;
		var rad = plasmid_radius; // current radius
		var rx = 1;               // radius counter
		var max_rad = plasmid_radius;

		function push(winner, loser) {
			// Record that the push happened
			winner.pushed_features.push(loser); 
			conflicts++;

			// Do it
			loser.radius = new_rad; 

			if (_debug) console.warn(loser.name() + " pushed by " + winner.name());
			
			// Since loser was pushed, un-push all the 
			// features it caused to be pushed, as long as
			// those features were not in conflict with others.
			for (var pfx in loser.pushed_features) {
				var pf = loser.pushed_features[pfx];
				// Check for conflict with other the winner feature itself.
				// If there's no conflict, we can pushh it back safely.
				if (pf.start_degrees() - winner.end_degrees() <= min_overlap_cutoff ||
					winner.start_degrees() - pf.end_degrees() <= min_overlap_cutoff) {
					if (_debug)
						console.warn(pf.name() + "unpushed, because " 
							+ loser.name() + " pushed by " + winner.name());
					pf.radius = rad;
				}
			}
		}

		do {
			// Keep alternating between inside and outside the plasmid.
			var new_rad = rad + Math.pow(-1, rx) * rx * radius_spacing;

			conflicts = 0; // Assume you have no conflicts until you find some

			// Clear the record of who pushed whom
			for (var fx in features) {
				features[fx].pushed_features = [];
			}

			var biggest_size = 0;
			var biggest_feature;
			var furthest_point = 90; // Start at the top of the circle
			for (var fx in features) {
				var f = features[fx];
				if (f.radius == rad && f.type() != ft_enzyme) { 
					var new_size = f.size_degrees();
					var overlap = -(furthest_point - f.start_degrees());
					if (overlap <= min_overlap_cutoff) { 
						// We've cleared all potential conflicts: reset
						// the indicators
						biggest_size = new_size;
						biggest_feature = f;
						furthest_point = f.end_degrees();
					} else if (biggest_size > min_overlap_feature_size &&
						       new_size > min_overlap_feature_size &&
						      (overlap <= 0 || 
						      (overlap/biggest_size > min_overlap_pct &&
						       overlap/new_size > min_overlap_pct))) {
						// Overlap: conflict!
						if (new_size > biggest_size) { // This feature is top dog,
						                               // move the original to the
						                               // new radius
							push(f, biggest_feature);

							// Update the new top dog
							biggest_size = new_size;
							biggest_feature = f;
							furthest_point = f.end_degrees();

						} else { // The original feature is top dog. move the new
						         // feature to the new radius

							push(biggest_feature, f);
						}

					}
				}
			}

			// Keep track of the biggest radius reached
			if (rad > max_rad)
				max_rad = rad;

			// Move on to the next radius
			rad = new_rad;
			rx++;

			
		} while (conflicts > 0); // Keep adding levels of resolution

		return max_rad;
	}

	function draw_features(features) {
		for (var fx in features) {
			features[fx].draw();
		}
	}

	function label_features(features, max_radius) {
		var label_radius = max_radius + label_radius_offset; 

		// Iterate counterclockwise
		for (var fx = features.length - 1; fx >= 0; fx--) {
			if (features[fx].type() != ft_enzyme) 
				features[fx].draw_label(label_radius);
		}
	}

	function register_handlers(features) {
		for (var fx in features) {
			features[fx].register_handlers();
		}
	}


	///////////////////////////////////////////////////////////////////
	// MAIN ENTRY POINT
	var paper = Raphael("plasmid-map", map_width, map_height);

	draw_plasmid();
	var features = parse_features();
	var max_radius = resolve_conflicts(features);
	draw_features(features);
	label_features(features, max_radius);
	register_handlers(features);
	
})
