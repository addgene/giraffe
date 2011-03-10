/////////////////////////////////////////////////////////////////////////////
//
// GiraffeDraw: Javascript API for drawing plasmid maps from 
//              JSON plasmid feature data
//
// 2011 Addgene
// Mikhail Wolfson (wolfsonm@addgene.org)
// Benjie Chen     (benjie@addgene.org)
//

/////////////////////////////////////////////////////////////////////////////
//
// Drawing API: caller should call the GiraffeDraw function to get a drawing
// and parsing object. The caller should use the object's read() method to
// read in a JSON list of features, then call one of the object's draw_<foo>_map
// methods, with a customized set of options.
//
// The read() method can be passed as the JSONP argument to the
// BLAT get API (i.e. API to retrieve array of features of a sequence
// from the server):
//
//    <script src="/headers/js/raphael-min.js"></script>
//    <script src="/headers/js/scale.raphael.js"></script>
//    <script src="http://host/api/js/draw.js"></script>
//    <script> var gd = GiraffeDraw(); </script>
//    <script src="http://host/blat/8de36469..../default?jsonp=gd.read">
//    </script>
//
// After the feature data is read in, one (or many) plasmid maps can be drawn.
// Options to the drawing function are passed in via a dictionary argument. For
// example:
//
//    <script> 
//    gd.draw_circular_map({ "map_dom_id" : "some_id", ... })
//    </script>
//
/// Available options are:
//
//  map_dom_id: ID of the DOM element that contains the map. Default
//  is "giraffe-draw-map"
//
//  fade_time: if non-zero, then highlight a feature will cause an
//  animated fade-in/out effect. Default is 0.
//
//  feature_opacity: opacity when feature is shown. if not 1.0, then
//  when feature is moused over or clicked on, the opacity will become
//  1.0. Default is 0.7.
//
//  map_width, map_height: default 640, 640.
//
//  plasmid_name: if given, show this name together with size of
//  sequence in the middle of the plasmid. Default is "".
//
//  label_offset: how far from the outter feature should we
//  start drawing labels. Default is "10".
//
//  cutters: which kinds of restriction enzymes to show, if any. This
//  list of integers is interpreted as follows: [1, 2]: show 1- and 2-
//  cut restriction enzymes. []: show nothing. etc. Default is [1].
//
//
//
// FRAMEWORK
// GiraffeDraw()
// |
// -- basic initialization of list of Features
// |
// -- read(JSON)
// |  feature parsing into "local" features object
// |
// -- draw_circular_map()
// |  circular feature drawing, which creates clones of the original Feature 
// |  objects with extended properties
// |
// -- draw_linear_map() (maybe) 
//    linear feature drawing, which makes its own extended clones of objects

// Protect scope, but ensure that GiraffeDraw() is global
(function(){window.GiraffeDraw = function () {

	///////////////////////////////////////////////////////////////////
	// Package-scope variables
	var _debug = false;
	var basic_features = [];
	var seq_length;

	// Feature Types
	var ft = { 
		feature:    1, promoter:   2, primer:        3,
		enzyme:     4, gene:       5, origin:        6,
		regulatory: 7, terminator: 8, exact_feature: 9,
        orf:       10
	};

	//// Package-scope utility functions
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

	// Feature Colors
	var colors = {
		bg_text: "#aaa",
		plasmid: "#000",
		feature: "#f00",
		primer:  "#090",
		origin:  "#333",
		enzyme:  "#00c",
		orf:     "#00c8c8",
	};


	///////////////////////////////////////////////////////////////////
	// Basic Feature class (private, internal class)
	function Feature(feat) {

		// Private data members, from the properties of the feature argument
		// since hash tables and objects are one and the same in js
		var _name = feat.feature;
		var _start = parseInt(feat.start);
		var _end = parseInt(feat.end);
		var _type = parseInt(feat.type_id);
        var _default_show_feature = 1;
        // Not all features have this option, e.g. annotated ones we
        // have to show the feature.
        if ('show_feature' in feat) { _default_show_feature = parseInt(feat.show_feature); }
		var _clockwise = feat.clockwise;
		var _cut = parseInt(feat.cut); // only for enzymes;
		var _other_cutters = []; // only for enzymes;

		// Accessors for private properties set at creation
		this.name = function() { return _name; };
		this.start = function() { return _start; };
		this.end = function() { return _end; };
		this.type = function() { return _type; };
		this.clockwise = function() { return _clockwise; };
        this.default_show_feature = function() { return _default_show_feature; }
		// returns - 1 if not enzyme		
		this.cut = function() { return _type == ft.enzyme ? _cut : -1 }; 

		this.crosses_boundary = function () { return _end < _start };

		// Enzyme-only data access methods
		// gives 0 if not enzyme
		this.cut_count = function() { return _other_cutters.length; };
		// gives [] if not enzyme
		this.other_cutters = function() { return _other_cutters; };
		// Mutator for other_cutters: only for enzymes
		this.set_other_cutters = function(c) {
			if (_type == ft.enzyme)
				_other_cutters = c;
		}

	}; // END Basic Feature Class

	///////////////////////////////////////////////////////////////////
	// JSON Parsing
	this.read = function (features_json) {
		basic_features = []; // package scope
		seq_length = features_json[0]; // package scope

		for (var ix = 1; ix < features_json.length; ix++) {
			basic_features.push(new Feature(features_json[ix]));
		}

		// Now that the features are parsed, calculate how many instances
		// of each cutter type there are.
		cut_counts(); 
	}

	// Private function
	// Calculate cut counts of all restriction enzyme
	function cut_counts() {
		var cut_counts = {}; 
		var f;
		// Calculate the counts
		for (var fx in basic_features) {
			f = basic_features[fx];
			if (f.type() == ft.enzyme) {
				// Store indices, not Feature objects, because
				// the objects will change, depending on
				// the kind of map is drawn
				if (f.name() in cut_counts)
					cut_counts[f.name()].push(fx);
				else
					cut_counts[f.name()] = [fx];
			}
		}
		// Store them for each enzyme feature
		for (var fx in basic_features) {
			f = basic_features[fx];
			if (f.type() == ft.enzyme)
				f.set_other_cutters(cut_counts[f.name()]);
		}
	}


	///////////////////////////////////////////////////////////////////
	// Circular Map Drawing Class
	this.CircularMap = function(options) {

		// Map-specific canvas element
		var paper;

		// Map-specific feature list
		var features = [];

		// Paper setup - not the final width, but how we will draw the
		// map, we will scale later on
		var map_width = 800;
		var map_height = 800;
		var cx = map_width/2;
		var cy = map_height/2;

		// Where to draw the map
		var map_dom_id = 'giraffe-draw-map';
		if ('map_dom_id' in options) {
			map_dom_id = options['map_dom_id'];
		}

		// Final size
		var final_map_width = 640;
		var final_map_height = 640;
		if ('map_width' in options) {
			final_map_width = parseInt(options['map_width'])
		}
		if ('map_height' in options) {
			final_map_height = parseInt(options['map_height'])
		}

		// Global plasmid info
		var plasmid_start = 90; // degrees

		// Loop radii
		var radius_spacing = 20; // spacing
		var plasmid_radius = 200;
		var inner_radius = plasmid_radius - radius_spacing; 
		var outer_radius = plasmid_radius + radius_spacing;
		var label_radius_offset = 10;
		if ('label_offset' in options) {
			label_radius_offset = parseInt(options['label_offset']);
		}

		// Feature visual properties
		var feature_width = 15;
		var enzyme_width = 25;
		var enzyme_weight = 1; // Restriction enzymes are drawn differently
							   // This controls their "thickness" on the map
		var enzyme_bold_weight = 3; // How thick when highlighted
		var feature_opacity = 0.7;
		var enzyme_opacity = 0.7;
		if ('opacity' in options) {
			feature_opacity = parseFloat(options['opacity']);
			enzyme_opacity = parseFloat(options['opacity']);
		}
		var bold_opacity = 1.0;
		var head_width = 25;
		var head_length = 7;

		// Cutters to show
		var cutters_to_show = [1];
		if ('cutters' in options) {
			cutters_to_show = options['cutters'];
		}

		// Animation properties
		var fade_time = 0;
		if ('fade_time' in options) {
			fade_time = parseInt(options['fade_time'])
		}

		// Overlaps
		var min_overlap_cutoff = -0.1;// in degrees
		var min_overlap_pct = 0.01;
		var min_overlap_feature_size = 0.5; // in degrees
		
		// Tic marks
		var tic_mark_length = 15;
		var tic_mark_radius = inner_radius - tic_mark_length/2;
		var tic_label_radius = tic_mark_radius - 1.5*tic_mark_length;

		// Labels and other text
		var label_line_weight = 1;
		var label_line_bold_weight = 1.5 * label_line_weight;
		var label_font_size = '13pt';
		var plasmid_font_size = '16pt';
		var plasmid_name = '';
		if ('plasmid_name' in options) {
			plasmid_name = options['plasmid_name'];
		}

        var label_letter_height = 0;
        var label_letter_width = 0;

		// This is based on using 8 label lists, which is pretty much hard
		// coded in a few places, so don't change this unless you figure
		// out how to change the number of label lists.
		var label_section_degree = 45;

		///////////////////////////////////////////////////////////////////
		// Internals start here

		// Conversion Object: local to this function because it assumes a 
		// circular geometry
		// Groups conversions into handy section
		var convert = {
			pos_to_angle: function (p) {
				//     start at the top of the circle
				return plasmid_start - (p/seq_length) * 360;
			},
			seq_length_to_angle: function (l) {
				//     just like pos_to_angle, but without caring about the start
				return (l/seq_length) * 360;
			},
			angle_to_pos: function (a) {
				//     start at the top of the circle
				return Math.round(1 + ((seq_length - 1)/360) * ((360 + plasmid_start - a) % 360));
			},
			/* Angle a is in degrees from the horizontal, counterclockwise, and
			 * r is relative to the center of the paper */
            polar_to_rect_center_at_zero: function (r, a) {
				var rect = {};
				rect.x = r * Math.cos(Raphael.rad(a));
				rect.y = r * Math.sin(Raphael.rad(a));
				return rect;
            },
			polar_to_rect: function (r, a) {
			    var rect = convert.polar_to_rect_center_at_zero(r,a);
				rect.x += cx;
				// Coordinates increase as you go down, so y is flipped.
                rect.y = cy - rect.y;
				return rect;
			}
		};

        function label_list_section_angle(section) {
		    return plasmid_start-label_section_degree/2.0-section*label_section_degree;
        };

		///////////////////////////////////////////////////////////////////
		// Circular Feature class
		function CircularFeature(basic_feature) {
			// Clone the basic feature, to extend it later
			function Clone() {}; // empty function to use as a hanger for prototype
			Clone.prototype = basic_feature;
			var _this = new Clone(); // Make a new "this" pointer to return at
			                         // the end
			// The result of this function will be a CircularFeature object
			
			// Visual properties
			var _visible = true;
			var _labeled = true;

			// Type-based property selection
			var _color = colors.feature;
			var _width = feature_width;
			var _draw_head = false;
			var _opacity = feature_opacity;
			var _opaque = false; // holds opacity for clicks
			switch(_this.type()) {
				case ft.promoter:
				case ft.primer:
					_draw_head = true; // Promotors and primers are the only primer-colored
				case ft.terminator:   // features with heads
					_color = colors.primer;
					break;
				case ft.regulatory:
				case ft.origin:
					_color = colors.origin;
					break;
				case ft.enzyme:
					_color = colors.enzyme;
					_width = enzyme_width;
					_opacity = enzyme_opacity;
					break;
                case ft.orf:
                    _color = colors.orf;
				case ft.gene:
					_draw_head = true;
					break;
			}

			// Radius is public, unlike other properties, which are permanent
			_this.radius = plasmid_radius; // Default to plasmid radius, can be changed
										  // later on by other methods

			// Accessors for private properties set at creation
			_this.visible = function() { return _visible; };
			_this.labeled = function() { return _labeled; };

			// The visual object to modify when accessing the feature.
			var _feature_set;
			var _arrow_set;
			var _label_set;
			var _label_drawn = false;

			// Check to see if label has been drawn yet
			_this.label_drawn = function() { return _label_drawn; };
			_this.feature_set = function() { return _feature_set };
			_this.label_set = function() { return _label_set };

            _this.initialize = function() {
			    _feature_set = paper.set();
			    _arrow_set = paper.set();
			    _label_set = paper.set();
            }

			// Calculated properties

			// Degree conversion, for overlap calculation:
			// for these functions, the sequence starts at 90 degrees and goes down.
			_this.start_degrees = function() {
				var sd;
				// Take the minimum head size into account. Only need to do this 
				// when the head is drawn and pointing clockwise, to
				// "push the start back."
				if (_draw_head && _this.clockwise()) { 
					sd = convert.pos_to_angle(_this.end()) + _this.size_degrees();
				} else { // Headless feature, or head is pointing the wrong way.
						 // Just give its typical start position
					sd = convert.pos_to_angle(_this.start());
				}
				return sd;
			};

			_this.end_degrees = function() {
				var ed;
				// Take the minimum head size into account. Only need to do this 
				// when the head is drawn and pointing counterclockwise, to 
				// "push the end forward."
				if (_draw_head && !_this.clockwise()) { // Take the minimum head size into account
					ed = convert.pos_to_angle(_this.start()) - _this.size_degrees();
				} else { // Headless feature, or head is pointing the wrong way.
						 // Just give its typical end position
					ed = convert.pos_to_angle(_this.end());
				}

				// For boundaries that cross the 0 point, ensure that they return
				// a number in the same range as everyone else
				if (_this.crosses_boundary()) 
					ed += 360;

				return ed;
			};

			_this.size_degrees = function() {
				var szd; // size in degrees
				// Normal definition of size
				if (_this.crosses_boundary()) 
					// Start and end are flipped here: non-intuitive
					szd = convert.seq_length_to_angle(seq_length - _this.start() + _this.end() + 1);
				else
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
			};

			// Actions for interactivity
			// Generic fading animation/property setting mechanism
			var _fade = function (props, line_props) {
				var sets = paper.set();
				var lines = paper.set();
				sets.push(_feature_set);
				lines.push(_label_set[0]); // label line

                // Cutters: highlight other examples of this enzyme if
                // it's a multi-cutter
                if (_this.type() == ft.enzyme) {
					for (var fx in _this.other_cutters()) {
						var f = features[_this.other_cutters()[fx]];
						sets.push(f.feature_set());
						lines.push(f.label_set()[0]);
					}
				}

				if (fade_time) { 
					sets.animate(props, fade_time); 
					lines.animateWith(sets, line_props, fade_time); 
				} else { 
					sets.attr(props);
					lines.attr(line_props); 
				}
			}

			var _bolder = function () {
				var props = {"opacity": bold_opacity,
                             "font-weight": "bold" };
                if (_this.type() == ft.enzyme) {
				    props["stroke-width"] = enzyme_bold_weight;
                }
				var line_props = {"stroke": colors.plasmid, 
				                  "stroke-width": label_line_bold_weight};
				_fade(props, line_props);
			}

			var _lighter = function () {
				var props = {"opacity": _opacity,
                             "font-weight":"normal"};
                if (_this.type() == ft.enzyme) {
				    props["stroke-width"] = enzyme_weight;
                }
				var line_props = {"stroke": colors.bg_text,
				                  "stroke-width": label_line_weight};
				_fade(props, line_props);
			}

			// Toggle solid/light upon click
			var _click = function (event) {
				if (_opaque) {
					_lighter();
					_opaque = false;
				} else {
					_bolder();
					_opaque = true;
				}
			};

			// Hovering: solid/light upon mouseover
			var _mouse_over = function (event) {
				if (!_opaque)
					_bolder();
			};
			var _mouse_up = function (event) {
				if (!_opaque)
					_lighter();
			};

			// Feature drawing
			_this.draw = function () {
                if (!_visible) { return; }

				// Convert from sequence positions to angles
				var a0 = convert.pos_to_angle(_this.start());
				var a1 = convert.pos_to_angle(_this.end());

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
					if (_this.clockwise()) {
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
					_arrow_set.push(head);
				}

				// Arc drawing
				if ((_this.crosses_boundary() || a1 < a0) && _this.type() != ft.enzyme) { 
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

					_arrow_set.push(arc);
				} else if (_this.type() == ft.enzyme) { 
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

					_arrow_set.push(arc);
				}

				_arrow_set.click(_click);
				_arrow_set.hover(_mouse_over, _mouse_up);

				_feature_set.push(_arrow_set);

				// Apply the feature-wide properties to the whole feature
				_feature_set.attr({"stroke":         _color,
								   "stroke-linecap": "butt",
								   "opacity":        _opacity,
								   "title":          _this.name()});

			} // END CircularFeature::draw()

			// Should we draw the label?
			_this.should_draw_label = function () {
				// Don't bother unless we need to
				if (!_visible || !_labeled) 
					return false;
				return true;
			} // END CircularFeature::should_draw_label()

            // What label should we draw?
            _this.label_name = function () {
				var label_name = _this.name();
				if (_this.type() == ft.enzyme) {
					label_name += " (" + _this.cut() + ")";
				}
                return label_name;
            } // END CircularFeature::label_name()

			// Set which label list this feature's label should be in
			_this.set_label_list = function () {
				if (!_this.should_draw_label()) { return; }

				// Figure out the center of the feature
				var a_c = (_this.start_degrees()+_this.end_degrees())/2.0;
				var adjust_a_c = a_c;
				if (adjust_a_c < 0) { adjust_a_c += 360; }

                // Figure out which section this label is in: divide
                // the grid up into eight sections.
				var section = Math.floor((plasmid_start-a_c)/label_section_degree);

				var l = label_f_c[section].length;
				label_f_c[section][l] = [adjust_a_c,_this.label_name()];

			} // END CircularFeature::set_label_list()

			// Draw the label associated with that feature
			_this.draw_label = function () {
				if (!_this.should_draw_label()) { return; }

				if (_label_drawn)
					_this.clear_label();

				// Figure out the center of the feature
				var a_c = (_this.start_degrees()+_this.end_degrees())/2.0;
				var adjust_a_c = a_c;
				if (adjust_a_c < 0) { adjust_a_c += 360; }

				var xy0 = convert.polar_to_rect(_this.radius, a_c);
				
                // Figure out which section this label is in: divide
                // the grid up into eight sections.
				var section = Math.floor((plasmid_start - a_c)/label_section_degree);
				// Figure out position in the label list.
				var pos_ls = 0
				for (pos_ls=0; pos_ls<label_f_c[section].length; pos_ls++) {
					if (label_f_c[section][pos_ls][0] == adjust_a_c) {
						// Claim this spot, 999 is not a valid value for
						// adjust_a_c.
						label_f_c[section][pos_ls][0] = 999;
						break;
					}
				}
				var sec_labels = label_f_c[section].length;
				var y_shift = pos_ls*label_letter_height;
				var xy1 = {}
				xy1.x = label_list_pos[section][0]
				xy1.y = label_list_pos[section][1]

				// We want to minimize the number of label lines that
				// cross. Which means depends on which section we are in,
				// we draw labels in different orders. See draw_labels on
                // how the positions of each label section are
                // computed. Remember, because we are sorted by
                // label_f_c, we are going counter clockwise on the
                // circle, drawing label for the feature with higher
                // bp first. label_f_c has the lower x and y
                // coordinates of each section.
				if (section == 0 || section == 1) {
					// upper right, higher bp at bottom
					xy1.y -= y_shift;
				}
				else if (section == 2 || section == 3) {
					// lower right, higher bp at bottom
					xy1.y -= y_shift;
				}
				else if (section == 4 || section == 5) {
					// lower left, higher bp on top
					xy1.y = xy1.y - sec_labels*label_letter_height + y_shift;
				}
				else if (section == 6 || section == 7) {
					// upper left, high bp on top
					xy1.y = xy1.y - sec_labels*label_letter_height + y_shift;
				}

				// Draw the line to the label position
				var label_line = paper.path(svg.move(xy0.x, xy0.y) +
											svg.line(xy1.x, xy1.y));
				label_line.attr({"stroke": colors.bg_text,
				                 "stroke-width": label_line_weight,
								 "opacity": 0.5 });

				// Enzymes show their cut sites in the label
				var label_name = _this.label_name();
				var label = paper.text(xy1.x, xy1.y, label_name);

                if (section > 3) {
					// Left half of wheel: align right
					label.attr({"text-anchor": "end"});
				}
                else {
					// Right half of wheel: align left
					label.attr({"text-anchor": "start"});
				}

				label.attr({"fill": _color,
							"font-size": label_font_size,
							"opacity": 1.0 });

				_label_set.push(label_line);
				_label_set.push(label);

				// Handlers
				_label_set.click(_click);
				_label_set.hover(_mouse_over, _mouse_up);

                // Only push label_line, so when we fade in and out,
                // we don't also fade the label.
				_feature_set.push(label_line);

				_labeled = true;
				_label_drawn = true;
			} // END CircularFeature::draw_label()

			_this.hide = function () {
				if (_visible) {
                    if (_feature_set) { _feature_set.hide(); }
					_visible = false;
					_labeled = false;
				}
			}; // END CircularFeature::hide()

			_this.show = function () {
				if (!_visible) {
					if (_feature_set) { _feature_set.show(); }
					if (!_labeled) {
                        if (_label_set) { _label_set.hide(); }
                    }
					_visible = true;
				}
			}; // END CircularFeature::show()

			_this.hide_label = function () {
				if (_labeled) {
                    if (_label_set) { _label_set.hide(); }
					_labeled = false;
				}
			}; // END CircularFeature::hide_label()

			_this.show_label = function () {
				if (!_labeled) {
                    if (_label_set) { _label_set.show(); }
					_labeled = true;
				}
			}; // END CircularFeature::show_label()

			_this.clear_label = function () {
				if (_label_drawn) {
                    if (_label_set) {
					    _label_set.unclick(_click);
					    _label_set.unhover(_mouse_over, _mouse_up);
					    _label_set.remove();
					    _label_set = paper.set();
                    }
					_labeled = false;
					_label_drawn = false;
				}
			}; // END CircularFeature::clear_label()

			return _this;
		}; // END CircularFeature Class

		///////////////////////////////////////////////////////////////////
		// Drawing utility functions

		// Circle setup
		function draw_plasmid() {
			function draw_tic_mark(a) {
				var r0 = tic_mark_radius - tic_mark_length/2;
				var r1 = tic_mark_radius + tic_mark_length/2;
				var xy0 = convert.polar_to_rect(r0,a);
				var xy1 = convert.polar_to_rect(r1,a);
				var tic = paper.path(svg.move(xy0.x, xy0.y) +
									 svg.line(xy1.x, xy1.y));
				tic.attr({"stroke": colors.bg_text});

				var xyl = convert.polar_to_rect(tic_label_radius, a);
				var label = paper.text(xyl.x, xyl.y, String(convert.angle_to_pos(a)));
				label.attr({"fill": colors.bg_text});
				if (a < plasmid_start || a > 360 - plasmid_start) { // Right half of wheel: align right
					label.attr({"text-anchor": "end"});
				} else if (a > plasmid_start && a < 360 - plasmid_start) { // Left half of wheel: align left
					label.attr({"text-anchor": "start"});
				} // Top and bottom default to middle, which is correct
			}

			var plasmid = paper.circle(cx, cy, plasmid_radius);
			plasmid.attr("stroke", colors.plasmid);
			var title = seq_length + ' bp';
			if (plasmid_name != "") {
				title = plasmid_name + "\n\n" + title;
			}
			var plasmid_label = paper.text(cx, cy, title);
			plasmid_label.attr({"fill":      colors.plasmid,
								"font-size": plasmid_font_size, });

			for (var ang = 0; ang < 360; ang += 30) {
				draw_tic_mark(ang);
			}
		}

		// Move features that overlap to other radii.
		function resolve_conflicts() {
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
				// those features were not in conflict with the winner
				for (var pfx in loser.pushed_features) {
					var pf = loser.pushed_features[pfx];
					// Check for conflict with other the winner feature itself.
					// If there's no conflict, we can pushh it back safely.
					if (pf.start_degrees() - winner.end_degrees() <= min_overlap_cutoff ||
						winner.start_degrees() - pf.end_degrees() <= min_overlap_cutoff) {
						if (_debug)
							console.warn(pf.name() + " unpushed, because " 
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
				var furthest_point = plasmid_start; // Start at a complete lower bound

				// Go through the feature list twice, to make sure that features
				// that cross the boundary are resolved
				for (var fx = 0; fx < 2 * features.length; fx++) {
					var f = features[fx % features.length];
					if (f.radius == rad && f.type() != ft.enzyme) { 
						var new_size = f.size_degrees();
						var overlap = -(furthest_point - f.start_degrees());
						if (overlap <= min_overlap_cutoff) { 
							// We've cleared all potential conflicts: reset
							// the indicators
							biggest_size = new_size;
							biggest_feature = f;
							furthest_point = f.end_degrees();
						// since we go around twice, it is now possible
						// for a feature to "conflict with itself," so we
						// explicitly prevent this
						} else if ( !(biggest_feature === f) && 
								   biggest_size > min_overlap_feature_size &&
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

		function draw_features() {
			for (var fx in features) {
				features[fx].draw();
			}
		}

		// Make sure that the appropriate cutters are shown
		function show_hide_cutters() {
			for (var fx in features) {
				var f = features[fx];
                if (f.default_show_feature()) {
                    // Only draw enzymes if they are in the list of
                    // cutters to show - i.e. 1 cutter, 2 cutters,
                    // etc.
                    if (f.type() == ft.enzyme) {
					    if (cutters_to_show.indexOf(f.cut_count()) < 0) {
						    f.hide();
						    f.clear_label();
					    } else {
						    f.show();
						    f.show_label();
					    }
				    }
                }
                else {
                    // If the enzyme is not set to be shown by
                    // default, don't show it
                    f.hide();
                    f.clear_label();
                }
			}
		}

		// Global label list: keeps track of which label should be in each
		// label list, and where the label line should point to, so we can
		// use this information to figure out where to put the label and
		// minimize the number of lines that intersect.
		var label_f_c = new Array(8);
		var label_list_pos = new Array(8);

        function set_label_lists() {
            // Global: keeps track of feature centers for each label
            // list, we need this to compute exactly where a label
            // should be within a label list, so to minimize
            // intersecting lines.
			label_f_c = [[], [], [], [], [], [], [], []];
			for (var fx = features.length - 1; fx >= 0; fx--) {
				features[fx].set_label_list();
			}
        }

		function draw_labels(label_radius) {

            // lower x, y starting position for each label list
            label_list_pos = [[0,0], [0,0], [0,0], [0,0],
							  [0,0], [0,0], [0,0], [0,0]];
		
			// Iterate counterclockwise, first get counts
			// Sort feature center list for each label list, and also
			// figure out where each label list should start
			for (var i=0; i<label_f_c.length; i++) {
				label_f_c[i].sort(function(a,b){return (a[0]-b[0])})
				var section_angle = label_list_section_angle(i);
				var xy1 = convert.polar_to_rect(label_radius, section_angle);

                // for each section, we also shift the x coordinate to
                // be further away from the circle, so that as much as
                // possible, the angle between a) line from the label
                // to the feature and b) the label is more than 90
                // degrees (i.e. visually, you don't have lines going
                // "backward").
                //
                // we also compute the lower y coordinate of each
                // label list below.
				if (i == 0 || i == 1) {
					xy1.x += 60;
				}
				else if (i == 2 || i == 3) {
					xy1.y += label_f_c[i].length*label_letter_height;
					xy1.x += 60;
				}
				else if (i == 4 || i == 5) {
					xy1.y += label_f_c[i].length*label_letter_height;
					xy1.x -= 60;
				}
				else if (i == 6 || i == 7) {
					xy1.x -= 60;
				}
				label_list_pos[i][0] = xy1.x;
				label_list_pos[i][1] = xy1.y;
			}
	 
			// Finally draw labels
			for (var fx = features.length - 1; fx >= 0; fx--) {
				features[fx].draw_label();
			}
		}

		function extend_features() {
			for (var bfx = 0; bfx < basic_features.length; bfx++) {
				features.push(new CircularFeature(basic_features[bfx]));
			}
		}

		function draw() { // Draw the circular map
            // Extend basic features to get list of circular features
			extend_features();
            // Hide the right cutters
            show_hide_cutters();
            // Resolve conflicts on the circle, push some overlapping
            // features to other radii
			var max_radius = resolve_conflicts();
			var label_radius = max_radius + label_radius_offset; 
            // Determine which labels are in which lists
            set_label_lists();

            // Figure out outter edge of label lists
            //
            // Just an educated guess based on 13pt font. we will use
            // this to compute height of label lists. These are
            // conservative.
            label_letter_height = 15;
            label_letter_width = 12;
          
            var min_x = map_width/2;
            var max_x = map_width/2;
            var min_y = map_width/2;
            var max_y = map_width/2;
			var label_radius = max_radius + label_radius_offset; 
            for (var section=0; section<label_f_c.length; section++) {
                var list_max_letters = 0;
                for (var i=0; i<label_f_c[section].length; i++) {
                    if (label_f_c[section][i][1].length > list_max_letters) {
                        list_max_letters = label_f_c[section][i][1].length;
                    }
                }
                // +1 - i am not sure why we need it, but otherwise it
                // crops the last label
                var list_height = (label_f_c[section].length+1)*label_letter_height;
                var list_width = list_max_letters*label_letter_width;
				var section_angle = label_list_section_angle(section);
				var xy = convert.polar_to_rect(label_radius,section_angle);

				if (section == 0 || section == 1) {
                    // upper right
                    if (min_y > xy.y-list_height) { min_y = xy.y-list_height; }
                    if (max_x < xy.x+list_width) { max_x = xy.x+list_width; }
				}
				else if (section == 2 || section == 3) {
					// lower right
                    if (max_y < xy.y+list_height) { max_y = xy.y+list_height; }
                    if (max_x < xy.x+list_width) { max_x = xy.x+list_width; }
				}
				else if (section == 4 || section == 5) {
					// lower left
                    if (max_y < xy.y+list_height) { max_y = xy.y+list_height; }
                    if (min_x > xy.x-list_width) { min_x = xy.x-list_width; }
				}
				else if (section == 6 || section == 7) {
					// upper left
                    if (min_y > xy.y-list_height) { min_y = xy.y-list_height; }
                    if (min_x > xy.x-list_width) { min_x = xy.x-list_width; }
                }
            }

            // Now we have a new bounding box: min_x,min_y to max_x,max_y

            var right_x_extend = max_x-cx;
            var left_x_extend = cx-min_x;
            var top_y_extend = cy-min_y;
            var bot_y_extend = max_y-cy;
            var bb_width = max_x-min_x;
            var bb_height = max_y-min_y;

		    map_width = bb_width;
		    map_height = bb_height;
            cx = left_x_extend;
            cy = top_y_extend;

			paper = ScaleRaphael(map_dom_id, map_width, map_height); // global
            for (var fx in features) {
                features[fx].initialize();
            }

            // figure out the real height of labels
			var label = paper.text(0,0,'M');
			label.attr({"font-size": label_font_size});
			label_letter_height = label.getBBox().height; // global
			label_letter_width = label.getBBox().width; // global
			paper.clear();

			draw_plasmid();
			draw_features(); // Draw all the features initially
			draw_labels(label_radius); // Draw only the necessary labels

			// Rescale
			if (final_map_width != map_width ||
				final_map_height != map_height) {
				// "center" parameter just adds unnecessary CSS to the container
				// object to give it an absolute position: not what we need
				paper.changeSize(final_map_width,final_map_height,false,false)
			}
		}

		///////////////////////////////////////////////////////////////////
		// Main entry point.
		draw();

		// Export the main properties as part of the CircularMap object
		this.paper = paper;
		this.draw = draw;
		this.features = features;
	}; // End CircularMap()

	///////////////////////////////////////////////////////////////////
	// Linear Map Drawing Class
	this.LinearMap = function (options) {
	
		// Map-specific canvas element
		var paper;

		// Map-specific feature list
		var features = [];

		// Paper setup - not the final width, but how we will draw the
		// map, we will scale later on
		var map_width = 800;
		var map_height = 800;
		var cx = map_width/2;
		var cy = map_height/2;

		var plasmid_y = cy;
		var plasmid_width = map_width * 0.9;
		var plasmid_left = (map_width - plasmid_width) / 2;
		var plasmid_right = plasmid_left + plasmid_width;

		// Where to draw the map
		var map_dom_id = 'giraffe-draw-map';
		if ('map_dom_id' in options) {
			map_dom_id = options['map_dom_id'];
		}

		// Final size
		var final_map_width = 640;
		var final_map_height = 640;
		if ('map_width' in options) {
			final_map_width = parseInt(options['map_width'])
		}
		if ('map_height' in options) {
			final_map_height = parseInt(options['map_height'])
		}

		// Heights of levels
		var y_spacing = 20; // spacing
		var label_y_offset = 50;
		if ('label_offset' in options) {
			label_y_offset = parseInt(options['label_offset']);
		}

		// Feature visual properties
		var feature_width = 15;
		var enzyme_width = 25;
		var enzyme_weight = 1; // Restriction enzymes are drawn differently
							   // This controls their "thickness" on the map
		var enzyme_bold_weight = 3; // How thick when highlighted
		var feature_opacity = 0.7;
		var enzyme_opacity = 0.7;
		if ('opacity' in options) {
			feature_opacity = parseFloat(options['opacity']);
			enzyme_opacity = parseFloat(options['opacity']);
		}
		var bold_opacity = 1.0;
		var head_width = 25;
		var head_length = 5;

		// Cutters to show
		var cutters_to_show = [1];
		if ('cutters' in options) {
			cutters_to_show = options['cutters'];
		}

		// Animation properties
		var fade_time = 0;
		if ('fade_time' in options) {
			fade_time = parseInt(options['fade_time'])
		}

		// Overlaps
		var min_overlap_cutoff = -1;// in pixels
		var min_overlap_pct = 0;
		var min_overlap_feature_size = 0; // in pixels
		
		// Labels and other text
		var label_line_weight = 1;
		var label_line_bold_weight = 1.5 * label_line_weight;
		var label_font_size = '13pt';
		var plasmid_font_size = '16pt';
		var plasmid_name = '';
		if ('plasmid_name' in options) {
			plasmid_name = options['plasmid_name'];
		}

        var label_letter_height = 0;
        var label_letter_width = 0;

		var convert = {
			pos_to_x: function (p) {
				return plasmid_left + (p/seq_length) * plasmid_width;
			}
		};

		// TODO: MAJOR CODE REORGANIZATION: MERGE COMMON ELEMENTS INTO ONE CLASS
		function LinearFeature(basic_feature) {
			// Clone the basic feature, to extend it later
			function Clone() {}; // empty function to use as a hanger for prototype
			Clone.prototype = basic_feature;
			var _this = new Clone(); // Make a new "this" pointer to return at
			                         // the end
			// The result of this function will be a LinearFeature object

			// The visual object to modify when accessing the feature.
			var _feature_set;
			var _arrow_set;
			var _label_set;
			var _label_drawn = false;

			// Visual properties
			var _visible = true;
			var _labeled = true;

			// Type-based property selection
			var _color = colors.feature;
			var _width = feature_width;
			var _draw_head = false;
			var _opacity = feature_opacity;
			var _opaque = false; // holds opacity for clicks
			switch(_this.type()) {
				case ft.promoter:
				case ft.primer:
					_draw_head = true; // Promotors and primers are the only primer-colored
				case ft.terminator:   // features with heads
					_color = colors.primer;
					break;
				case ft.regulatory:
				case ft.origin:
					_color = colors.origin;
					break;
				case ft.enzyme:
					_color = colors.enzyme;
					_width = enzyme_width;
					_opacity = enzyme_opacity;
					break;
                case ft.orf:
                    _color = colors.orf;
				case ft.gene:
					_draw_head = true;
					break;
			}

			// Accessors for private properties set at creation
			_this.visible = function() { return _visible; };
			_this.labeled = function() { return _labeled; };

			_this.y = 0; // default to plasmid height (y is an offset)

			// Check to see if label has been drawn yet
			_this.label_drawn = function() { return _label_drawn; };
			_this.feature_set = function() { return _feature_set };
			_this.label_set = function() { return _label_set };

            _this.initialize = function() {
			    _feature_set = paper.set();
			    _arrow_set = paper.set();
			    _label_set = paper.set();
            }

			_this.real_center = function() {
				return (this.real_start() + this.real_end()) / 2.0;
			}

			// Degree conversion, for overlap calculation:
			// for these functions, the sequence starts at 90 degrees and goes down.
			_this.real_start = function() {
				var rs;
				// Take the minimum head size into account. Only need to do this 
				// when the head is drawn and pointing clockwise, to
				// "push the start back."
				if (_draw_head && this.clockwise()) { 
					rs = convert.pos_to_x(this.end()) - this.real_size();
				} else { // Headless feature, or head is pointing the wrong way.
						 // Just give its typical start position
					rs = convert.pos_to_x(this.start());
				}
				return rs;
			};

			_this.real_end = function() {
				var re;
				// Take the minimum head size into account. Only need to do this 
				// when the head is drawn and pointing counterclockwise, to 
				// "push the end forward."
				if (_draw_head && !this.clockwise()) { // Take the minimum head size into account
					re = convert.pos_to_x(this.start()) + this.real_size();
				} else { // Headless feature, or head is pointing the wrong way.
						 // Just give its typical end position
					re = convert.pos_to_x(this.end());
				}

				return re;
			};

			_this.real_size = function() {
				var rsz; 
				// Normal definition of size
				rsz = convert.pos_to_x(_this.end()) -
				      convert.pos_to_x(_this.start());

				// Head size: return this if it's bigger
				if (_draw_head && head_length > rsz)
						rsz = head_length;

				return rsz;
			};

			// Actions for interactivity
			// Generic fading animation/property setting mechanism
			var _fade = function (props, line_props) {
				var sets = paper.set();
				var lines = paper.set();
				sets.push(_feature_set);
				lines.push(_label_set[0]); // label line

                // Cutters: highlight other examples of this enzyme if
                // it's a multi-cutter
                if (_this.type() == ft.enzyme) {
					for (var fx in _this.other_cutters()) {
						var f = features[_this.other_cutters()[fx]];
						sets.push(f.feature_set());
						lines.push(f.label_set()[0]);
					}
				}

				if (fade_time) { 
					sets.animate(props, fade_time); 
					lines.animateWith(sets, line_props, fade_time); 
				} else { 
					sets.attr(props);
					lines.attr(line_props); 
				}
			}

			var _bolder = function () {
				var props = {"opacity": bold_opacity,
                             "font-weight": "bold" };
                if (_this.type() == ft.enzyme) {
				    props["stroke-width"] = enzyme_bold_weight;
                }
				var line_props = {"stroke": colors.plasmid, 
				                  "stroke-width": label_line_bold_weight};
				_fade(props, line_props);
			}

			var _lighter = function () {
				var props = {"opacity": _opacity,
                             "font-weight":"normal"};
                if (_this.type() == ft.enzyme) {
				    props["stroke-width"] = enzyme_weight;
                }
				var line_props = {"stroke": colors.bg_text,
				                  "stroke-width": label_line_weight};
				_fade(props, line_props);
			}

			// Toggle solid/light upon click
			var _click = function (event) {
				if (_opaque) {
					_lighter();
					_opaque = false;
				} else {
					_bolder();
					_opaque = true;
				}
			};

			// Hovering: solid/light upon mouseover
			var _mouse_over = function (event) {
				if (!_opaque)
					_bolder();
			};
			var _mouse_up = function (event) {
				if (!_opaque)
					_lighter();
			};

			_this.draw = function () {
				// Don't draw features that cross the boundary, as this is not
				// a circular plasmid
                if (!_visible || _this.crosses_boundary()) { return; }

				// Convert from sequence positions to x-coords
				var x0 = convert.pos_to_x(_this.start());
				var x1 = convert.pos_to_x(_this.end());

				var y = plasmid_y + _this.y;

				// Arrowhead drawing, if needed
				if (_draw_head) {
					var hx_tip, hx_back;
					if (_this.clockwise()) {
						hx_tip = x1;
						x1 -= head_length;
						hx_back = x1;
					} else {
						hx_tip = x0;
						x0 += head_length;
						hx_back = x0;
					}

					// Unlike the body, the head is traced with a line, and
					// then created entirely with the fill color
					var head = paper.path(svg.move(hx_tip, y) +
					                 svg.line(hx_back, y - head_width/2.0) +
					                 svg.line(hx_back, y + head_width/2.0) +
					                 svg.close());
					head.attr({"stroke-width": 0,
							   "fill":         _color});
					_arrow_set.push(head);
				}

				// Body drawing
				if (x0 < x1 && _this.type() != ft.enzyme) { 
					// Compensating for the head may have "taken up" all
					// the room on the plasmid, in which case no arc needs
					// to be drawn

					// The body has no fill-color: it's just a thick line
					var body = paper.path(svg.move(x0, y) +
						  				  svg.line(x1, y));
					body.attr({"stroke-width": _width});

					_arrow_set.push(body);
				} else if (_this.type() == ft.enzyme) { 
					// Restriction enzymes get drawn on their own
					var x_m = (x0 + x1)/2;

					var body = paper.path(svg.move(x_m, y - _width/2.0) +
					                      svg.line(x_m, y + _width/2.0));
					body.attr({"stroke-width": enzyme_weight});
					body.toBack();

					_arrow_set.push(body);
				}

				_arrow_set.click(_click);
				_arrow_set.hover(_mouse_over, _mouse_up);

				_feature_set.push(_arrow_set);

				// Apply the feature-wide properties to the whole feature
				_feature_set.attr({"stroke":         _color,
								   "stroke-linecap": "butt",
								   "opacity":        _opacity,
								   "title":          _this.name()});

			} // END LinearFeature::draw()

			// Should we draw the label?
			_this.should_draw_label = function () {
				// Don't bother unless we need to
				if (!_visible || !_labeled || this.crosses_boundary()) 
					return false;
				return true;
			} // END LinearFeature::should_draw_label()

            // What label should we draw?
            _this.label_name = function () {
				var label_name = _this.name();
				if (_this.type() == ft.enzyme) {
					label_name += " (" + _this.cut() + ")";
				}
                return label_name;
            } // END LinearFeature::label_name()

			_this.label_size = function() {
				var fake_label = paper.text(0, 0, this.label_name());
				var w = fake_label.getBBox().width;
				var h = fake_label.getBBox().height;
				fake_label.remove();
				
				return { width: w, height: h };
			}

			_this.draw_label = function (height, pos) {
				if (!this.should_draw_label()) { return; }

				if (_label_drawn)
					this.clear_label();

				// Figure out the center of the feature
				var x_c = this.real_center();

				// Enzymes show their cut sites in the label
				var label_name = _this.label_name();
				var label = paper.text(pos, height, label_name);
				
				// Below, right-justify. Above, left-justify.
				var anchor = (height >= plasmid_y) ? "end" : "start";
				label.attr({"fill": _color,
							"text-anchor": anchor,
							"font-size": label_font_size,
							"opacity": 1.0 });

				// Draw the line to the label position
				var label_line = paper.path(svg.move(x_c, plasmid_y + this.y) +
											svg.line(pos, height));
				label_line.attr({"stroke": colors.bg_text,
				                 "stroke-width": label_line_weight,
								 "opacity": 0.5 });

				_label_set.push(label_line);
				_label_set.push(label);

				// Handlers
				_label_set.click(_click);
				_label_set.hover(_mouse_over, _mouse_up);

                // Only push label_line, so when we fade in and out,
                // we don't also fade the label.
				_feature_set.push(label_line);

				_labeled = true;
				_label_drawn = true;

				return label.getBBox();
			} // END LinearFeature::draw_label()

			_this.hide = function () {
				if (_visible) {
                    if (_feature_set) { _feature_set.hide(); }
					_visible = false;
					_labeled = false;
				}
			}; // END LinearFeature::hide()

			_this.show = function () {
				if (!_visible) {
					if (_feature_set) { _feature_set.show(); }
					if (!_labeled) {
                        if (_label_set) { _label_set.hide(); }
                    }
					_visible = true;
				}
			}; // END LinearFeature::show()

			_this.hide_label = function () {
				if (_labeled) {
                    if (_label_set) { _label_set.hide(); }
					_labeled = false;
				}
			}; // END LinearFeature::hide_label()

			_this.show_label = function () {
				if (!_labeled) {
                    if (_label_set) { _label_set.show(); }
					_labeled = true;
				}
			}; // END LinearFeature::show_label()

			_this.clear_label = function () {
				if (_label_drawn) {
                    if (_label_set) {
					    _label_set.unclick(_click);
					    _label_set.unhover(_mouse_over, _mouse_up);
					    _label_set.remove();
					    _label_set = paper.set();
                    }
					_labeled = false;
					_label_drawn = false;
				}
			}; // END LinearFeature::clear_label()

			return _this;
		}; // END LinearFeature Class

		function draw_plasmid() {

			// Title
			var title_y = 1.5 * y_spacing;

			// Tic marks
			var tic_mark_length = 15;
			var tic_mark_y = plasmid_y + 2* y_spacing;
			var tic_label_y = tic_mark_y + 1.5*tic_mark_length;

			function draw_tic_mark(p) {
				var x = convert.pos_to_x(p);
				var y0 = tic_mark_y - tic_mark_length/2;
				var y1 = tic_mark_y + tic_mark_length/2;
				var tic = paper.path(svg.move(x, y0) +
									 svg.line(x, y1));
				tic.attr({"stroke": colors.bg_text});

				var label = paper.text(x, tic_label_y, String(p));
				label.attr({"fill": colors.bg_text});
			}

			var plasmid = paper.path(svg.move(plasmid_left,  plasmid_y) +
									 svg.line(plasmid_right, plasmid_y));

			plasmid.attr("stroke", colors.plasmid);
			var title = seq_length + ' bp';
			if (plasmid_name != "") {
				title = plasmid_name + ": " + title;
			}
			var plasmid_label = paper.text(cx, title_y, title);
			plasmid_label.attr({"fill":      colors.plasmid,
								"font-size": plasmid_font_size, });

			// Set the scale to be the order of magnitude of seq_length
			// i.e. 100, 1000, 10, etc.
			var scale = Math.pow(10, Math.floor(Math.log(seq_length)/Math.log(10)))
			for (var xx = scale; xx <= seq_length; xx += scale) {
				draw_tic_mark(xx);
			}

		}

		function resolve_conflicts() {
			var conflicts;
			var y = 0; // current radius
			var yx = 1;               // radius counter
			var max_dist = 0;

			function push(winner, loser) {
				// Record that the push happened
				winner.pushed_features.push(loser); 
				conflicts++;

				// Do it
				loser.y = new_y; 

				if (_debug) console.warn(loser.name() + " pushed by " + winner.name());
				
				// Since loser was pushed, un-push all the 
				// features it caused to be pushed, as long as
				// those features were not in conflict with the winner
				for (var pfx in loser.pushed_features) {
					var pf = loser.pushed_features[pfx];
					// Check for conflict with other the winner feature itself.
					// If there's no conflict, we can pushh it back safely.
					if (winner.real_end() - pf.real_start() <= min_overlap_cutoff ||
						pf.real_end() - winner.real_start() <= min_overlap_cutoff) {
						if (_debug)
							console.warn(pf.name() + " unpushed, because " 
								+ loser.name() + " pushed by " + winner.name());
						pf.y = y;
					}
				}
			}

			do {
				// Keep alternating between inside and outside the plasmid.
				var new_y = y + Math.pow(-1, yx) * yx * y_spacing;

				conflicts = 0; // Assume you have no conflicts until you find some

				// Clear the record of who pushed whom
				for (var fx in features) {
					features[fx].pushed_features = [];
				}

				var biggest_size = 0;
				var biggest_feature;
				var furthest_point = plasmid_left; // Start at a complete lower bound

				for (var fx = 0; fx < features.length; fx++) {
					var f = features[fx];
					if (f.y == y && f.type() != ft.enzyme) { 
						var new_size = f.real_size();
						var overlap = furthest_point - f.real_start();
						if (overlap <= min_overlap_cutoff) { 
							// We've cleared all potential conflicts: reset
							// the indicators
							biggest_size = new_size;
							biggest_feature = f;
							furthest_point = f.real_end();
						// explicitly prevent conflicts with self
						} else if ( !(biggest_feature === f) && 
								   biggest_size > min_overlap_feature_size &&
								   new_size > min_overlap_feature_size &&
								  (overlap <= 0 || 
								  (overlap/biggest_size > min_overlap_pct &&
								   overlap/new_size > min_overlap_pct))) {
							// Overlap: conflict!
							if (new_size > biggest_size) { // This feature is top dog,
														   // move the original to the
														   // new height
								push(f, biggest_feature);

								// Update the new top dog
								biggest_size = new_size;
								biggest_feature = f;
								furthest_point = f.real_end();

							} else { // The original feature is top dog. move the new
									 // feature to the new height

								push(biggest_feature, f);
							}

						}
					}
				}

				// Keep track of the biggest distance from the plasmid 
				// reached
				var new_dist = Math.abs(y);
				if (new_dist > max_dist)
					max_dist = new_dist;

				// Move on to the next radius
				y = new_y;
				yx++;
				
			} while (conflicts > 0); // Keep adding levels of resolution

			return max_dist;
		}

		function extend_features() {
			for (var bfx = 0; bfx < basic_features.length; bfx++) {
				features.push(new LinearFeature(basic_features[bfx]));
			}
		}

		// XXX copy-and-pasted: eliminate redundancy
		function draw_features() {
			for (var fx in features) {
				features[fx].draw();
			}
		}

		// Make sure that the appropriate cutters are shown
		function show_hide_cutters() {
			for (var fx in features) {
				var f = features[fx];
                if (f.default_show_feature()) {
                    // Only draw enzymes if they are in the list of
                    // cutters to show - i.e. 1 cutter, 2 cutters,
                    // etc.
                    if (f.type() == ft.enzyme) {
					    if (cutters_to_show.indexOf(f.cut_count()) < 0) {
						    f.hide();
						    f.clear_label();
					    } else {
						    f.show();
						    f.show_label();
					    }
				    }
                } else {
                    // If the enzyme is not set to be shown by
                    // default, don't show it
                    f.hide();
                    f.clear_label();
                }
			}
		}

		var label_pos, label_lists;
		function assign_label_lists() {
			var label_overlap_cutoff = -1; // pixel

			var nlists = 6;
			//                   top                bottom
            label_pos    = [ new Array(nlists), new Array(nlists)];
            label_lists  = [ new Array(nlists), new Array(nlists)];
                                
			// Calculate positions of label lists
			var list_spread = plasmid_width / (nlists + 1);
			//                 top  bottom
			var offset_frac = [1, -1];
			for (var ix = 0; ix < nlists; ix++) {
				for (var lx = 0; lx < 2; lx++) {
					// offset to the right or left by offset_frac
					// depending on if it's above or below
					label_pos[lx][ix] = (ix + 1 + offset_frac[lx]) * list_spread;
					label_lists[lx][ix] = [];
				}
			}

			// Figure out which list each feature goes in
            for (var fx in features) {
				var f = features[fx];

				// Which quarter of the plasmid is the feature in?
				var section = Math.floor(nlists*(f.real_center() - plasmid_left)/
				                                 plasmid_width);
				// Is it in the top or bottom?
				var bottom = section % 2;

				if (f.should_draw_label()) {
					// push it on the appropriate list
					label_lists[bottom][section].push(f);
				}
            }

		}

		function draw_labels(height) {
			// Calculate the height of a label
			var label_leading = 1.3;
			var label_height = features[0].label_size().height * label_leading;

			// Finally, draw all the features
            for (var sx = 0; sx < label_pos[0].length; sx++) {
				for (var lx = 0; lx < 2; lx++) {
					var ll = label_lists[lx][sx];

					// Sort the list by center position
					var comp_factor = 0.75;
					ll.sort(function (a,b) {
						// Make some compensation for height as well
						function key(feat) {
							// On top, heigts closer to you are 
							// negative, so we flip the sign depending on
							// whether or not lx is 0 (top) or 1 (bottom)
							return feat.real_center() + 
								comp_factor * feat.y;
						}
						return key(a) - key(b);
					})

					var num_labels = ll.length;

					// Iterate over every label in the list
					for (var ix = 0; ix < num_labels; ix++) {
						var curr_height;
						if (lx) { // Bottom list: top to bottom
							curr_height = plasmid_y + height + 
								(ix) * label_height;
						} else { // Top list: bottom to top
							curr_height = plasmid_y - height - 
								(num_labels - 1 - ix) * label_height;
						}
						ll[ix].draw_label(curr_height, label_pos[lx][sx]);
					}
				}
			}
		}

		function set_bounding_box(height) {
            // Figure out outer edge of label lists
            //
            // Just an educated guess based on 13pt font. we will use
            // this to compute height of label lists. These are
            // conservative.
            label_letter_height = 15;
            label_letter_width = 8;
          
            var min_y = map_height/2;
            var max_y = map_height/2;
			// By default, plasmid will scale to width of map, so unless we
			// actually have lists that go off the page, no reason to adjust
			// them. i.e., never make them smaller than they orignially were,
			// because there is never a need to "zoom in."
            var min_x = 0;
            var max_x = map_width;

			// Iterate over every list in a level
            for (var sx = 0; sx < label_pos[0].length; sx++) {
				for (var lx = 0; lx < 2; lx++) {
					var ll = label_lists[lx][sx];

					var list_max_letters = 0;
					for (var ix = 0; ix < ll.length; ix++) {
						var num_letts = ll[ix].label_name().length; 
						if (num_letts > list_max_letters) {
							list_max_letters = num_letts;
						}
					}

					if (lx == 0) { // Top lists: move top and right
						var list_top = plasmid_y - height - 
							label_letter_height * (ll.length + 1);
						if (list_top < min_y)
							min_y = list_top;
						var list_right =  label_pos[lx][sx] + 
							label_letter_width * list_max_letters;
						if (list_right > max_x)
							max_x = list_right;

					} else if (lx == 1) { // Bot lists: move bot and left
						var list_bot = plasmid_y + height + 
							label_letter_height * (ll.length + 1);
						if (list_bot > max_y)
							max_y = list_bot;
						var list_left =  label_pos[lx][sx] - 
							label_letter_width * list_max_letters;
						if (list_left < min_x)
							min_x = list_left;
					}
				}
			}

            // Now we have a new bounding box (height only): min_y to max_y
			
			// Extend or compress the box dimensions to encompas this new size
		    map_width = max_x - min_x;
		    map_height = max_y - min_y;

			// Shift all the reference points to compensate for the re-zooming
            cy -= min_y;
            cx -= min_x;
			plasmid_y = cy;
			plasmid_left -= min_x;
			plasmid_right -= min_x;

			// Shift the label positions as well, to compensate
			for (var lx = 0; lx < label_pos.length; lx++) {
				for (var ix = 0; ix < label_pos[lx].length; ix++) {
					label_pos[lx][ix] -= min_x;
				}
			}

		}

		function draw() { // Draw the linear map
            // Extend basic features to get list of linear features
			extend_features();

            // Hide the right cutters
            show_hide_cutters();
            // Resolve conflicts on the line, push some overlapping
            // features to other radii
			var max_height = resolve_conflicts();
			var label_height = max_height + label_y_offset; 

			assign_label_lists();
			set_bounding_box(label_height);
        
			paper = ScaleRaphael(map_dom_id, map_width, map_height); // global
            for (var fx in features) {
                features[fx].initialize();
            }

			draw_plasmid();
			draw_features(); // Draw all the features initially
			draw_labels(label_height); // Draw only the necessary labels

			// Rescale
			if (final_map_width != map_width ||
				final_map_height != map_height) {
				
				// Make sure not to add additional height to the map, once we've
				// trimmed it off
				final_map_height = final_map_width * (map_height/map_width);

				// "center" parameter just adds unnecessary CSS to the container
				// object to give it an absolute position: not what we need
				paper.changeSize(final_map_width,final_map_height,false,true)
			}
		}
		
		
		///////////////////////////////////////////////////////////////////
		// Main entry point.
		draw();

		// Export the main properties as part of the LinearMap object
		this.paper = paper;
		this.draw = draw;
		this.features = features;
	}; // End LinearMap()

	return this;
}})();
