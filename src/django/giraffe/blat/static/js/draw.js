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
(function(){
 	// XXX: ONLY PROTOTYPAL INHERITANCE 


	///////////////////////////////////////////////////////////////////
	// Prototypal inheritance utilities
	// h/t: crockford
	if (typeof Object.extend !== 'function') {
		Object.extend = function (original) {
			function Extender() {}
			Extender.prototype = o;
			return new Extender();
		};
	}

	///////////////////////////////////////////////////////////////////
	// Package-scope settings
	var _debug = false;

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
		orf:     "#00c8c8"
	};

	// Sequence data
	var sequence = {
		length: 0,
		features: []
	};

	///////////////////////////////////////////////////////////////////
	// BasicFeature class (private, internal class)
	// holds only the most essential information about a feature
	function BasicFeature(feat) {

		// Private data members, from the properties of the feature argument
		// since hash tables and objects are one and the same in js
		this.name = feat.feature;
		this.start = parseInt(feat.start);
		this.end = parseInt(feat.end);
		this.type = parseInt(feat.type_id);
        this.default_show_feature = 1;
		this.clockwise = feat.clockwise;
		this.cut = parseInt(feat.cut); // only for enzymes;

		this.other_cutters = []; // only for enzymes;

        // Not all features have this option, e.g. annotated ones we
        // have to show the feature.
        if ('show_feature' in feat) { 
			this.default_show_feature = parseInt(feat.show_feature); 
		}

		// Enzyme-only data access methods
		this.crosses_boundary = function () { return this.end < this.start };
		this.cut_count = function() { return this.other_cutters.length; }; // gives 0 if not enzyme

	} // End BasicFeature class

	///////////////////////////////////////////////////////////////////
	// Feature class (private, internal class)
	// Holds elements and properties common to all features which
	// are drawn
	function Feature(basic_feature) {

		var _this = Object.extend(basic_feature);

		// Type-based property selection
		_this.color = colors.feature;
		_this.width = Feature.feature_width;
		_this.draw_head = false;

		// Type-based property selection
		_this.opacity = Feature.feature_opacity;
		_this.opaque = false; // holds opacity for clicks

		switch(_this.type) {
			case ft.promoter:
			case ft.primer:
				_this.draw_head = true; // Promotors and primers are the only primer-colored
			case ft.terminator:   // features with heads
				_this.color = colors.primer;
				break;
			case ft.regulatory:
			case ft.origin:
				_this.color = colors.origin;
				break;
			case ft.enzyme:
				_this.color = colors.enzyme;
				_this.width = Feature.enzyme_width;
				_this.opacity = Feature.enzyme_opacity;
				break;
			case ft.orf:
				_this.color = colors.orf;
			case ft.gene:
				_this.draw_head = true;
				break;
		}

		// Radius or height at which to draw feature
		_this.position = 0;

		// The visual object to modify when accessing the feature.
		_this.feature_set = {};
		_this.arrow_set = {};
		_this.label_set = {};
		_this.label_drawn = false;

		// Actions for interactivity
		// Generic fading animation/property setting mechanism
		_this.fade = function (props, line_props) {
			var sets = this.paper.set();
			var lines = this.paper.set();
			sets.push(this.feature_set);
			lines.push(this.label_set[0]); // label line

			// Cutters: highlight other examples of this enzyme if
			// it's a multi-cutter
			if (this.type() == ft.enzyme) {
				for (var fx in this.other_cutters()) {
					var f = features[this.other_cutters()[fx]];
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

		_this.bolder = function () {
			var props = {"opacity": bold_opacity,
						 "font-weight": "bold" };
			if (this.type() == ft.enzyme) {
				props["stroke-width"] = enzyme_bold_weight;
			}
			var line_props = {"stroke": colors.plasmid, 
							  "stroke-width": label_line_bold_weight};
			this.fade(props, line_props);
		}

		_this.lighter = function () {
			var props = {"opacity": this.opacity,
						 "font-weight":"normal"};
			if (this.type() == ft.enzyme) {
				props["stroke-width"] = enzyme_weight;
			}
			var line_props = {"stroke": colors.bg_text,
							  "stroke-width": label_line_weight};
			this.fade(props, line_props);
		}

		// Toggle solid/light upon click
		_this.click = function (event) {
			if (_opaque) {
				this.lighter();
				this.opaque = false;
			} else {
				this.bolder();
				this.opaque = true;
			}
		};

		// Hovering: solid/light upon mouseover
		_this.mouse_over = function (event) {
			if (!this.opaque)
				this.bolder();
		}; // END Feature::mouse_over()

		_this.mouse_up = function (event) {
			if (!this.opaque)
				this.lighter();
		}; // END Feature::mouse_up()

		_this.hide = function () {
			if (this.visible) {
				if (this.feature_set) { this.feature_set.hide(); }
				this.visible = false;
				this.labeled = false;
			}
		}; // END Feature::hide()

		_this.show = function () {
			if (!this.visible) {
				if (this.feature_set) { this.feature_set.show(); }
				if (!this.labeled) {
					if (this.label_set) { this.label_set.hide(); }
				}
				this.visible = true;
			}
		}; // END Feature::show()

		_this.hide_label = function () {
			if (this.labeled) {
				if (this.label_set) { this.label_set.hide(); }
				this.labeled = false;
			}
		}; // END Feature::hide_label()

		_this.show_label = function () {
			if (!this.labeled) {
				if (this.label_set) { this.label_set.show(); }
				this.labeled = true;
			}
		}; // END Feature::show_label()

		_this.clear_label = function () {
			if (this.label_drawn) {
				if (this.label_set) {
					this.label_set.unclick(this.click);
					this.label_set.unhover(this.mouse_over, this.mouse_up);
					this.label_set.remove();
					this.label_set = paper.set();
				}
				this.labeled = false;
				this.label_drawn = false;
			}
		}; // END Feature::clear_label()

		// Should we draw the label?
		_this.should_draw_label = function () {
			// Don't bother unless we need to
			if (!this.visible || !this.labeled) 
				return false;
			return true;
		} // END Feature::should_draw_label()

		// What label should we draw?
		_this.label_name = function () {
			var label_name = _this.name();
			if (this.type == ft.enzyme) {
				label_name += " (" + this.cut + ")";
			}
			return label_name;
		} // END Feature::label_name()

		return _this;
	}; // END Feature Class

	// Feature property defaults are static members of the Feature class
	Feature.feature_width = 15;
	Feature.enzyme_width = 25;

	Feature.enzyme_weight = 1; // Restriction enzymes are drawn differently
							   // This controls their "thickness" on the map
	Feature.enzyme_bold_weight = 3; // How thick when highlighted
	Feature.feature_opacity = 0.7;
	Feature.enzyme_opacity = 0.7;
	Feature.bold_opacity = 1.0;
	Feature.head_width = 25;
	Feature.head_length = 7;


	///////////////////////////////////////////////////////////////////
	// Circular Feature class
	function CircularFeature(feature, map) {
		var _this = Object.extend(feature);

		// Default to plasmid radius, can be changed later on by other methods
		_this.position = map.plasmid_pos; 

		_this.feature_set = map.paper.set();
		_this.arrow_set = map.paper.set();
		_this.label_set = map.paper.set();

		// Calculated properties

		// Degree conversion, for overlap calculation:
		// for these functions, the sequence starts at 90 degrees and goes down.
		_this.start_degrees = function() {
			var sd;
			// Take the minimum head size into account. Only need to do this 
			// when the head is drawn and pointing clockwise, to
			// "push the start back."
			if (this.draw_head && this.clockwise) { 
				sd = map.convert.pos_to_angle(this.end) + this.size_degrees();
			} else { // Headless feature, or head is pointing the wrong way.
					 // Just give its typical start position
				sd = map.convert.pos_to_angle(this.start);
			}
			return sd;
		};

		_this.end_degrees = function() {
			var ed;
			// Take the minimum head size into account. Only need to do this 
			// when the head is drawn and pointing counterclockwise, to 
			// "push the end forward."
			if (this.draw_head && !this.clockwise) { // Take the minimum head size into account
				ed = map.convert.pos_to_angle(this.start) - this.size_degrees();
			} else { // Headless feature, or head is pointing the wrong way.
					 // Just give its typical end position
				ed = map.convert.pos_to_angle(this.end);
			}

			// For boundaries that cross the 0 point, ensure that they return
			// a number in the same range as everyone else
			if (this.crosses_boundary()) 
				ed += 360;

			return ed;
		};

		_this.size_degrees = function() {
			var szd; // size in degrees
			// Normal definition of size
			if (this.crosses_boundary()) 
				// Start and end are flipped here: non-intuitive
				szd = map.convert.seq_length_to_angle(sequence.length - this.start + this.end + 1);
			else
				szd = map.convert.seq_length_to_angle(this.end - this.start + 1);

			// Head size: return this if it's bigger
			if (this.draw_head) {
				// Convert the head length into degrees, just as you do
				// in the draw() method. Must recalcualte every time, as
				// radius may have changed
				var r_p = Math.sqrt(this.position * this.position + 
						Feature.head_length*Feature.head_length);
				var hszd = Raphael.deg(Math.asin(Feature.head_length/r_p));
				if (hszd > szd)
					szd = hszd;
			}

			return szd;
		};

		// Feature drawing
		_this.draw = function () {
			if (!_visible) { return; }

			// Convert from sequence positions to angles
			var a0 = map.convert.pos_to_angle(_this.start);
			var a1 = map.convert.pos_to_angle(_this.end);

			// Arrowhead drawing, if needed
			if (_this.draw_head) {

				// Arrow tip point lines up with a0 or a1 and it points
				// tangent to the circle.
				// We need to figure out how many radians the arrow takes up
				// in order to adjust a0 or a1 by that amount, and to set the 
				// base of the triangle even with that angle
				var r_p = Math.sqrt(_this.radius*_this.radius + 
						Feature.head_length*Feature.head_length);
				// "height" of the arrowhead, in degrees
				var a_b;
				var a_p = Raphael.deg(Math.asin(Feature.head_length/r_p));
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
				var xy_p = map.convert.polar_to_rect(_this.radius, a_p);
				
				// bottom and top points, rectangular
				var xy_b = map.convert.polar_to_rect(_this.radius - Feature.head_width/2.0, a_b);
				var xy_t = map.convert.polar_to_rect(_this.radius + Feature.head_width/2.0, a_b);

				// Unlike the arc, the head is traced with a line, and
				// then created entirely with the fill color
				var head = paper.path(svg.move(xy_p.x, xy_p.y) +
									  svg.line(xy_b.x, xy_b.y) +
									  svg.line(xy_t.x, xy_t.y) + 
									  svg.close());
				head.attr({"stroke-width": 0,
						   "fill":         _this.color()});
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
				var xy0 = map.convert.polar_to_rect(_this.radius, a1);
				var xy1 = map.convert.polar_to_rect(_this.radius, a0);

				// The arc has no fill-color: it's just a thick line
				var arc = paper.path(svg.move(xy0.x, xy0.y) +
									 svg.arc(_this.radius, xy1.x, xy1.y));
				arc.attr({"stroke-width": _this.width()});

				_arrow_set.push(arc);
			} else if (_this.type() == ft.enzyme) { 
				// Restriction enzymes get drawn on their own
				var xy0 = map.convert.polar_to_rect(_this.radius - _this.width()/2.0, 
						(a0+a1)/2.0);
				var xy1 = map.convert.polar_to_rect(_this.radius + _this.width()/2.0, 
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
			_feature_set.attr({"stroke":         _this.color(),
							   "stroke-linecap": "butt",
							   "opacity":        _opacity,
							   "title":          _this.name()});

		} // END CircularFeature::draw()

		// Set which label list this feature's label should be in
		_this.set_label_list = function () {
			if (!_this.should_draw_label()) { return; }

			// Figure out the center of the feature
			var a_c = (_this.start_degrees()+_this.end_degrees())/2.0;
			var adjust_a_c = a_c;
			if (adjust_a_c < 0) { adjust_a_c += 360; }

			// Figure out which section this label is in: divide
			// the grid up into eight sections.
			var section = Math.floor((_this.plasmid_start-a_c)/label_section_degree);

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

			var xy0 = map.convert.polar_to_rect(_this.radius, a_c);
			
			// Figure out which section this label is in: divide
			// the grid up into eight sections.
			var section = Math.floor((_this.plasmid_start - a_c)/label_section_degree);
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

			label.attr({"fill": _this.color(),
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

		return _this;
	}; // END CircularFeature Derived Class


	///////////////////////////////////////////////////////////////////
	// Map Class 
	function Map(options) {
		// Map-specific canvas element
		this.paper = {};

		// Map-specific feature list
		this.features = [];

		// Paper setup - not the final width, but how we will draw the
		// map, we will scale later on
		this.width = 800;
		this.height = 800;
		this.cx = this.width/2;
		this.cy = this.height/2;

		// Where to draw the map
		this.dom_id = 'giraffe-draw-map';
		if ('map_dom_id' in options) {
			this.dom_id = options['map_dom_id'];
		}

		// Final size
		this.final_width = 640;
		this.final_height = 640;
		if ('map_width' in options) {
			this.final_width = parseInt(options['map_width'])
		}
		if ('map_height' in options) {
			this.final_height = parseInt(options['map_height'])
		}

		// Global plasmid info
		this.plasmid_start;
		this.plasmid_end;

		// Levels
		this.spacing = 20;
		this.plasmid_pos;
		this.label_offset = 10;
		if ('label_offset' in options) {
			this.label_offset = parseInt(options['label_offset']);
		}

		// Cutters to show
		this.cutters_to_show = [1];
		if ('cutters' in options) {
			this.cuters_to_show = options['cutters'];
		}

		// Animation properties
		this.fade_time = 0;
		if ('fade_time' in options) {
			fade_time = parseInt(options['fade_time'])
		}

		// Tic marks
		this.tic_mark_length = 15;

		// Labels and other text
		this.label_line_weight = 1;
		this.label_line_bold_weight = 1.5 * label_line_weight;
		this.label_font_size = '13pt';
		this.plasmid_font_size = '16pt';
		this.plasmid_name = '';
		if ('plasmid_name' in options) {
			plasmid_name = options['plasmid_name'];
		}

		if ('opacity' in options) {
			Feature.feature_opacity = parseFloat(options['opacity']);
			Feature.enzyme_opacity = parseFloat(options['opacity']);
		}

        this.label_letter_height = 0;
        this.label_letter_width = 0;

		// Save the options
		this.options = options;


		///////////////////////////////////////////////////////////////////
		// Methods
		this.extend_basic_features = function (basic_features) {
			for (var bfx = 0; bfx < basic_features.length; bfx++) {
				var extFeat = Feature(basic_features[bfx]);
				// Feature-specific option
				this.features.push(extFeat);
			}
		}

		this.draw_features = function () {
			for (var fx in this.features) {
				this.features[fx].draw();
			}
		}

		// Move features that overlap to other radii.
		this.resolve_conflicts = function () {
			// Overlaps
			var min_overlap_cutoff = -0.1; // in degrees
			var min_overlap_pct = 0.01;
			var min_overlap_feature_size = 0.5; // in degrees
		
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
							console.warn(pf.name() + "unpushed, because " 
								+ loser.name() + " pushed by " + winner.name());
						pf.radius = rad;
					}
				}
			}

			do {
				// Keep alternating between inside and outside the plasmid.
				var new_rad = rad + Math.pow(-1, rx) * rx * _this.spacing;

				conflicts = 0; // Assume you have no conflicts until you find some

				// Clear the record of who pushed whom
				for (var fx in features) {
					features[fx].pushed_features = [];
				}

				var biggest_size = 0;
				var biggest_feature;
				var furthest_point = _this.plasmid_start; // Start at a complete lower bound

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

		// Make sure that the appropriate cutters are shown
		this.show_hide_cutters = function () {
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


	} // End Map Class()

	///////////////////////////////////////////////////////////////////
	// Circular Map Drawing Class
	this.CircularMap = function(options) {

		// Extend a new base Map object
		var _this = Object.extend(new Map(options))

		// Global plasmid info
		_this.plasmid_start = 90; // degrees
		_this.plasmid_end = _this.plasmid_start - 360; // degrees

		// Conversion Object: local to this function because it assumes a 
		// circular geometry
		// Groups conversions into handy section
		_this.convert = {
			pos_to_angle: function (p) {
				//     start at the top of the circle
				return _this.plasmid_start - (p/sequence.length ) * 360;
			},
			seq_length_to_angle: function (l) {
				//     just like pos_to_angle, but without caring about the start
				return (l/sequence.length ) * 360;
			},
			angle_to_pos: function (a) {
				//     start at the top of the circle
				return Math.round(1 + ((sequence.length - 1)/360) * 
						((360 + _this.plasmid_start - a) % 360));
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
				rect.x += _this.cx;
				// Coordinates increase as you go down, so y is flipped.
                rect.y = _this.cy - rect.y;
				return rect;
			}
		};

		// This is based on using 8 label lists, which is pretty much hard
		// coded in a few places, so don't change this unless you figure
		// out how to change the number of label lists.
		var label_section_degree = 45;

        function label_list_section_angle(section) {
		    return _this.plasmid_start - label_section_degree/2.0 - 
				section * label_section_degree;
        };

		///////////////////////////////////////////////////////////////////
		// Drawing utility functions

		// Circle setup
		function draw_plasmid() {
			function draw_tic_mark(a) {
				var tic_mark_radius = _this.plasmid_pos - _this.spacing - _this.tic_mark_length/2;
				var tic_label_radius = tic_mark_radius - 1.5*_this.tic_mark_length;

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
				if (a < _this.plasmid_start || a > 360 - _this.plasmid_start) { // Right half of wheel: align right
					label.attr({"text-anchor": "end"});
				} else if (a > _this.plasmid_start && a < 360 - _this.plasmid_start) { // Left half of wheel: align left
					label.attr({"text-anchor": "start"});
				} // Top and bottom default to middle, which is correct
			}

			var plasmid = paper.circle(cx, cy, plasmid_radius);
			plasmid.attr("stroke", colors.plasmid);
			var title = sequence.length + ' bp';
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
          
            var min_x = _this.width/2;
            var max_x = _this.width/2;
            var min_y = _this.width/2;
            var max_y = _this.width/2;
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

		    _this.width = bb_width;
		    _this.height = bb_height;
            cx = left_x_extend;
            cy = top_y_extend;

			paper = ScaleRaphael(_this.dom_id, _this.width, _this.height); // global
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
			if (_this.final_width != _this.width ||
				_this.final_height != _this.height) {
				// "center" parameter just adds unnecessary CSS to the container
				// object to give it an absolute position: not what we need
				paper.changeSize(_this.final_width,_this.final_height,false,false)
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
	
		var plasmid_y = cy;
		var plasmid_width = _this.width * 0.9;
		var plasmid_left = (_this.width - plasmid_width) / 2;
		var plasmid_right = plasmid_left + plasmid_width;

		// Heights of levels
		var y_spacing = 20; // spacing
		var label_y_offset = 10;
		if ('label_offset' in options) {
			label_y_offset = parseInt(options['label_offset']);
		}

		// Feature visual properties
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

		// Tic marks
		var tic_mark_length = 15;
		// Title
		var title_y = 3 * y_spacing;

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
				return plasmid_left + (p/sequence.length ) * plasmid_width;
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
			var _opacity = feature_opacity;
			var _opaque = false; // holds opacity for clicks
			switch(_this.type()) {
				case ft.enzyme:
					_opacity = enzyme_opacity;
					break;
			}

			_this.y = plasmid_y; // default to plasmid height

			// Check to see if label has been drawn yet
			_this.label_drawn = function() { return _label_drawn; };
			_this.feature_set = function() { return _feature_set };
			_this.label_set = function() { return _label_set };

            _this.initialize = function() {
			    _feature_set = paper.set();
			    _arrow_set = paper.set();
			    _label_set = paper.set();
            }

			_this.draw = function () {
				// Don't draw features that cross the boundary, as this is not
				// a circular plasmid
                if (!_visible || _this.crosses_boundary()) { return; }

				// Convert from sequence positions to x-coords
				var x0 = convert.pos_to_x(_this.start());
				var x1 = convert.pos_to_x(_this.end());

				// Arrowhead drawing, if needed
				if (_this.draw_head()) {
					var hx_tip, hx_back;
					if (_this.clockwise()) {
						hx_tip = x1;
						x1 -= Feature.head_length;
						hx_back = x1;
					} else {
						hx_tip = x0;
						x0 += Feature.head_length;
						hx_back = x0;
					}

					// Unlike the body, the head is traced with a line, and
					// then created entirely with the fill color
					var head = paper.path(svg.move(hx_tip, _this.y) +
					                 svg.line(hx_back, _this.y - Feature.head_width/2.0) +
					                 svg.line(hx_back, _this.y + Feature.head_width/2.0) +
					                 svg.close());
					head.attr({"stroke-width": 0,
							   "fill":         _this.color()});
					_arrow_set.push(head);
				}

				// Body drawing
				if (x0 < x1 && _this.type() != ft.enzyme) { 
					// Compensating for the head may have "taken up" all
					// the room on the plasmid, in which case no arc needs
					// to be drawn

					// The body has no fill-color: it's just a thick line
					var body = paper.path(svg.move(x0, _this.y) +
						  				  svg.line(x1, _this.y));
					body.attr({"stroke-width": _this.width()});

					_arrow_set.push(body);
				} else if (_this.type() == ft.enzyme) { 
					// Restriction enzymes get drawn on their own
					var x_m = (x0 + x1)/2;

					var body = paper.path(svg.move(x_m, _this.y - _this.width()/2.0) +
					                      svg.line(x_m, _this.y + _this.width()/2.0));
					body.attr({"stroke-width": enzyme_weight});
					body.toBack();

					_arrow_set.push(body);
				}

				//_arrow_set.click(_click);
				//_arrow_set.hover(_mouse_over, _mouse_up);

				_feature_set.push(_arrow_set);

				// Apply the feature-wide properties to the whole feature
				_feature_set.attr({"stroke":         _this.color(),
								   "stroke-linecap": "butt",
								   "opacity":        _opacity,
								   "title":          _this.name()});

			} // END LinearFeature::draw()
			return _this;
		}; // END LinearFeature Class

		function draw_plasmid() {
			function draw_tic_mark(p) {
				var tic_mark_y = _this.plasmid_pos + _this.spacing;
				var tic_label_y = tic_mark_y + 1.5*_this.tic_mark_length;

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
			var title = sequence.length + ' bp';
			if (plasmid_name != "") {
				title = plasmid_name + "\n\n" + title;
			}
			var plasmid_label = paper.text(cx, title_y, title);
			plasmid_label.attr({"fill":      colors.plasmid,
								"font-size": plasmid_font_size, });

			// Set the scale to be the order of magnitude of seq_length
			// i.e. 100, 1000, 10, etc.
			var scale = Math.pow(10, Math.floor(Math.log(sequence.length )/Math.log(10)))
			for (var xx = scale; xx <=sequence.length ; xx += scale) {
				draw_tic_mark(xx);
			}

		}


		function draw() { // Draw the linear map
            // Extend basic features to get list of linear features
			extend_features();
            // Hide the right cutters
            //show_hide_cutters();
            // Resolve conflicts on the line, push some overlapping
            // features to other radii
			//var max_height = resolve_conflicts();
			//var label_height = max_height + label_height_offset; 
        
			paper = ScaleRaphael(_this.dom_id, _this.width, _this.height); // global
            for (var fx in features) {
                features[fx].initialize();
            }

			draw_plasmid();
			draw_features(); // Draw all the features initially
			//draw_labels(label_height); // Draw only the necessary labels

			// Rescale
			if (_this.final_width != _this.width ||
				_this.final_height != _this.height) {
				// "center" parameter just adds unnecessary CSS to the container
				// object to give it an absolute position: not what we need
				paper.changeSize(_this.final_width,_this.final_height,false,false)
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


	///////////////////////////////////////////////////////////////////
	// JSON Parsing and initialization
	window.GiraffeDraw = function (features_json) {
		sequence.features = []; // package scope
		sequence.length = features_json[0]; // package scope

		for (var ix = 1; ix < features_json.length; ix++) {
			sequence.features.push(new Feature(features_json[ix]));
		}

		// Now that the features are parsed, calculate how many instances
		// of each cutter type there are.
		cut_counts(sequence.features); 
	}

	// Private function
	// Calculate cut counts of all restriction enzyme
	function cut_counts(features) {
		var cut_counts = {}; 
		var f;
		// Calculate the counts
		for (var fx in features) {
			f = features[fx];
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
		for (var fx in features) {
			f = features[fx];
			if (f.type() == ft.enzyme)
				f.set_other_cutters(cut_counts[f.name()]);
		}
	}

})();
