/******************************************************************************
 *
 * GiraffeDraw: Javascript API for drawing plasmid maps from
 *              JSON plasmid feature data
 *
 * Copyright 2011 Addgene, Inc
 *
 * @version 0.1
 * @author  Mikhail Wolfson (wolfsonm@addgene.org)
 * @author  Benjie Chen     (benjie@addgene.org)
 *
 */

/////////////////////////////////////////////////////////////////////////////
//
// Drawing API: caller should call the GiraffeDraw function to get a drawing
// and parsing object. The caller should use the object's read() method to
// read in a JSON list of features, then call one of the object's <Foo>Map
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
//    gd.CircularMap({ "map_dom_id" : "some_id", ... })
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
//  opacity: opacity when features/enzymes is shown. if not 1.0, then
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
//  feature_click_callback: a callback that gets called when a feature
//  is clicked on. Argument to this callback is a Feature object.
//
//  digest: draw the map like a restriction digest: features are very
//  transparent and have no labels, and enzymes are drawn like normal.
//
//  digest_fade_factor: fade the features by this much (multiplicative, in [0,1.0])
//  Defaults to 0.5.
//
//  draw_tic_mark: true or false.
//
//  region_start_offset: for bp tick marks, start at this offset
//
//
// FRAMEWORK
// GiraffeDraw()
// |
// |- basic initialization of list of Features
// |
// |- read(JSON)
// |  feature parsing into "local" features object
// |
// |- CircularMap()
// |  circular feature drawing, which creates clones of the original Feature
// |  objects with extended properties
// |
// `- LinearMap()
//    linear feature drawing, which makes its own extended clones of objects

// Protect scope, but ensure that GiraffeDraw() is global
(function () {

    if (typeof(Object.create) !== 'function') {

        /**
         * Crockford-stype prototypal inheritence.
         *
         * @constructor Create a new object that inherits all public properites
         *              of the given object.
         *
         * @param   {Object} o the object to extend
         * @returns {o}        a new object whose __proto__ property points to o
         * @ignore
         */
        Object.create = function (o) {
            // Empty function object, to be used as a placeholder
            function F() {}

            // Set the prototype property of this function to point to o
            F.prototype = o;

            // Now, all objects created with F as a constructor will have their
            // internal __proto__ prototype pointer pointing to o, which means
            // that they will, effectively, use o as a "parent" for any properies
            // that they do not explicitly define.

            // return just such an object
            return new F();
        };
    }

    // Function.prototype.bind polyfill
    if (typeof(Function.bind) !== 'function') {

      Function.prototype.bind = function (obj) {
        var slice = [].slice,
            args = slice.call(arguments, 1),
            self = this,
            nop = function () {},
            bound = function () {
              return self.apply( this instanceof nop ? this : ( obj || {} ),
                                  args.concat( slice.call(arguments) ) );
            };

        nop.prototype = self.prototype;

        bound.prototype = new nop();

        return bound;
      };
    }

    // Implement Array.indexOf
    // courtesy: http://soledadpenades.com/2007/05/17/arrayindexof-in-internet-explorer/
    if (typeof(Array.indexOf) !== 'function') {
        Array.prototype.indexOf = function(obj) {
            for(var i= 0; i < this.length; i++){
                if(this[i] === obj){
                    return i;
                }
            }

            return -1;
        };
    }


    ///////////////////////////////////////////////////////////////////
    // Package-scope variables
    var _debug,
        ft,
        svg,
        colors;

    _debug = false;

    // Feature types: numerical codes to names
    ft = {
        feature:    1, promoter:   2, primer:        3,
        enzyme:     4, gene:       5, origin:        6,
        regulatory: 7, terminator: 8, exact_feature: 9,
        orf:       10
    };

    // For dealing with the SVG syntax
    svg = {
        move: function (x, y) {
            return this.to_path(['M', x, y]);
        },
        arc: function (r, x, y, l) {
            if (arguments.length < 4) {
                l = 0;
            }
            return this.to_path(['A', r, r, 0, l, 0, x, y]);
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
    // TODO: Make this map-specific and an option for this
    colors = {
        bg_text: "#aaa",
        plasmid: "#000",
        feature: "#f00",
        primer:  "#090",
        origin:  "#333",
        enzyme:  "#00c",
        orf:     "#00c8c8"
    };

/**
 * @class JSON parsing and feature drawing capabilities for a single feature.
 *        Only one GiraffeDraw per sequence.
 */
window.GiraffeDraw = function () {

    ///////////////////////////////////////////////////////////////////
    // GD-scope variables
    var seq_length,
        full_sequence = '',
        // Helpers to keep track of frequently used features
        all_features = [],
        orf_features = [],
        enzyme_features = [],
        std_features = [];


    // XXX hack to provide link to GiraffeDraw to the map functions.
    var gd = this;

    ///////////////////////////////////////////////////////////////////
    // Public data members

    // Expose the feature types
    this.Feature_Type = ft;


    /**
     * Parse JSON list of features
     *
     */
    this.read = function (json) {
        all_features = []; // package scope
        seq_length = json[0]; // package scope
        features_json = json[1];
        if (json.length > 2) { this.sequence = full_sequence = json[2]; }
        for (var ix = 0; ix < features_json.length; ix++) {
            var f = new Feature(features_json[ix]);
            all_features.push(f);
            if (f.is_enzyme()) {
                enzyme_features.push(f);
            } else if (f.is_orf()) {
                orf_features.push(f);
            } else {
                std_features.push(f);
            }
        }

        // These stuff are only available if we read the feature list
        this.all_features = all_features;
        this.enzyme_features = enzyme_features;
        this.orf_features = orf_features;
        this.std_features = std_features;

        // Now that the features are parsed, calculate how many instances
        // of each cutter type there are.
        cut_counts();
    };

    // Private function
    // Calculate cut counts of all restriction enzyme
    function cut_counts() {
        var cc = {},
            f, fx;

        // Calculate the counts
        for (fx = 0; fx < all_features.length; fx++) {
            f = all_features[fx];
            if (f.type() == ft.enzyme) {
                // Store indices, not Feature objects, because
                // the objects will change, depending on
                // the kind of map is drawn
                if (f.name() in cc) {
                    cc[f.name()].push(fx);
                } else {
                    cc[f.name()] = [fx];
                }
            }
        }
        // Store them for each enzyme feature
        for (fx = 0; fx < all_features.length; fx++) {
            f = all_features[fx];
            if (f.type() == ft.enzyme) {
                f.set_other_cutters(cc[f.name()]);
            }
        }
    }

    ///////////////////////////////////////////////////////////////////
    // Basic Feature class (private, internal class)
    function Feature(feat) {

        // Private data members, from the properties of the feature argument
        // since hash tables and objects are one and the same in js
        var _name = feat.feature;
        var _start = parseInt(feat.start, 10);
        var _end = parseInt(feat.end, 10);
        var _type = parseInt(feat.type_id, 10);
        var _default_show_feature = 1;
        // Not all features have this option, e.g. annotated ones we
        // have to show the feature.
        if ('show_feature' in feat) {
            _default_show_feature = parseInt(feat.show_feature, 10);
        }
        var _clockwise = feat.clockwise;
        var _cut = parseInt(feat.cut, 10); // only for enzymes;
        var _other_cutters = []; // only for enzymes;
        var _id = Feature.nfeat;

        // Accessors for private properties set at creation
        this.name = function() { return _name; };
        this.start = function() { return _start; };
        this.end = function() { return _end; };
        this.type = function() { return _type; };
        this.clockwise = function() { return _clockwise; };
        this.default_show_feature = function() { return _default_show_feature; };
        // returns - 1 if not enzyme
        this.cut = function() { return _type == ft.enzyme ? (_cut ? _cut : _start) : -1; };
        this.actually_have_cut = function() {
            return _type == ft.enzyme ? (_cut ? true : false) : false;
        };

        this.is_enzyme = function() { return _type == ft.enzyme; };
        this.is_orf = function() { return _type == ft.orf; };

        this.crosses_boundary = function () { return _end < _start; };
        this.id = function () { return _id; };

        // Enzyme-only data access methods
        // gives 0 if not enzyme
        this.cut_count = function() { return _other_cutters.length; };
        // gives [] if not enzyme
        this.other_cutters = function() { return _other_cutters; };
        // Mutator for other_cutters: only for enzymes
        this.set_other_cutters = function(c) {
            if (_type == ft.enzyme) {
                _other_cutters = c;
            }
        };

        this.clockwise_sequence = function() {
            if (full_sequence.length < this.start() ||
                full_sequence.length < this.end()) { return ''; }
            var s = '';
            // positions we read from JSON are 1-based
            if (this.end() >= this.start()) {
                s = full_sequence.substring(this.start()-1,this.end()-1+1);
            }
            else {
                s = full_sequence.substring(this.start()-1)+
                    full_sequence.substring(0,this.end()-1+1);
            }
            return s;
        };

        Feature.nfeat++;
        return this;
    } // END Basic Feature Class

    Feature.nfeat = 0;


    ///////////////////////////////////////////////////////////////////
    // Generic Features that are drawn on a map
    function DrawnFeature(basic_feature, map) {
        var thi$ = Object.create(basic_feature);

        thi$.label_drawn = false;
        thi$.feature_set = undefined;
        thi$.arrow_set = undefined;
        thi$.label_set = undefined;

        // Visual properties
        thi$.visible = true;
        thi$.labeled = true;

        // Easy access to drawing objects
        thi$.map = map;

        // Defaults
        thi$.color = colors.feature;
        thi$.width = map.feature_width();
        thi$.draw_head = false;
        thi$.opacity = map.feature_opacity();
        thi$.opaque = false; // holds opacity for clicks
        set_properties();

        // Type-based property selection
        function set_properties() {

            switch(thi$.type()) {
                case ft.promoter:
                case ft.primer:
                    thi$.draw_head = true;
                    thi$.color = colors.primer;
                    break;
                case ft.terminator:
                    thi$.color = colors.primer;
                    break;
                case ft.regulatory:
                case ft.origin:
                    thi$.color = colors.origin;
                    break;
                case ft.enzyme:
                    thi$.color = colors.enzyme;
                    thi$.width = map.enzyme_width();
                    thi$.opacity = map.enzyme_opacity();
                    break;
                case ft.orf:
                    thi$.color = colors.orf;
                    thi$.draw_head = true;
                    break;
                case ft.gene:
                    thi$.draw_head = true;
                    break;
                default:
                    break;
            }
        }

        thi$.initialize = function() {
            this.feature_set = this.map.paper.set();
            this.arrow_set = this.map.paper.set();
            this.label_set = this.map.paper.set();
        };

        // Actions for interactivity
        // Generic fading animation/property setting mechanism
        thi$.fade = function (props, line_props) {
            var sets = this.map.paper.set(),
                lines = this.map.paper.set(),
                fx, f;

            sets.push(this.feature_set);
            lines.push(this.label_set[0]); // label line

            // Cutters: highlight other examples of this enzyme if
            // it's a multi-cutter
            if (this.type() == ft.enzyme) {
                for (fx = 0; fx < this.other_cutters().length; fx++) {
                    f = this.map.features[this.other_cutters()[fx]];
                    sets.push(f.feature_set);
                    lines.push(f.label_set[0]);
                }
            }

            if (this.map.fade_time) {
                sets.animate(props, this.map.fade_time);
                lines.animateWith(sets, line_props, this.map.fade_time);
            } else {
                sets.attr(props);
                lines.attr(line_props);
            }
        };

        thi$.bolder = function () {
            var props = {"opacity": this.map.bold_opacity(),
                         "font-weight": "bold" };
            if (this.type() == ft.enzyme) {
                props["stroke-width"] = this.map.enzyme_bold_weight;
            }
            var line_props = {"stroke": colors.plasmid,
                              "stroke-width": this.map.label_line_bold_weight()};
            this.fade(props, line_props);
        };

        thi$.lighter = function () {
            var props = {"opacity": this.opacity,
                         "font-weight":"normal"};
            if (this.type() == ft.enzyme) {
                props["stroke-width"] = this.map.enzyme_weight();
            }
            var line_props = {"stroke": colors.bg_text,
                              "stroke-width": this.map.label_line_weight()};
            this.fade(props, line_props);
        };

        // Toggle solid/light upon click
        thi$.click = function(event) {
            if (this.map.feature_click_callback) {
                this.map.feature_click_callback(basic_feature);
            } 
            /* benjie: don't do this for now...
            else {
               
                if (this.opaque) {
                    this.lighter();
                    this.opaque = false;
                } else {
                    this.bolder();
                    this.opaque = true;
                }
                
            }*/
        };

        // Hovering: solid/light upon mouseover
        thi$.mouse_over = function (event) {
            // this object context will change when called, so we use thi$
            if (!this.opaque) {
                this.bolder();
            }
        };

        thi$.mouse_up = function (event) {
            // this object context will change when called, so we use thi$
            if (!this.opaque) {
                this.lighter();
            }
        };

        // Showing and hiding

        thi$.hide = function () {
            if (this.visible) {
                if (this.feature_set) { this.feature_set.hide(); }
                this.visible = false;
                this.labeled = false;
            }
        }; // END DrawnFeature::hide()

        thi$.show = function () {
            if (!this.visible) {
                if (this.feature_set) { this.feature_set.show(); }
                if (!this.labeled) {
                    if (this.label_set) { this.label_set.hide(); }
                }
                this.visible = true;
            }
        }; // END DrawnFeature::show()

        thi$.hide_label = function () {
            if (this.labeled) {
                if (this.label_set) { this.label_set.hide(); }
                this.labeled = false;
            }
        }; // END DrawnFeature::hide_label()

        thi$.show_label = function () {
            if (!this.labeled) {
                if (this.label_set) { this.label_set.show(); }
                this.labeled = true;
            }
        }; // END DrawnFeature::show_label()

        thi$.clear_label = function () {
            if (this.label_drawn) {
                if (this.label_set) {
                    this.label_set.unclick(this.click);
                    this.label_set.unhover(this.mouse_over, this.mouse_up);
                    this.label_set.remove();
                    this.label_set = this.map.paper.set();
                }
                this.labeled = false;
                this.label_drawn = false;
            }
        }; // END DrawnFeature::clear_label()

        // Common properties

        /** Automatically provide the cut site in cutter labels */
        thi$.label_name = function () {
            var label_name = this.name();
            if (this.type() == ft.enzyme && this.actually_have_cut()) {
                label_name += " (" + this.cut() + ")";
            }
            return label_name;
        }; // END DrawnFeature::label_name()

        return thi$;
    } // END DrawnFeature()



    ///////////////////////////////////////////////////////////////////
    // Generic Map prototype class
    function Map(options) {

        var _cutters_to_show,
            _is_digest,
            _digest_fade_factor,
            _feature_opacity, _enzyme_opacity, _bold_opacity,
            _feature_width, _enzyme_width,
            _enzyme_weight, enzyme_bold_weight,
            _label_line_weight, _label_line_bold_weight, _label_font_size,
            _plasmid_font_size, _plasmid_name, _draw_tic_mark, _draw_plasmid_size,
            _region_start_offset, thi$ = {};

        init();

        function init() {
            // Digest flag: read-only, set once
            _is_digest = false;
            if ('digest' in options) {
                _is_digest = options.digest;
            }

            _digest_fade_factor = 0.15;
            if ('digest_fade_factor' in options) {
                _digest_fade_factor = options.digest_fade_factor;
            }

            // Feature visual properties
            _feature_width = 15;
            _enzyme_width = 25;
            _enzyme_weight = 2; // Restriction enzymes are drawn differently
                                // This controls their "thickness" on the map
            _enzyme_bold_weight = 3; // How thick when highlighted

            // Labels and other text
            _label_line_weight = 1;
            _label_line_bold_weight = 1.5 * _label_line_weight;
            _label_font_size = '13pt';

            // Plasmid text
            _plasmid_font_size = '16pt';
            _plasmid_name = '';
            if ('plasmid_name' in options) {
                _plasmid_name = options.plasmid_name;
            }
            _region_start_offset = 0;
            if ('region_start_offset' in options) { _region_start_offset = parseInt(options['region_start_offset']); }

            // Draw plasmid size
            _draw_plasmid_size = true;
            if ('draw_plasmid_size' in options &&
                (options.draw_plasmid_size === 0 ||
                 options.draw_plasmid_size === '0' ||
                 options.draw_plasmid_size === false ||
                 options.draw_plasmid_size === 'false')) {
                _draw_plasmid_size = false;
            }

            // Draw tic mark
            _draw_tic_mark = true;
            if ('draw_tic_mark' in options &&
                (options.draw_tic_mark === 0 ||
                 options.draw_tic_mark === '0' ||
                 options.draw_tic_mark === false ||
                 options.draw_tic_mark === 'false')) {
                _draw_tic_mark = false;
            }

            // Opacities
            _feature_opacity = 0.7;
            _enzyme_opacity = 0.7;
            if ('opacity' in options) {
                _feature_opacity = parseFloat(options.opacity);
                _enzyme_opacity  = parseFloat(options.opacity);
            }
            _bold_opacity = 1.0;

            if (_is_digest) {
                _feature_opacity = _feature_opacity * _digest_fade_factor;
            }

            // Cutters to show
            _cutters_to_show = [1];
            if ('cutters' in options) {
                _cutters_to_show = options.cutters;
            }

            // Where to draw the map
            thi$.map_dom_id = 'giraffe-draw-map';
            if ('map_dom_id' in options) {
                thi$.map_dom_id = options.map_dom_id;
            }

            thi$.paper  = undefined; // To be used for RaphaelJS;
            thi$.label_offset = 0;

            thi$.features = [];
            thi$.show_all_features = false;

            thi$.width = 800;
            thi$.height = 800;

            // Final size: to be scaled down to this at the end
            thi$.final_width = 640;
            thi$.final_height = 640;
            if ('map_width' in options) {
                thi$.final_width = parseInt(options.map_width, 10);
            }
            if ('map_height' in options) {
                thi$.final_height = parseInt(options.map_height, 10);
            }

            // Animation properties
            thi$.fade_time = 0;
            if ('fade_time' in options) {
                thi$.fade_time = parseInt(options.fade_time, 10);
            }

            // Callback
            thi$.feature_click_callback = undefined;
            if ('feature_click_callback' in options) {
                thi$.feature_click_callback = options.feature_click_callback;
            }
        }

        thi$.is_digest = function () { return _is_digest; };
        thi$.feature_width = function () { return _feature_width; };
        thi$.enzyme_width = function () { return _enzyme_width; };
        thi$.enzyme_weight = function () { return _enzyme_weight; };
        thi$.enzyme_bold_weight = function () { return _enzyme_bold_weight; };

        thi$.feature_opacity = function () { return _feature_opacity; };
        thi$.enzyme_opacity = function () { return _enzyme_opacity; };
        thi$.bold_opacity = function () { return _bold_opacity; };

        thi$.label_line_weight = function() { return  _label_line_weight; };
        thi$.label_line_bold_weight = function() { return  _label_line_bold_weight; };
        thi$.label_font_size = function() { return  _label_font_size; };
        thi$.plasmid_font_size = function() { return  _plasmid_font_size; };
        thi$.plasmid_name = function() { return  _plasmid_name; };
        thi$.region_start_offset = function() { return  _region_start_offset; };
        thi$.draw_tic_mark = function() { return  _draw_tic_mark; };
        thi$.draw_plasmid_size = function() { return  _draw_plasmid_size; };

        thi$.extend_features = function () {
            var bfx,
                df, ef;

            for (bfx = 0; bfx < all_features.length; bfx++) {
                df = new DrawnFeature(all_features[bfx], this);
                ef = new this.FeatureType(df);
                this.features.push(ef);
            }
        };

        thi$.initialize_features = function() {
            var fx;

            for (fx = 0; fx < this.features.length; fx++) {
                this.features[fx].initialize();
            }
        };

        // Make sure that the appropriate cutters are shown
        thi$.show_hide_cutters = function () {
            var fx, f;

            for (fx = 0; fx< this.features.length; fx++) {
                f = this.features[fx];
                if (f.default_show_feature() || this.show_all_features) {
                    // Only draw enzymes if they are in the list of
                    // cutters to show - i.e. 1 cutter, 2 cutters,
                    // etc.
                    // If the list is undefined, draw all cutters
                    if (f.type() == ft.enzyme) {
                        if (typeof(_cutters_to_show) !== 'undefined' &&
                            _cutters_to_show.indexOf(f.cut_count()) < 0) {
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
        };

        thi$.redraw_cutters = function (new_cutters_to_show) {
            if (arguments.length > 0 &&
                typeof(new_cutters_to_show) != 'undefined') {
                _cutters_to_show = new_cutters_to_show;
            } else {
                // Undefined cutter list: hide all cutters
                _cutters_to_show = [];
            }
            this.redraw(false);
        };

        thi$.show_extra_features = function () {
            this.show_all_features = true;
            this.redraw(true);
        };

        thi$.hide_extra_features = function () {
            this.show_all_features = false;
            this.redraw(true);
        };

        // any number of arguments will work. accessed via arguments object
        function apply_to_feature_type(recalc) {
            var funcs_to_call = arguments;
            return (function (type_name) {
                var fx, ax,
                    type_code = ft[type_name];

                for (fx = 0; fx < this.features.length; fx++) {
                    if (this.features[fx].type() == type_code) {
                        // All arguments after the first are functions to call
                        for (ax = 1; ax < funcs_to_call.length; ax++) {
                            this.features[fx][funcs_to_call[ax]]();
                        }
                    }
                }

                this.redraw(recalc);
            });
        }

        thi$.hide_feature_type = apply_to_feature_type(true, "hide");
        thi$.show_feature_type = apply_to_feature_type(true, "show", "show_label");
        thi$.show_feature_label_type = apply_to_feature_type(false, "show_label");
        thi$.hide_feature_label_type = apply_to_feature_type(false, "hide_label");

        // any number of arguments will work. accessed via arguments object
        function apply_to_feature() {
            var funcs_to_call = arguments;

            return (function (id) {
                var ax;


                // All arguments after the first are functions to call
                for (ax = 0; ax < funcs_to_call.length; ax++) {
                    this.features[id][funcs_to_call[ax]]();
                }
            });
        }

        thi$.hide_feature = apply_to_feature("hide_label", "hide");
        thi$.show_feature = apply_to_feature("show", "show_label");
        thi$.show_feature_label = apply_to_feature("show_label");
        thi$.hide_feature_label = apply_to_feature("hide_label");


        thi$.draw_features = function () {
            var fx;

            for (fx = 0; fx < this.features.length; fx++) {
                this.features[fx].draw();
            }
        };

        thi$.draw = function() { // Draw the map
            // Extend basic features to get list of circular features:
            // do this only once, the first time the map is created.
            this.extend_features();

            this.update(true);
        };

        thi$.clear = function () {
            var map_dom,
                kids,
                kx;

            // Since we don't want to depend on jQuery, use the raw DOM (ugh!)
            map_dom = document.getElementById(this.map_dom_id);
            kids = map_dom.childNodes;

            for (kx = 0; kx < kids.length; kx++) {
                map_dom.removeChild(kids[kx]);
            }
        };

        thi$.redraw = function (recalc) {
            this.clear();
            this.update(recalc);
        };

        thi$.recalculate_positions = function() { // Redraw the map
            // Resolve conflicts on the circle, push some overlapping
            // features to other radii
            this.max_extent = this.resolve_conflicts();
            this.label_pos = this.max_extent + this.label_offset;
        };

        thi$.update = function(recalc) { // Redraw the map
            var fx, f;

            if (arguments.length < 1) {
                recalc = true;
            }

            if (recalc) {
                this.recalculate_positions();
            }

            // For digest maps, make features more transparent, and hide their labels
            if (this.is_digest()) {
                for (fx = 0; fx < this.features.length; fx++) {
                    f = this.features[fx];
                    if (f.type() != ft.enzyme) {
                        f.hide_label();
                    }
                }
            }

            // Hide the right cutters
            this.show_hide_cutters();

            // Determine which labels are in which lists
            this.set_label_lists();

            this.set_bounding_box();

            this.paper = ScaleRaphael(this.map_dom_id, this.width, this.height); // global
            this.initialize_features();

            this.draw_plasmid();
            this.draw_features(); // Draw all the features initially
            this.draw_labels(); // Draw only the necessary labels

            this.rescale();
        };

        // Centralized mechanism for exposing public properties of maps to the ouside user
        // Do not expose any of the mechanics, just an API of utility methods useful
        // in turning on and off features, and the gd object itself, for static things
        // like feature types, etc.
        thi$.expose = function() {

            // Export the main properties as part of a Map-like object
            // descendant classes inheret this ability, and because of the
            // this pointer, will return /Their own/ draw and features objects
            return {
                redraw_cutters: this.redraw_cutters.bind(this),
                show_feature_type: this.show_feature_type.bind(this),
                hide_feature_type: this.hide_feature_type.bind(this),
                show_feature_label_type: this.show_feature_label_type.bind(this),
                hide_feature_label_type: this.hide_feature_label_type.bind(this),
                show_feature: this.show_feature.bind(this),
                hide_feature: this.hide_feature.bind(this),
                show_feature_label: this.show_feature_label.bind(this),
                hide_feature_label: this.hide_feature_label.bind(this),
                show_extra_features: this.show_extra_features.bind(this),
                hide_extra_features: this.hide_extra_features.bind(this),
                is_digest: this.is_digest,
                gd: gd
            };
        };

        return thi$;
    }

    // Utility function for creating wrapper closures around Map classes
    Map.make = function (MapType) {

        return function(options) {
            // Create a map object
            var thi$ = MapType(options);

            // Draw it
            thi$.draw();

            // Expose only the parts of it we care about
            return thi$.expose();
        };
    };

    this.CircularMap = Map.make(CircMap);
    this.LinearMap = Map.make(LinMap);

    ///////////////////////////////////////////////////////////////////
    // Circular Map Drawing Class
    function CircMap(options) {

        // Inherit the common Map functions
        var thi$ = Object.create(new Map(options));
        thi$.FeatureType = CircularFeature;
        thi$.x_shift_on_labels = 60;

        // Paper setup - not the final width, but how we will draw the
        // map, we will scale later on
        var cx = thi$.width/2;
        var cy = thi$.height/2;

        // Global plasmid info
        var plasmid_start = 90; // degrees

        // Loop radii
        var radius_spacing = 20; // spacing
        var plasmid_radius = 200;
        var inner_radius = plasmid_radius - radius_spacing;
        var outer_radius = plasmid_radius + radius_spacing;
        thi$.label_offset = 10;
        if ('label_offset' in options) {
            thi$.label_offset = parseInt(options.label_offset, 10);
        }

        var head_width = 25;
        var head_length = 7;

        // Overlaps
        var min_overlap_cutoff = -0.1;// in degrees
        var min_overlap_pct = 0.01;
        var min_overlap_feature_size = 0.5; // in degrees

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

        thi$.label_list_section_angle = function (section) {
            return plasmid_start-label_section_degree/2.0-section*label_section_degree;
        };

        ///////////////////////////////////////////////////////////////////
        // Circular Feature class
        function CircularFeature(basic_feature) {

            // Create a prototypal descendent of the basic_feature to expand
            var thi$ = Object.create(basic_feature);

            // The result of this function will be a CircularFeature object


            // Radius is public, unlike other properties, which are permanent
            thi$.radius = plasmid_radius; // Default to plasmid radius, can be changed
                                          // later on by other methods

            // Calculated properties

            function normalize(deg) {
                if (deg < plasmid_start - 360) {
                    deg += 360;
                } else if (deg > plasmid_start) {
                    deg -= 360;
                }

                return deg;
            }

            // Degree conversion, for overlap calculation:
            // for these functions, the sequence starts at 90 degrees and goes down.
            thi$.real_start = function() {
                var sd;
                // Take the minimum head size into account. Only need to do this
                // when the head is drawn and pointing clockwise, to
                // "push the start back."
                if (this.draw_head && thi$.clockwise()) {
                    sd = convert.pos_to_angle(thi$.end()) + thi$.real_size();
                } else { // Headless feature, or head is pointing the wrong way.
                         // Just give its typical start position
                    sd = convert.pos_to_angle(thi$.start());
                }

                return normalize(sd);
            };

            thi$.real_center = function() {
                var cd;

                // Remember: minus means going along the circle
                cd = thi$.real_start() - thi$.real_size()/2.0;

                return normalize(cd);
            };

            thi$.real_end = function() {
                var ed;
                // Take the minimum head size into account. Only need to do this
                // when the head is drawn and pointing counterclockwise, to
                // "push the end forward."
                if (this.draw_head && !thi$.clockwise()) { // Take the minimum head size into account
                    ed = convert.pos_to_angle(thi$.start()) - thi$.real_size();
                } else { // Headless feature, or head is pointing the wrong way.
                         // Just give its typical end position
                    ed = convert.pos_to_angle(thi$.end());
                }

                return normalize(ed);
            };

            thi$.bp_size = function() {
                if (thi$.crosses_boundary()) {
                    return thi$.end()+seq_length-thi$.start()+1;
                }
                else {
                    return thi$.end()-thi$.start()+1;
                }
            }

            thi$.real_size = function() {
                var szd; // size in degrees
                // Normal definition of size
                if (thi$.crosses_boundary()) {
                    // Start and end are flipped here: non-intuitive
                    szd = convert.seq_length_to_angle(seq_length - thi$.start() + thi$.end() + 1);
                } else {
                    szd = convert.seq_length_to_angle(thi$.end() - thi$.start() + 1);
                }

                // Head size: return this if it's bigger
                if (this.draw_head) {
                    // Convert the head length into degrees, just as you do
                    // in the draw() method. Must recalcualte every time, as
                    // radius may have changed
                    var r_p = Math.sqrt(thi$.radius*thi$.radius +
                            head_length*head_length);
                    var hszd = Raphael.deg(Math.asin(head_length/r_p));
                    if (hszd > szd) {
                        szd = hszd;
                    }
                }

                return szd;
            };


            // Feature drawing
            thi$.draw = function () {
                if (!this.visible) { return; }

                    // Convert from sequence positions to angles
                var a0 = convert.pos_to_angle(this.start()),
                    a1 = convert.pos_to_angle(this.end()),
                    // Head dimensions
                    r_p, a_b, a_p, xy_p, xy_b, xy_t, head,
                    // Arc dimensions
                    xy0, xy1, arc, large_angle, a_cut;

                // Create the draw feature, a set which will have the head
                // and arc pushed onto it as necessary.

                // Arrowhead drawing, if needed
                if (this.draw_head) {

                    // Arrow tip point lines up with a0 or a1 and it points
                    // tangent to the circle.
                    // We need to figure out how many radians the arrow takes up
                    // in order to adjust a0 or a1 by that amount, and to set the
                    // base of the triangle even with that angle
                    r_p = Math.sqrt(this.radius*this.radius +
                            head_length*head_length);
                    // "height" of the arrowhead, in degrees
                    a_p = Raphael.deg(Math.asin(head_length/r_p));
                    // Adjust the appropriate edge to compensate for the arrowhead
                    if (this.clockwise()) {
                        a_b = (a1 + a_p) % 360 ; // base angle
                        a_p = a1;       // point angle
                        a1  = a_b;      // adjust arc edge
                    } else {
                        a_b = (a0 - a_p) % 360 ; // base angle
                        a_p = a0;       // point angle
                        a0  = a_b;      // adjust arc edge
                    }
                    xy_p = convert.polar_to_rect(this.radius, a_p);

                    // bottom and top points, rectangular
                    xy_b = convert.polar_to_rect(this.radius - head_width/2.0, a_b);
                    xy_t = convert.polar_to_rect(this.radius + head_width/2.0, a_b);

                    // Unlike the arc, the head is traced with a line, and
                    // then created entirely with the fill color
                    head = this.map.paper.path(svg.move(xy_p.x, xy_p.y) +
                                          svg.line(xy_b.x, xy_b.y) +
                                          svg.line(xy_t.x, xy_t.y) +
                                          svg.close());
                    head.attr({"stroke-width": 0,
                               "fill":         this.color});
                    this.arrow_set.push(head);
                }

                // Arc drawing
                if ((this.crosses_boundary() || a1 < a0) && this.type() != ft.enzyme) {
                    // Compensating for the head may have "taken up" all
                    // the room on the plasmid, in which case no arc needs
                    // to be drawn

                    // Rectangular coordinates of the edges of the arc:
                    // arcs are drawn counterclockwise, even though the plasmid
                    // sequence increases clockwise, so we flip the
                    // indices
                    xy0 = convert.polar_to_rect(this.radius, a1);
                    xy1 = convert.polar_to_rect(this.radius, a0);

                    // The arc has no fill-color: it's just a thick line
                    large_angle = this.real_size() > 180 ? 1 : 0;

                    arc = this.map.paper.path(svg.move(xy0.x, xy0.y) +
                                         svg.arc(this.radius, xy1.x, xy1.y,
                                                 large_angle));
                    arc.attr({"stroke-width": this.width});

                    this.arrow_set.push(arc);
                } else if (this.type() == ft.enzyme) {
                    a_cut = convert.pos_to_angle(this.cut());

                    // Restriction enzymes get drawn on their own
                    xy0 = convert.polar_to_rect(this.radius - this.map.enzyme_width()/2.0, a_cut);
                    xy1 = convert.polar_to_rect(this.radius + this.map.enzyme_width()/2.0, a_cut);

                    // Not really an arc, just a line, but left this way
                    // for consistency
                    arc = this.map.paper.path(svg.move(xy0.x, xy0.y) +
                                         svg.line(xy1.x, xy1.y));
                    arc.attr({"stroke-width": this.map.enzyme_weight()});
                    arc.toBack();

                    this.arrow_set.push(arc);
                }

                this.arrow_set.click(this.click.bind(this));
                this.arrow_set.hover(this.mouse_over.bind(this), this.mouse_up.bind(this));

                this.feature_set.push(this.arrow_set);

                // Apply the feature-wide properties to the whole feature
                this.feature_set.attr({"stroke":         this.color,
                                   "stroke-linecap": "butt",
                                   "opacity":        this.opacity});

                if (this.map.is_digest() && this.type() != ft.enzyme) {
                    this.feature_set.attr("title", this.label_name());
                }

            }; // END CircularFeature::draw()

            // Should we draw the label?
            thi$.should_draw_label = function () {
                // Don't bother unless we need to
                if (!this.visible || !this.labeled) {
                    return false;
                }
                return true;
            }; // END CircularFeature::should_draw_label()

            // Set which label list this feature's label should be in
            thi$.set_label_list = function () {
                if (!thi$.should_draw_label()) { return; }

                // Figure out the center of the feature
                var a_c;
                if (this.type() == ft.enzyme) {
                    a_c = convert.pos_to_angle(this.cut());
                } else {
                    a_c = this.real_center();
                }

                var adjust_a_c = a_c;
                if (adjust_a_c < 0) { adjust_a_c += 360; }

                // Figure out which section this label is in: divide
                // the grid up into eight sections.
                var section = Math.floor((plasmid_start-a_c)/label_section_degree);
                section = section % label_f_c.length;

                var l = label_f_c[section].length;
                label_f_c[section][l] = [adjust_a_c,thi$.label_name()];

            }; // END CircularFeature::set_label_list()

            // Draw the label associated with that feature
            thi$.draw_label = function () {
                if (!thi$.should_draw_label()) { return; }

                if (this.label_drawn) {
                    thi$.clear_label();
                }

                // Figure out the center of the feature
                var a_c;
                if (this.type() == ft.enzyme) {
                    a_c = convert.pos_to_angle(this.cut());
                } else {
                    a_c = this.real_center();
                }

                var adjust_a_c = a_c;
                if (adjust_a_c < 0) { adjust_a_c += 360; }

                var xy0 = convert.polar_to_rect(thi$.radius, a_c);

                // Figure out which section this label is in: divide
                // the grid up into eight sections.
                var section = Math.floor((plasmid_start - a_c)/label_section_degree);
                section = section % label_f_c.length;

                // Figure out position in the label list.
                var pos_ls = 0;

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
                var xy1 = {};
                xy1.x = label_list_pos[section][0];
                xy1.y = label_list_pos[section][1];

                // We want to minimize the number of label lines that
                // cross. Which means depends on which section we are in,
                // we draw labels in different orders. See draw_labels on
                // how the positions of each label section are
                // computed. Remember, because we are sorted by
                // label_f_c, we are going counter clockwise on the
                // circle, drawing label for the feature with higher
                // bp first. label_f_c has the lower x and y
                // coordinates of each section.
                if (section === 0 || section == 1) {
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
                var label_line = this.map.paper.path(
                    svg.move(xy0.x, xy0.y) + svg.line(xy1.x, xy1.y)
                );
                label_line.attr({"stroke": colors.bg_text,
                                 "stroke-width": this.map.label_line_weight(),
                                 "opacity": 0.5 });

                // Enzymes show their cut sites in the label
                var label_name = thi$.label_name();
                var label = this.map.paper.text(xy1.x, xy1.y, label_name);

                if (section > 3) {
                    // Left half of wheel: align right
                    label.attr({"text-anchor": "end"});
                }
                else {
                    // Right half of wheel: align left
                    label.attr({"text-anchor": "start"});
                }

                label.attr({"fill": this.color,
                            "font-size": this.map.label_font_size(),
                            "opacity": 1.0 });

                this.label_set.push(label_line);
                this.label_set.push(label);

                // Handlers
                this.label_set.click(this.click.bind(this));
                this.label_set.hover(this.mouse_over.bind(this), this.mouse_up.bind(this));

                // Only push label_line, so when we fade in and out,
                // we don't also fade the label.
                this.feature_set.push(label_line);

                this.labeled = true;
                this.label_drawn = true;
            }; // END CircularFeature::draw_label()

            return thi$;
        } // END CircularFeature Class

        ///////////////////////////////////////////////////////////////////
        // Drawing utility functions

        // Circle setup
        thi$.draw_plasmid = function () {
            // Tic marks
            var tic_mark_length = 15;
            var tic_mark_radius = inner_radius - tic_mark_length/2;
            var tic_label_radius = tic_mark_radius - 1.5*tic_mark_length;

            function draw_tic_mark(a) {
                var r0, r1, xy0, xy1, tic,
                    xyl, label_pos, label;
 
                // use thi$ for paper because `this` is killed by local scope
                r0 = tic_mark_radius - tic_mark_length/2;
                r1 = tic_mark_radius + tic_mark_length/2;
                xy0 = convert.polar_to_rect(r0,a);
                xy1 = convert.polar_to_rect(r1,a);
                tic = thi$.paper.path(svg.move(xy0.x, xy0.y) +
                                     svg.line(xy1.x, xy1.y));
                tic.attr({"stroke": colors.bg_text});

                xyl = convert.polar_to_rect(tic_label_radius, a);

                label_pos = convert.angle_to_pos(a);
                if (label_pos == 1) {
                    label_pos = seq_length;
                }
                label = thi$.paper.text(xyl.x, xyl.y, String(label_pos));
                label.attr({"fill": colors.bg_text});

                if (a < plasmid_start || a > 360 - plasmid_start) { // Right half of wheel: align right
                    label.attr({"text-anchor": "end"});
                } else if (a > plasmid_start && a < 360 - plasmid_start) { // Left half of wheel: align left
                    label.attr({"text-anchor": "start"});
                } // Top and bottom default to middle, which is correct
            }

            var plasmid = this.paper.circle(cx, cy, plasmid_radius);
            plasmid.attr("stroke", colors.plasmid);
            var title = '';
            if (this.draw_plasmid_size()) { title = seq_length + ' bp'; }
            if (this.plasmid_name() !== "") {
                title = this.plasmid_name() + "\n\n" + title;
            }

            var plasmid_label = this.paper.text(cx, cy, title);
            plasmid_label.attr({"fill":      colors.plasmid,
                                "font-size": this.plasmid_font_size() });

            if (this.draw_tic_mark()) {
                for (var ang = 0; ang < 360; ang += 30) {
                    draw_tic_mark(ang);
                }
            }
        };

        // Move features that overlap to other radii.
        thi$.resolve_conflicts = function () {
            var conflicts = 0,
                rad = plasmid_radius, // current radius
                rx = 1,               // radius counter
                max_rad = plasmid_radius,
                fx, f,
                biggest_size, biggest_feature, furthest_point,
                new_rad, new_size,
                eff_furthest_point, overlap;

            // Utility functions
            function push(winner, loser) {
                var ppfx, pfx, pf,
                    can_push;

                // Record that the push happened
                winner.pushed_features.push(loser);
                conflicts++;

                // Do it
                loser.radius = new_rad;

                if (_debug) {
                    console.warn(loser.name() +
                        " (" + loser.real_start() + ", " + loser.real_end() + ")" +
                        " pushed by " +
                        winner.name() +
                        " (" + winner.real_start() + ", " + winner.real_end() + ")");
                }

                // Since loser was pushed, un-push all the
                // features that it pushed, as long as
                // those features are not in conflict with the winner,
                // or with their own, previously pushed features, which are
                // now unpushed
                for (pfx = 0; pfx < loser.pushed_features.length; pfx++) {
                    pf = loser.pushed_features[pfx];

                    // Check for conflict with the winner feature itself
                    // If there's no conflict, we can push it back safely.
                    if (pf.real_start() - winner.real_end() <= min_overlap_cutoff ||
                        winner.real_start() - pf.real_end() <= min_overlap_cutoff) {

                        // Check for conflict with previously pushed features
                        // that may have been unpushed
                        can_push = true;
                        for (ppfx = 0; ppfx < pf.pushed_features.length; ppfx++) {
                            if (pf.pushed_features[ppfx].radius == rad) {
                                can_push = false;
                                break;
                            }
                        }

                        // Finally!
                        if (can_push) {
                            if (_debug) {
                                console.warn(pf.name() + " unpushed, because " +
                                             loser.name() + " pushed by " + winner.name());
                            }
                            pf.radius = rad;
                        }
                    }
                }
            }

            // Main method body
            // reset radii in case this is being done any time but the first
            for (fx = 0; fx < this.features.length; fx++) {
                this.features[fx].radius = plasmid_radius;
            }

            // Don't try to resolve conflicts too many times
            var max_tries = 11;

            do {
                max_tries--;

                // Keep alternating between inside and outside the plasmid.
                delta = rx *radius_spacing;

                if (rx %2 === 0) { // Even rx
                    new_rad = rad + delta;
                } else {
                    new_rad = rad - delta;
                }

                conflicts = 0; // Assume you have no conflicts until you find some

                // Clear the record of who pushed whom
                for (fx = 0; fx < this.features.length; fx++) {
                    this.features[fx].pushed_features = [];
                }

                biggest_size = 0;
                biggest_feature = undefined;
                furthest_point = 0; // Start at a complete lower bound

                // Go through the feature list twice, to make sure that features
                // that cross the boundary are resolved
                for (fx = 0; fx < this.features.length; fx++) {
                    f = this.features[fx % this.features.length];

                    // Stop when we are no longer working with overlapped features
                    if (fx > this.features.length &&
                        biggest_feature.start() < biggest_feature.end()) { break; }

                    if (f.visible && f.type() != ft.enzyme && f.radius == rad) {

                        // first one at this radius, so this is the biggest
                        // feature at the moment
                        if (biggest_feature === undefined) {
                            biggest_feature = f;
                            biggest_size = f.bp_size();
                            furthest_point = f.end();
                            if (f.end() < f.start()) { furthest_point+=seq_length; }
                            continue;
                        }

                        new_size = f.bp_size();
                        overlap = furthest_point-f.start();

                        if (overlap <= min_overlap_cutoff) {
                            // We've cleared all potential conflicts: reset
                            // the indicators
                            biggest_size = new_size;
                            biggest_feature = f;
                            furthest_point = f.end();
                            if (f.end() < f.start()) { furthest_point+=seq_length; }

                        // since we go around twice, it is now possible
                        // for a feature to "conflict with itself," so we
                        // explicitly prevent this
                        } else if (biggest_feature &&
                                   biggest_feature != f &&
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
                                if (_debug) {
                                    console.warn(fx + ", " + rx + ": overlap=" + overlap);
                                }

                                // Update the new top dog
                                biggest_size = new_size;
                                biggest_feature = f;
                                furthest_point = f.end();
                                if (f.end() < f.start()) { furthest_point+=seq_length; }

                            } else { // The original feature is top dog. move the new
                                     // feature to the new radius
                                push(biggest_feature, f);
                            }

                        }
                    }
                }

                // Keep track of the biggest radius reached
                if (rad > max_rad) {
                    max_rad = rad;
                }

                // Move on to the next radius
                rad = new_rad;
                rx++;

            } while (conflicts > 0 && max_tries > 0); // Keep adding levels of resolution

            return max_rad;


        };

        // Global label list: keeps track of which label should be in each
        // label list, and where the label line should point to, so we can
        // use this information to figure out where to put the label and
        // minimize the number of lines that intersect.
        var n_sections = 8;
        var label_f_c = new Array(n_sections);
        var label_list_pos = new Array(n_sections);

        thi$.set_label_lists = function () {
            var lx, fx;
            // Global: keeps track of feature centers for each label
            // list, we need this to compute exactly where a label
            // should be within a label list, so to minimize
            // intersecting lines.
            for (lx = 0; lx < n_sections; lx++) {
                label_f_c[lx] = [];
            }

            for (fx = thi$.features.length - 1; fx >= 0; fx--) {
                thi$.features[fx].set_label_list();
            }
        };

        thi$.draw_labels = function() {

            // lower x, y starting position for each label list
            label_list_pos = [[0,0], [0,0], [0,0], [0,0],
                              [0,0], [0,0], [0,0], [0,0]];

            // Iterate counterclockwise, first get counts
            // Sort feature center list for each label list, and also
            // figure out where each label list should start
            for (var i=0; i<label_f_c.length; i++) {
                label_f_c[i].sort(function (a,b) {return (a[0]-b[0]);});
                var section_angle = thi$.label_list_section_angle(i);
                var xy1 = convert.polar_to_rect(thi$.label_pos, section_angle);

                // for each section, we also shift the x coordinate to
                // be further away from the circle, so that as much as
                // possible, the angle between a) line from the label
                // to the feature and b) the label is more than 90
                // degrees (i.e. visually, you don't have lines going
                // "backward").
                //
                // we also compute the lower y coordinate of each
                // label list below.
                if (i === 0 || i == 1) {
                    xy1.x += thi$.x_shift_on_labels;
                }
                else if (i == 2 || i == 3) {
                    xy1.y += label_f_c[i].length*label_letter_height;
                    xy1.x += thi$.x_shift_on_labels;
                }
                else if (i == 4 || i == 5) {
                    xy1.y += label_f_c[i].length*label_letter_height;
                    xy1.x -= thi$.x_shift_on_labels;
                }
                else if (i == 6 || i == 7) {
                    xy1.x -= thi$.x_shift_on_labels;
                }
                label_list_pos[i][0] = xy1.x;
                label_list_pos[i][1] = xy1.y;
            }

            // Finally draw labels
            for (var fx = thi$.features.length - 1; fx >= 0; fx--) {
                thi$.features[fx].draw_label();
            }
        };

        thi$.set_bounding_box = function () {
            // Figure out outter edge of label lists
            //
            // Just an educated guess based on 13pt font. we will use
            // this to compute height of label lists. These are
            // conservative.
            label_letter_height = 15;
            label_letter_width = 12;

            var min_x = thi$.width/2;
            var max_x = thi$.width/2;
            var min_y = thi$.width/2;
            var max_y = thi$.width/2;
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
                var section_angle = thi$.label_list_section_angle(section);
                var xy = convert.polar_to_rect(thi$.label_pos,section_angle);

                if (section === 0 || section == 1) { xy.x += thi$.x_shift_on_labels; }
                else if (section == 2 || section == 3) { xy.x += thi$.x_shift_on_labels; }
                else if (section == 4 || section == 5) { xy.x -= thi$.x_shift_on_labels; }
                else if (section == 6 || section == 7) { xy.x -= thi$.x_shift_on_labels; }

                if (section === 0 || section == 1) {
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

            if (max_x < cx+outer_radius+thi$.label_offset) {
                max_x = cx+outer_radius+thi$.label_offset;
            }
            if (min_x > cx-outer_radius-thi$.label_offset) {
                min_x = cx-outer_radius-thi$.label_offset;
            }
            if (max_y < cy+outer_radius+thi$.label_offset) {
                max_y = cy+outer_radius+thi$.label_offset;
            }
            if (min_y > cy-outer_radius-thi$.label_offset) {
                min_y = cy-outer_radius-thi$.label_offset;
            }

            // Now we have a new bounding box: min_x,min_y to max_x,max_y

            var right_x_extend = max_x-cx;
            var left_x_extend = cx-min_x;
            var top_y_extend = cy-min_y;
            var bot_y_extend = max_y-cy;
            var bb_width = max_x-min_x;
            var bb_height = max_y-min_y;

            thi$.width = bb_width;
            thi$.height = bb_height;
            cx = left_x_extend;
            cy = top_y_extend;
        };

        thi$.rescale = function() { // Rescale
            if (thi$.final_width != thi$.width ||
                thi$.final_height != thi$.height) {
                // "center" parameter just adds unnecessary CSS to the container
                // object to give it an absolute position: not what we need
                this.paper.changeSize(thi$.final_width,thi$.final_height,false,false);
            }
        };

        return thi$;
    } // End CircMap()

    ///////////////////////////////////////////////////////////////////
    // Linear Map Drawing Class
    function LinMap(options) {

        // Inherit the common Map functions
        var thi$ = Object.create(new Map(options));
        thi$.FeatureType = LinearFeature;

        // Paper setup - not the final width, but how we will draw the
        // map, we will scale later on
        var cx = thi$.width/2;
        var cy = thi$.height/2;

        var plasmid_y = cy;
        var plasmid_width = thi$.width * 0.9;
        var plasmid_left = (thi$.width - plasmid_width) / 2;
        var plasmid_right = plasmid_left + plasmid_width;

        // Where to draw the map
        thi$.map_dom_id = 'giraffe-draw-map';
        if ('map_dom_id' in options) {
            thi$.map_dom_id = options.map_dom_id;
        }

        // Heights of levels
        var y_spacing = 20; // spacing
        thi$.label_offset = 50;
        if ('label_offset' in options) {
            thi$.label_offset = parseInt(options.label_offset, 10);
        }

        var head_width = 25;
        var head_length = 5;

        // Overlaps
        var min_overlap_cutoff = -1;// in pixels
        var min_overlap_pct = 0;
        var min_overlap_feature_size = 0; // in pixels

        var label_letter_height = 0;
        var label_letter_width = 0;

        var convert = {
            pos_to_x: function (p) {
                return plasmid_left + (p/seq_length) * plasmid_width;
            }
        };

        ///////////////////////////////////////////////////////////////////
        // Circular Feature class
        function LinearFeature(basic_feature) {

            // Create a prototypal descendent of the basic_feature to expand
            var thi$ = Object.create(basic_feature);

            // The result of this function will be a LinearFeature object


            thi$.y = 0; // default to plasmid height (y is an offset)

            thi$.real_center = function() {
                return (this.real_start() + this.real_end()) / 2.0;
            };

            // Degree conversion, for overlap calculation:
            // for these functions, the sequence starts at 90 degrees and goes down.
            thi$.real_start = function() {
                var rs;
                // Take the minimum head size into account. Only need to do this
                // when the head is drawn and pointing clockwise, to
                // "push the start back."
                if (this.draw_head && this.clockwise()) {
                    rs = convert.pos_to_x(this.end()) - this.real_size();
                } else { // Headless feature, or head is pointing the wrong way.
                         // Just give its typical start position
                    rs = convert.pos_to_x(this.start());
                }
                return rs;
            };

            thi$.real_end = function() {
                var re;
                // Take the minimum head size into account. Only need to do this
                // when the head is drawn and pointing counterclockwise, to
                // "push the end forward."
                if (this.draw_head && !this.clockwise()) { // Take the minimum head size into account
                    re = convert.pos_to_x(this.start()) + this.real_size();
                    if (re > plasmid_width+plasmid_left) { re = plasmid_left+plasmid_width; }
                } else { // Headless feature, or head is pointing the wrong way.
                         // Just give its typical end position
                    re = convert.pos_to_x(this.end());
                }

                return re;
            };

            thi$.real_size = function() {
                var rsz;
                // Normal definition of size
                rsz = convert.pos_to_x(thi$.end()) -
                      convert.pos_to_x(thi$.start());

                // Head size: return this if it's bigger
                if (this.draw_head && head_length > rsz) {
                        rsz = head_length;
                }

                return rsz;
            };


            thi$.draw = function () {
                var x0, x1, y,
                    // Arrowhead drawing
                    hx_tip, hx_back, head,
                    // Body drawing
                    body, x_m;


                // Don't draw features that cross the boundary, as this is not
                // a circular plasmid
                if (! this.visible || this.crosses_boundary()) { return; }

                // Convert from sequence positions to x-coords
                x0 = convert.pos_to_x(this.start());
                x1 = convert.pos_to_x(this.end());

                y = plasmid_y + this.y;

                // Arrowhead drawing, if needed
                if (this.draw_head) {
                    if (this.clockwise()) {
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
                    head = this.map.paper.path(svg.move(hx_tip, y) +
                                     svg.line(hx_back, y - head_width/2.0) +
                                     svg.line(hx_back, y + head_width/2.0) +
                                     svg.close());
                    head.attr({"stroke-width": 0,
                               "fill":         this.color});
                    this.arrow_set.push(head);
                }

                // Body drawing
                if (x0 < x1 && this.type() != ft.enzyme) {
                    // Compensating for the head may have "taken up" all
                    // the room on the plasmid, in which case no arc needs
                    // to be drawn

                    // The body has no fill-color: it's just a thick line
                    body = this.map.paper.path(svg.move(x0, y) +
                                          svg.line(x1, y));
                    body.attr({"stroke-width": this.width});

                    this.arrow_set.push(body);
                } else if (this.type() == ft.enzyme) {
                    // Restriction enzymes get drawn on their own
                    x_m = convert.pos_to_x(this.cut());

                    body = this.map.paper.path(svg.move(x_m, y - this.width/2.0) +
                                          svg.line(x_m, y + this.width/2.0));
                    body.attr({"stroke-width": this.map.enzyme_weight()});
                    body.toBack();

                    this.arrow_set.push(body);
                }

                this.arrow_set.click(this.click.bind(this));
                this.arrow_set.hover(this.mouse_over.bind(this), this.mouse_up.bind(this));

                this.feature_set.push(this.arrow_set);

                // Apply the feature-wide properties to the whole feature
                this.feature_set.attr({"stroke": this.color,
                                   "stroke-linecap": "butt",
                                   "opacity": this.opacity});

                if (this.map.is_digest() && this.type() != ft.enzyme) {
                    this.feature_set.attr("title", this.label_name());
                }

            }; // END LinearFeature::draw()

            // Should we draw the label?
            thi$.should_draw_label = function () {
                return this.visible && this.labeled && !this.crosses_boundary();
            }; // END LinearFeature::should_draw_label()

            thi$.label_size = function() {
                var fake_label = this.map.paper.text(0, 0, this.label_name());
                var w = fake_label.getBBox().width;
                var h = fake_label.getBBox().height;
                fake_label.remove();

                return { width: w, height: h };
            };

            thi$.draw_label = function (height, pos) {
                if (!this.should_draw_label()) { return undefined; }

                if (this.label_drawn) {
                    this.clear_label();
                }

                // Figure out the center of the feature
                var x_c;
                if (this.type() == ft.enzyme) {
                    x_c = convert.pos_to_x(this.cut());
                } else {
                    x_c = this.real_center();
                }

                // Enzymes show their cut sites in the label
                var label_name = thi$.label_name();
                var label = this.map.paper.text(pos, height, label_name);

                // Below, right-justify. Above, left-justify.
                var anchor = (height >= plasmid_y) ? "end" : "start";
                label.attr({"fill": this.color,
                            "text-anchor": anchor,
                            "font-size": this.map.label_font_size(),
                            "opacity": 1.0 });

                // Draw the line to the label position
                var label_line = this.map.paper.path(svg.move(x_c, plasmid_y + this.y) +
                                            svg.line(pos, height));
                label_line.attr({"stroke": colors.bg_text,
                                 "stroke-width": this.map.label_line_weight(),
                                 "opacity": 0.5 });

                this.label_set.push(label_line);
                this.label_set.push(label);

                // Handlers
                this.label_set.click(this.click.bind(this));
                this.label_set.hover(this.mouse_over.bind(this), this.mouse_up.bind(this));

                // Only push label_line, so when we fade in and out,
                // we don't also fade the label.
                this.feature_set.push(label_line);

                this.labeled = true;
                this.label_drawn = true;

                return label.getBBox();
            }; // END LinearFeature::draw_label()

            return thi$;
        } // END LinearFeature Class

        thi$.draw_plasmid = function() {

            // Title
            var title_y = 1.5 * y_spacing;

            // Tic marks
            var tic_mark_length = 15;
            var tic_mark_y = plasmid_y + 2* y_spacing;
            var tic_label_y = tic_mark_y + 1.5*tic_mark_length;

            function draw_tic_mark(p,t) {
                var x = convert.pos_to_x(p);
                var y0 = tic_mark_y - tic_mark_length/2;
                var y1 = tic_mark_y + tic_mark_length/2;
                var tic = thi$.paper.path(svg.move(x, y0) +
                                     svg.line(x, y1));
                tic.attr({"stroke": colors.bg_text});

                var label = thi$.paper.text(x, tic_label_y, String(t));
                label.attr({"fill": colors.bg_text});
            }

            var plasmid = this.paper.path(svg.move(plasmid_left,  plasmid_y) +
                                     svg.line(plasmid_right, plasmid_y));

            plasmid.attr("stroke", colors.plasmid);
            var title = seq_length + ' bp';
            if (this.region_start_offset() != 0) { title = (1+this.region_start_offset())+' - '+(this.region_start_offset()+seq_length); }
            if (this.plasmid_name() !== "") {
                title = this.plasmid_name() + ": " + title;
            }
            var plasmid_label = this.paper.text(cx, title_y, title);
            plasmid_label.attr({"fill":      colors.plasmid,
                                "font-size": this.plasmid_font_size() });

            // Set the scale to be the order of magnitude of seq_length
            // i.e. 100, 1000, 10, etc.
            if (this.draw_tic_mark()) {
                var scale = Math.pow(10, Math.floor(Math.log(seq_length)/Math.log(10)));
                draw_tic_mark(1,this.region_start_offset()+1);
                for (var xx = scale; xx <= seq_length; xx += scale) {
                    draw_tic_mark(xx,xx+this.region_start_offset());
                }
                draw_tic_mark(seq_length,seq_length+this.region_start_offset());
            }
        };

        thi$.resolve_conflicts = function () {
            var conflicts,
                y = 0, // current radius
                yx = 1,               // radius counter
                max_dist = 0, new_dist, new_y,
                biggest_size, biggest_feature, furthest_point,
                new_size, overlap,
                fx, f;

            // Utility functions
            function push(winner, loser) {
                var ppfx, pfx, pf,
                    can_push;


                // Record that the push happened
                winner.pushed_features.push(loser);
                conflicts++;

                // Do it
                loser.y = new_y;

                if (_debug) {
                    console.warn(loser.name() + " pushed by " + winner.name());
                }

                // Since loser was pushed, un-push all the
                // features that it pushed, as long as
                // those features are not in conflict with the winner,
                // or with their own, previously pushed features, which are
                // now unpushed
                for (pfx = 0; pfx < loser.pushed_features.length; pfx++) {
                    pf = loser.pushed_features[pfx];
                    // Check for conflict with the winner feature itself
                    // If there's no conflict, we can push it back safely.
                    if (winner.real_end() - pf.real_start() <= min_overlap_cutoff ||
                        pf.real_end() - winner.real_start() <= min_overlap_cutoff) {

                        // Check for conflict with previously pushed features
                        // that may have been unpushed
                        can_push = true;
                        for (ppfx = 0; ppfx < pf.pushed_features.length; ppfx++) {
                            if (pf.pushed_features[ppfx].y == y) {
                                can_push = false;
                                break;
                            }
                        }

                        // Finally!
                        if (can_push) {
                            if (_debug) {
                                console.warn(pf.name() + " unpushed, because " +
                                             loser.name() + " pushed by " + winner.name());
                            }
                            pf.y = y;
                        }
                    }
                }
            }

            // Main body

            // reset radii in case this is being done any time but the first
            for (fx = 0; fx < this.features.length; fx++) {
                this.features[fx].y = 0;
            }

            var max_tries = 11;

            do {
                max_tries--;

                // Keep alternating between inside and outside the plasmid.
                new_y = y + Math.pow(-1, yx) * yx * y_spacing;

                conflicts = 0; // Assume you have no conflicts until you find some

                // Clear the record of who pushed whom
                for (fx = 0; fx < this.features.length; fx++) {
                    this.features[fx].pushed_features = [];
                }

                biggest_size = 0;
                biggest_feature = undefined;
                furthest_point = plasmid_left; // Start at a complete lower bound

                for (fx = 0; fx < this.features.length; fx++) {
                    f = this.features[fx];
                    if (f.y == y && f.type() != ft.enzyme && f.visible) {
                        new_size = f.real_size();
                        // first feature in this level, so this is the biggest
                        // one at this level at the moment
                        if (biggest_feature === undefined) {
                            biggest_size = new_size;
                            biggest_feature = f;
                            furthest_point = f.real_end();
                            continue;
                        }
                        overlap = furthest_point - f.real_start();
                        if (overlap <= min_overlap_cutoff) {
                            // We've cleared all potential conflicts: reset
                            // the indicators
                            biggest_size = new_size;
                            biggest_feature = f;
                            furthest_point = f.real_end();
                        // explicitly prevent conflicts with self
                        } else if (biggest_feature &&
                                   biggest_feature != f &&
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
                new_dist = Math.abs(y);
                if (new_dist > max_dist) {
                    max_dist = new_dist;
                }

                // Move on to the next radius
                y = new_y;
                yx++;

            } while (conflicts > 0 && max_tries > 0); // Keep adding levels of resolution

            return max_dist;

        };

        var label_pos, label_lists;

        thi$.set_label_lists = function () {
            var label_overlap_cutoff = -1, // pixel
                nlists = 6,
                ix, fx, f, feature_center,
                section, bottom,
                list_offset;

            //                   top                bottom
            label_pos    = [ new Array(nlists), new Array(nlists)];
            label_lists  = [ new Array(nlists), new Array(nlists)];

            // Initialize
            for (ix = 0; ix < nlists; ix++) {
                for (lx = 0; lx < 2; lx++) {
                    label_lists[lx][ix] = [];
                }
            }

            // Figure out which list each feature goes in
            for (fx = 0; fx < thi$.features.length; fx++) {
                f = thi$.features[fx];

                if (f.type() == ft.enzyme) {
                    // Restriction enzymes get drawn on their own
                    feature_center = convert.pos_to_x(f.cut());
                } else {
                    feature_center = f.real_center();
                }

                // Which nth of the plasmid is the feature in?
                section = Math.floor(nlists*(feature_center - plasmid_left)/
                                                 plasmid_width);
                // Is it in the top or bottom?
                bottom = section % 2;

                if (f.should_draw_label()) {
                    // push it on the appropriate list
                    label_lists[bottom][section].push(f);
                }
            }

            // Calculate positions of label lists
            //                 top  bottom
            list_offset = [20, -20];
            for (ix = 0; ix < nlists; ix++) {
                if (ix % 2 === 0) {

                    if (label_lists[0][ix].length >= 1) {
                        // Top label: just to the right of the last feature
                        label_pos[0][ix] =
                            label_lists[0][ix][label_lists[0][ix].length - 1].real_end() +
                            list_offset[0];
                    } else { // Empty list: just put it somewhere
                        label_pos[0][ix] = plasmid_left +
                            ix * plasmid_width / nlists + list_offset[0];
                    }
                } else {
                    if (label_lists[1][ix].length >= 1) {
                        // Bottom label: just to the left of the first feature
                        label_pos[1][ix] =
                            label_lists[1][ix][0].real_start() +
                            list_offset[1];
                    } else { // Empty list: just put it somewhere
                        label_pos[1][ix] = plasmid_left +
                            (ix + 0.5) * plasmid_width / nlists + list_offset[1];
                    }
                }
            }

        };

        thi$.draw_labels = function () {

            var label_height = 13; // Just an educated guess
            var label_leading = 1.3;

            // Finally, draw all the features
            for (var sx = 0; sx < label_pos[0].length; sx++) {
                for (var lx = 0; lx < 2; lx++) {
                    var ll = label_lists[lx][sx];

                    // Sort the list by center position
                    var comp_factor = 0.55;
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
                    });

                    var num_labels = ll.length;

                    // Iterate over every label in the list
                    for (var ix = 0; ix < num_labels; ix++) {
                        var curr_height;
                        if (lx) { // Bottom list: top to bottom
                            curr_height = plasmid_y + thi$.label_pos +
                                (ix) * label_height;
                        } else { // Top list: bottom to top
                            curr_height = plasmid_y - thi$.label_pos -
                                (num_labels - 1 - ix) * label_height;
                        }
                        ll[ix].draw_label(curr_height, label_pos[lx][sx]);
                    }
                }
            }
        };

        thi$.set_bounding_box = function () {
            var ix, lx, ll, sx;

            // Figure out outer edge of label lists
            //
            // Just an educated guess based on 13pt font. we will use
            // this to compute height of label lists. These are
            // conservative.
            label_letter_height = 17;
            label_letter_width = 8;

            var min_y = plasmid_y;
            var max_y = plasmid_y;
            // By default, plasmid will scale to width of map, so unless we
            // actually have lists that go off the page, no reason to adjust
            // them. i.e., never make them smaller than they orignially were,
            // because there is never a need to "zoom in."
            var min_x = 0;
            var max_x = thi$.width;

            // Iterate over every list in a level
            for (sx = 0; sx < label_pos[0].length; sx++) {
                for (lx = 0; lx < 2; lx++) {
                    ll = label_lists[lx][sx];

                    var list_max_letters = 0;
                    for (ix = 0; ix < ll.length; ix++) {
                        var num_letts = ll[ix].label_name().length;
                        if (num_letts > list_max_letters) {
                            list_max_letters = num_letts;
                        }
                    }

                    if (lx === 0) { // Top lists: move top and right
                        var list_top = plasmid_y - thi$.label_pos -
                            label_letter_height * (ll.length + 1);
                        if (list_top < min_y) {
                            min_y = list_top;
                        }
                        var list_right =  label_pos[lx][sx] +
                            label_letter_width * list_max_letters;
                        if (list_right > max_x) {
                            max_x = list_right;
                        }

                    } else if (lx == 1) { // Bot lists: move bot and left
                        var list_bot = plasmid_y + thi$.label_pos +
                            label_letter_height * (ll.length + 1);
                        if (list_bot > max_y) {
                            max_y = list_bot;
                        }
                        var list_left =  label_pos[lx][sx] -
                            label_letter_width * list_max_letters;
                        if (list_left < min_x) {
                            min_x = list_left;
                        }
                    }
                }
            }

            // Now we have a new bounding box (height only): min_y to max_y

            // Extend or compress the box dimensions to encompas this new size
            thi$.width = max_x - min_x;
            thi$.height = max_y - min_y;

            // Shift all the reference points to compensate for the re-zooming
            cy -= min_y;
            cx -= min_x;
            plasmid_y = cy;
            plasmid_left -= min_x;
            plasmid_right -= min_x;

            // Shift the label positions as well, to compensate
            for (lx = 0; lx < label_pos.length; lx++) {
                for (ix = 0; ix < label_pos[lx].length; ix++) {
                    label_pos[lx][ix] -= min_x;
                }
            }

        };

        thi$.rescale = function() {
            // Rescale
            if (thi$.final_width != thi$.width ||
                thi$.final_height != thi$.height) {

                // Make sure not to add additional height to the map, once we've
                // trimmed it off
                thi$.final_height = thi$.final_width * (thi$.height/thi$.width);

                // "center" parameter just adds unnecessary CSS to the container
                // object to give it an absolute position: not what we need
                this.paper.changeSize(thi$.final_width,thi$.final_height,false,true);
            }
        };

        return thi$;
    } // END LinMap()

    return this;
};// END GiraffeDraw
})(); // END GiraffeDraw Namespace

// vi: set expandtab:ts=4:sw=4:sts=4
