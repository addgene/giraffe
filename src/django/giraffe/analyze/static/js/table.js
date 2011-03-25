(function () { 
// Private utility function with no external or closure use
function type_code_to_name(code, gd) {
	for (var type in gd.Feature_Type) {
		if (code == gd.Feature_Type[type]) {
			return type[0].toUpperCase() + type.slice(1);
		}
	}

	return "";
}

window.GiraffeTable = function ($,gd,dom) {	
 	var enzyme_table,
	    feature_table,
		orf_table;

	// Read the features in for each table, only displaying the table
	// if there are features for it

	// Standard features
	if (gd.std_features.length > 0) {
		
		// Initialize table
		feature_table = $('<table></table>')
			.append('<thead><tr><th>Feature</th><th>Type</th><th>Start</th><th>End<th></tr></thead>')
			.append('<tbody></tbody>')
			.addClass('giraffe-table-feature');

		for (var fx = 0; fx < gd.std_features.length; fx++) {
			var f = gd.std_features[fx],
				row = feature_table.children('tbody').append('<tr></tr>');

			row.append('<td>' + f.name() + '</td>');
			row.append('<td>' + type_code_to_name(f.type(), gd) + '</td>');

			if (f.clockwise()) {
				row.append('<td>' + f.start() + '</td>');
				row.append('<td>' + f.end() + '</td>');
			} else {
				row.append('<td>' + f.end() + '</td>');
				row.append('<td>' + f.start() + '</td>');
			}
		}

		$(dom).append(feature_table);
	}

	// ORFs
	if (gd.orf_features.length > 0) {

		// Initialize table
		orf_table = $('<table></table>')
			.append('<thead><tr><th>ORF Frame</th><th>Start</th><th>End<th></tr></thead>')
			.append('<tbody></tbody>')
			.addClass('giraffe-table-orf');

		for (var fx = 0; fx < gd.orf_features.length; fx++) {

			var f = gd.orf_features[fx],
				row = orf_table.children('tbody').append('<tr></tr>');

			row.append('<td>' + f.name().replace(/ORF\s+[fF]rame\s+/, '') + '</td>');

			if (f.clockwise()) {
				row.append('<td>' + f.start() + '</td>');
				row.append('<td>' + f.end() + '</td>');
			} else {
				row.append('<td>' + f.end() + '</td>');
				row.append('<td>' + f.start() + '</td>');
			}
		}

		$(dom).append(orf_table);
	}

	// Enzymes
	if (gd.enzyme_features.length > 0) {

		enzyme_table = $('<table></table>')
			.append('<thead><tr><th>Enzyme</th><th>Cut</th></tr></thead>')
			.append('<tbody></tbody>')
			.addClass('giraffe-table-enzyme');

		for (var fx = 0; fx < gd.enzyme_features.length; fx++) {
			var f = gd.enzyme_features[fx],
				row;

			if (f.default_show_feature() && f.cut_count() == 1) {
				row = enzyme_table.children('tbody').append('<tr></tr>');
				row.append('<td>' + f.name() + '</td>');
				row.append('<td>' + f.cut() + '</td>');
			}
		}

		$(dom).append(enzyme_table);
	}

	// Set general appearance properties
	$(dom).children().addClass('giraffe-table');
	
}})();
