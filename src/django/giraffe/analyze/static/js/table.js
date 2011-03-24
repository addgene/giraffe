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

	// Initialize feature, ORF, and enzyme tables
 	feature_table = $('<table></table>');
 	enzyme_table = $('<table></table>');
 	orf_table = $('<table></table>');

	feature_table
		.append('<thead><tr><th>Feature</th><th>Type</th><th>Start</th><th>End<th></tr></thead>')
	orf_table
		.append('<thead><tr><th>ORF Frame</th><th>Start</th><th>End<th></tr></thead>')
	enzyme_table
		.append('<thead><tr><th>Enzyme</th><th>Cut</th></tr></thead>')

 	$(dom)
		.append(feature_table)
		.append(orf_table)
		.append(enzyme_table);

	$(dom).children().append('<tbody></tbody>');

	// Appearance and visuals
	$(dom).children().addClass('giraffe-table');
	feature_table.addClass('giraffe-table-feature');
	enzyme_table.addClass('giraffe-table-enzyme');
	orf_table.addClass('giraffe-table-orf');

	// Read the features in for each table

	// Standard features
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

	// ORFs
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

	// Enzymes
 	for (var fx = 0; fx < gd.enzyme_features.length; fx++) {
		var f = gd.enzyme_features[fx],
			row;

		if (f.default_show_feature() && f.cut_count() == 1) {
			row = enzyme_table.children('tbody').append('<tr></tr>');
			row.append('<td>' + f.name() + '</td>');
			row.append('<td>' + f.cut() + '</td>');
		}
	}


	
}})();
