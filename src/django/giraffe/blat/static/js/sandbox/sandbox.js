(function () {
	var one = 1;
	var listy = [5, 9];

	window.play = function () {
		one = 2;
		listy = [3, 4];
	};

	window.display = function () {
		document.write("<p>");
		document.write(one);
		document.write(" ");
		document.write(listy);
		document.write("</p>");
	}

	window.Something = function(origin) {
		function Clone () {};
		Clone.prototype = origin;

		var _this = new Clone();

		this.dispy = function () {
			document.write("<p>");
			document.write(this.a);
			document.write(" ");
			document.write(this.b);
			document.write("</p>");
		}

		//return _this;
	}

	Something.prototype.blart = function () {
		document.write("Blart!");
	}
	
})();

display();
play();
display();

var obj1 = { a: 23, b: 64 };
var osp = Something.prototype;
Something.prototype = obj1;
var obj2 = new Something();
Something.protptype = osp;
obj2.dispy();
obj2.blart();

