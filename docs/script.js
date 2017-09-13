var menuPos = 0, isFixed = false;
window.addEventListener('load', function(evt) {
	
	// Check whether we have a menuPos element
	menuPosEl = document.getElementById('menuPos');
	if(typeof(menuPosEl) == 'undefined') {
		return;
	}
	
	// Get position of menu
	menuPos = menuPosEl.getBoundingClientRect().top - document.getElementsByTagName('html')[0].getBoundingClientRect().top;
	
	// Fix menu when we scroll under it
	document.addEventListener('scroll', function(evt) {
		var isDown = document.body.scrollTop > menuPos || document.documentElement.scrollTop > menuPos;
		if (isDown && !isFixed) {
			isFixed = true;
			var menu = document.getElementById('menu');
			document.getElementById('menuPos').style.height = menu.offsetHeight+"px";
			menu.className = "fixed";
		}
		else if(!isDown && isFixed) {
			isFixed = false;
			document.getElementById('menuPos').style.height = "0px";
			document.getElementById('menu').className = "";
		}
	});

});
