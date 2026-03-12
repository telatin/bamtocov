(function () {
  var tocLinks = document.querySelectorAll('.toc-link');
  if (!tocLinks.length) return;

  var headings = [];
  tocLinks.forEach(function (link) {
    var id = link.getAttribute('href').slice(1);
    var el = document.getElementById(id);
    if (el) headings.push({ el: el, link: link });
  });

  if (!headings.length) return;

  var current = null;

  function updateActive() {
    var scrollTop = window.scrollY + 100;
    var active = null;

    for (var i = headings.length - 1; i >= 0; i--) {
      if (headings[i].el.offsetTop <= scrollTop) {
        active = headings[i];
        break;
      }
    }

    if (active && active !== current) {
      if (current) current.link.classList.remove('active');
      active.link.classList.add('active');
      current = active;
    } else if (!active && current) {
      current.link.classList.remove('active');
      current = null;
    }
  }

  var ticking = false;
  window.addEventListener('scroll', function () {
    if (!ticking) {
      requestAnimationFrame(function () {
        updateActive();
        ticking = false;
      });
      ticking = true;
    }
  });

  updateActive();
})();
