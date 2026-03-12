(function () {
  var toggle = document.getElementById('dark-mode-toggle');
  if (!toggle) return;

  toggle.addEventListener('click', function () {
    var html = document.documentElement;
    var current = html.getAttribute('data-theme');
    var next = current === 'dark' ? 'light' : 'dark';
    html.setAttribute('data-theme', next);
    localStorage.setItem('theme', next);
  });
})();
