(function () {
  // Collapsible categories
  var categoryBtns = document.querySelectorAll('.sidebar-category-btn, .sidebar-toggle-btn');
  categoryBtns.forEach(function (btn) {
    btn.addEventListener('click', function () {
      var category = btn.closest('.sidebar-category');
      if (category) {
        category.classList.toggle('expanded');
        var expanded = category.classList.contains('expanded');
        btn.setAttribute('aria-expanded', expanded);
      }
    });
  });

  // Mobile sidebar toggle
  var sidebarToggle = document.getElementById('sidebar-toggle');
  var sidebar = document.getElementById('sidebar');
  var sidebarClose = document.getElementById('sidebar-close');

  if (sidebarToggle && sidebar) {
    sidebarToggle.addEventListener('click', function () {
      sidebar.classList.toggle('open');
    });
  }

  if (sidebarClose && sidebar) {
    sidebarClose.addEventListener('click', function () {
      sidebar.classList.remove('open');
    });
  }

  // Close sidebar on click outside (mobile)
  document.addEventListener('click', function (e) {
    if (sidebar && sidebar.classList.contains('open') &&
        !sidebar.contains(e.target) &&
        e.target !== sidebarToggle &&
        !sidebarToggle.contains(e.target)) {
      sidebar.classList.remove('open');
    }
  });
})();
