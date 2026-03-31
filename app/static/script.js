function showTool(toolId, btn) {

    document.querySelectorAll('.tool').forEach(t => t.classList.remove('active'));

    document.getElementById(toolId).classList.add('active');

    document.querySelectorAll('.topnav button, .dropdown-content button').forEach(b => b.classList.remove('active'));

    btn.classList.add('active');
}