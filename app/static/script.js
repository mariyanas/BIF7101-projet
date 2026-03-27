function showTool(toolId, button) {
    const tools = document.querySelectorAll('.tool');
    tools.forEach(t => t.classList.remove('active'));
    const buttons = document.querySelectorAll('.sidebar button');
    buttons.forEach(b => b.classList.remove('active'));

    document.getElementById(toolId).classList.add('active');
    if (button) button.classList.add('active');
}

