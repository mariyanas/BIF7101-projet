// async function sendName() {
//   const your_name = document.getElementById("your_name").value;

//   const formData = new FormData();
//   formData.append("your_name", your_name);

//   const response = await fetch("/get_name", {
//     method: "POST",
//     body: formData
//   });

//   const message = await response.text();
//   document.getElementById("greeting").textContent = message;
// }

async function submitTree() {
  const newick = document.getElementById("newick").value; // gets the value from textarea id = "newick" (index.html)

  const formData = new FormData();
  formData.append("newick", newick);

  const response = await fetch("/draw_tree", {
    method: "POST",
    body: formData
  }); // sends a POST request to /draw_tree (main.py)

  const blob = await response.blob();
  const url = URL.createObjectURL(blob);

  document.getElementById("treeImage").src = url; // image to appear
}