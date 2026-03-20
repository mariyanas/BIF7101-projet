from flask import Flask, render_template, request
from Bio import Phylo
from io import StringIO

import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("Agg")

app = Flask(__name__)

# TOOL: newick_tree_viewer
@app.route("/", methods=["GET", "POST"])
def index():
    tree_image = None

    if request.method == "POST":
        newick_str = request.form["newick"]

        try:
            # parse Newick string
            handle = StringIO(newick_str)
            tree = Phylo.read(handle, "newick")

            # draw tree
            plt.figure(figsize=(8, 6))
            Phylo.draw(tree, do_show=False)

            # save image
            image_path = "static/tree.png"
            plt.savefig(image_path, bbox_inches="tight")
            plt.close()

            tree_image = image_path

        except Exception as e:
            tree_image = None
            print("Error:", e)

    return render_template("index.html", tree_image=tree_image)


if __name__ == "__main__":
    app.run(debug=True)