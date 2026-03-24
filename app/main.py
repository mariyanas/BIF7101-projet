from flask import Flask, render_template, request, Response, url_for
from Bio import Phylo, AlignIO
from io import StringIO

import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("Agg")

app = Flask(__name__)

# TOOL: newick_tree_viewer
@app.route("/", methods=["GET", "POST"])
def tree_viewer():
    tree_image = None
    error = None
    active_tab = request.form.get("active_tab") or request.args.get("tab") or "main-page"

    if request.method == "POST":
        newick_str = request.form.get("newick")
        active_tab = "newick_tree_viewer"  # keep tab active

        if not newick_str or newick_str.strip() == "":
            error = "Please provide a Newick string."
        else:
            try:
                handle = StringIO(newick_str)
                tree = Phylo.read(handle, "newick")

                plt.figure(figsize=(8, 6))
                Phylo.draw(tree, do_show=False)

                image_path = "static/tree.png"
                plt.savefig(image_path, bbox_inches="tight")
                plt.close()

                tree_image = url_for("static", filename="tree.png")

            except Exception as e:
                error = f"Failed to generate tree: {str(e)}"

    return render_template("index.html", tree_image=tree_image, error=error, active_tab=active_tab)

# TOOL: conversion
@app.route("/convert", methods=["POST"])
def convert():
    file = request.files.get("fasta_file")
    molecule_type = request.form.get("molecule_type")

    if not file or file.filename == "":
        return render_template("index.html", error="Please upload a FASTA file.")

    if not molecule_type:
        return render_template("index.html", error="Please select a molecule type.")

    try:
        fasta_text = file.read().decode("utf-8")

        fasta_io = StringIO(fasta_text)
        alignment = AlignIO.read(fasta_io, "fasta")

        for record in alignment:
            record.annotations["molecule_type"] = molecule_type

        nexus_io = StringIO()
        AlignIO.write(alignment, nexus_io, "nexus")

        result = nexus_io.getvalue()
        
        fasta_filename = file.filename.rsplit('.', 1)[0]
        nexus_filename = f"{fasta_filename}.nex"

        return Response(
            result,
            mimetype="text/plain",
            headers={
                "Content-Disposition": f"attachment; filename={nexus_filename}"
            }
        )
        
    except Exception as e:
        return render_template(
            "index.html",
            error=f"Conversion failed: {str(e)}"
        )

if __name__ == "__main__":
    app.run(debug=True)