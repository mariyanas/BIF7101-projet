from flask import Flask, render_template, request, Response, url_for, request, render_template, send_file, send_from_directory
from Bio import Phylo, AlignIO
from io import StringIO
from werkzeug.utils import secure_filename
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

import matplotlib
import matplotlib.pyplot as plt
import subprocess, time, threading, os
import os
import subprocess
import uuid

matplotlib.use("Agg")

app = Flask(__name__)

UPLOAD_FOLDER = "uploads"
os.makedirs(UPLOAD_FOLDER, exist_ok=True)

@app.route('/uploads/<filename>')
def uploaded_file(filename):
    return send_from_directory(UPLOAD_FOLDER, filename)

# TOOL: Newick Tree Viewer
@app.route("/", methods=["GET", "POST"])
def tree_viewer():
    tree_image = None
    tree_error = None
    active_tab = request.form.get("active_tab") or request.args.get("tab") or "main-page"

    if request.method == "POST":
        newick_str = request.form.get("newick")
        active_tab = "newick_tree_viewer"  # keep tab active
        
        width = float(request.form.get("width") or 6.4)
        height = float(request.form.get("height") or 4.8)

        if not newick_str or newick_str.strip() == "":
            tree_error = "Please provide a Newick string."
        else:
            try:
                handle = StringIO(newick_str)
                tree = Phylo.read(handle, "newick")
                fig = plt.figure(figsize=(width, height))
                ax = fig.add_subplot(1, 1, 1)
                Phylo.draw(tree, axes=ax, do_show=False)

                filename = f"tree_{uuid.uuid4().hex}.png"
                image_path = os.path.join(UPLOAD_FOLDER, filename)
                plt.savefig(image_path, bbox_inches="tight")
                plt.close()

                tree_image = url_for("uploaded_file", filename=filename)
                print(os.path.exists(image_path))

            except Exception as e:
                tree_error = f"Failed to generate tree: {str(e)}"

    return render_template("index.html", tree_image=tree_image, tree_error=tree_error, active_tab=active_tab)

# TOOL: Conversion
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

# https://ena01.uqam.ca/pluginfile.php/9775398/mod_resource/content/2/analyse_phylogenetique_alignement_distances_arbres.html?embed=1
def parse_iqtree_distance_matrix(file_path):
    with open(file_path, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
    
    n = int(lines[0])
    names = []
    matrix = []
    
    for i in range(1, n + 1):
        parts = lines[i].split()
        names.append(parts[0])
        distances = []
        for j in range(1, i + 1):
            distances.append(float(parts[j]))
        
        matrix.append(distances)
    
    return DistanceMatrix(names, matrix)

# TOOL: Distance (Neighbour-Joining and UPGMA)
@app.route("/distance-methods", methods=["POST"])
def NJ():
    active_tab = "distance"
    file = request.files.get("fasta_file_distance")
    method = request.form.get("method")
    
    if not file or file.filename == "":
        return render_template(
            "index.html",
            error="Please upload a file.",
            active_tab="distance"
        )
    
    try:
        unique_folder = os.path.join(UPLOAD_FOLDER, str(uuid.uuid4()))
        os.makedirs(unique_folder, exist_ok=True)

        safe_filename = secure_filename(file.filename)
        filepath = os.path.join(unique_folder, safe_filename)

        file.save(filepath)
        print("Saved file to:", filepath)

        # run IQ-TREE to find best model
        cmd = [
            "iqtree3",
            "-s", filepath,
            "-m", "MFP",
            "-nt", "AUTO",
            "-redo"
        ]     
        # print("Running command:", " ".join(cmd))

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True
        )

        # print("STDOUT:", result.stdout)
        # print("STDERR:", result.stderr)
        
        # read .iqtree file to get the best model
        iqtree_file = filepath + ".iqtree"
        best_model = None

        if os.path.exists(iqtree_file):
            with open(iqtree_file, 'r', encoding='utf-8') as f:
                for line in f:
                    if "Best-fit model according to BIC:" in line:
                        best_model = line.split(":")[1].strip()
                        break
        else :
            return render_template(
                    "index.html",
                    error="Could not find best-fit model",
                    active_tab="distance"
                )
        
        prefix = "best_model_tree"
        
        model_cmd = [
            "iqtree3",
            "-s", filepath,
            "-m", best_model,
            "-nt", "AUTO",
            "-redo", 
            "-pre", os.path.join(os.path.dirname(filepath), prefix)
        ]
        
        result2 = subprocess.run(
            model_cmd,
            capture_output=True,
            text=True
        )
        
        # print("STDOUT:", result2.stdout)
        # print("STDERR:", result2.stderr)
        
        mldist_file = prefix + ".mldist"
        filepath_mldist_file = os.path.join(unique_folder, mldist_file)
        
        if not os.path.exists(filepath_mldist_file):
            return render_template(
                "index.html",
                error="Distance matrix (.mldist) not found.",
                active_tab="distance"
            )

        dm = parse_iqtree_distance_matrix(filepath_mldist_file)

        constructor = DistanceTreeConstructor()
        if method == "nj":
            tree_distance = constructor.nj(dm)
            image_name = "tree_nj.png"
        elif method == "upgma":
            tree_distance = constructor.upgma(dm)
            image_name = "tree_upgma.png"
        else:
            raise ValueError("Invalid method selected")

        nj_tree_file = os.path.join(unique_folder, "nj_tree.nwk")
        Phylo.write(tree_distance, nj_tree_file, "newick")
        
        with open(nj_tree_file, "r") as f:
            newick_str = f.read().strip()
            
        # TODO: Add lengths, rename root and nodes
        tree = Phylo.read(StringIO(newick_str), "newick")
        num_leaves = tree.count_terminals()
        set_height = max(5, num_leaves * 1.0)
        fig = plt.figure(figsize=(10, set_height))
        ax = fig.add_subplot(1, 1, 1)
        Phylo.draw(tree, axes=ax, do_show=False)
        
        # TODO: eventually to be put into unique_folder
        image_path = os.path.join("static", image_name)
        plt.savefig(image_path, bbox_inches="tight")
        plt.close()

        tree_image_distance = url_for("static", filename=image_name)

    except Exception as e:
        return render_template("index.html", distance_error=str(e), active_tab="distance")
        
    return render_template("index.html", tree_image_distance=tree_image_distance, active_tab=active_tab)
        
# TOOL: Bayesian Inference
def run_mrbayes(nexus_path, working_dir):
    try:
        print("MrBayes start")
        result = subprocess.run(
            ["mb", os.path.basename(nexus_path)],
            cwd=working_dir,
            capture_output=True,
            text=True
        )
        print("MrBayes end")

        log_path = os.path.join(working_dir, "mrbayes_log.txt")
        with open(log_path, "w") as f:
            f.write(result.stdout + "\n\n" + result.stderr)

    except Exception as e:
        print("MrBayes error:", str(e))
        
        
@app.route("/bayesian_inference", methods=["POST"])
def bayesian_inference():
    file = request.files.get("sequence_file")
    
    print("File object:", file)
    print("Filename:", file.filename if file else "No file")

    print("ngen:", request.form.get("ngen"))
    print("samplefreq:", request.form.get("samplefreq"))
    print("printfreq:", request.form.get("printfreq"))
    print("burnin:", request.form.get("burnin"))

    if not file or file.filename == "":
        return render_template(
            "index.html",
            error="Please upload a file.",
            active_tab="bayesian-inference"
        )

    try:
        filename = file.filename.lower()
        file_text = file.read().decode("utf-8")

        # Default MrBayes parameters
        ngen = request.form.get("ngen", 1000000)
        samplefreq = request.form.get("samplefreq", 100)
        printfreq = request.form.get("printfreq", 100000)
        burnin = request.form.get("burnin", 2500)

        if filename.endswith(".nex"):
            nexus_text = file_text
        else:
            return render_template(
                "index.html",
                error="Unsupported file type.",
                active_tab="bayesian-inference"
            )

        # Add MrBayes block
        nexus_text += "\nBEGIN MRBAYES;\n"
        nexus_text += "  set autoclose=yes nowarn=yes;\n"
        nexus_text += "  lset nst=6 rates=invgamma;\n"
        nexus_text += f"  mcmc ngen={ngen} samplefreq={samplefreq} printfreq={printfreq} nchains=4 nruns=2;\n"
        nexus_text += f"  sumt burnin={burnin};\n"
        nexus_text += f"  sump burnin={burnin};\n"
        nexus_text += "END;\n"

        # Save file
        # nexus_filename = "static/mrbayes_input.nex"
        
        unique_folder = os.path.join(UPLOAD_FOLDER, str(uuid.uuid4()))
        os.makedirs(unique_folder, exist_ok=True)

        nexus_filename = os.path.join(unique_folder, "mrbayes.nex")
        
        with open(nexus_filename, "w") as f:
            f.write(nexus_text)

        # Run MrBayes in background
        threading.Thread(target=run_mrbayes, args=(nexus_filename, unique_folder)).start()
        
        # run_mrbayes(nexus_filename, unique_folder)

        return render_template(
            "index.html",
            message="MrBayes is running... refresh in a few seconds.",
            active_tab="bayesian-inference"
        )

    except Exception as e:
        return render_template(
            "index.html",
            error=f"Error: {str(e)}",
            active_tab="bayesian-inference"
        )

# if __name__ == "__main__":
#     app.run(debug=True)