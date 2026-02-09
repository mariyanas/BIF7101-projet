import matplotlib
matplotlib.use("Agg")

from fastapi import FastAPI, Form
from fastapi.responses import HTMLResponse, Response
from fastapi.staticfiles import StaticFiles
from Bio import Phylo
from io import StringIO, BytesIO
import matplotlib.pyplot as plt

app = FastAPI()

app.mount("/static", StaticFiles(directory = "static"), name = "static")

@app.get("/", response_class = HTMLResponse)
def home():
    with open("static/index.html") as f:
        return f.read()

# @app.post("/get_name")
# def say_hello(your_name: str = Form(...)):
#     return "Hello, " + your_name + "!"

@app.post("/draw_tree")
def render_tree(newick: str = Form(...)):
    tree = Phylo.read(StringIO(newick), "newick")

    fig = plt.figure(figsize=(10, 10))
    Phylo.draw(tree, do_show=False)

    buf = BytesIO()
    plt.savefig(buf, format="png", bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)

    return Response(buf.getvalue(), media_type="image/png")