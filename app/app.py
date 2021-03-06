from flask import (
    Flask,
    request,
    Response,
    redirect,
    url_for,
    render_template,
    jsonify,
    send_from_directory,
    flash,
)
from flask_cors import CORS
import random
import sys
import json

from hetbuilder.atom_checks import check_atoms

app = Flask(__name__)

app = Flask(__name__, static_url_path="", static_folder="client/public/")

# Path for our main Svelte page
@app.route("/")
def index():
    return app.send_static_file("index.html")


# Path for all the static files (compiled JS/CSS, etc.)
@app.route("/<path:path>")
def home(path):
    return send_from_directory("client/public", path)


@app.route("/rand")
def hello():
    return str(random.randint(0, 100))


@app.route("/post", methods=["POST"])
def post():
    if request.method == "POST":
        f1 = request.files.get("lower")
        f2 = request.files.get("upper")
        print(f1, f2)
        # data = jsonify(request.get_json())

        return '{"foo": 1}'


if __name__ == "__main__":
    app.run(debug=True)
