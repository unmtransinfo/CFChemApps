from flask import Flask

app = Flask(__name__)

@app.route("/")
def app_start():
    return "<h1>Hello</h1>"