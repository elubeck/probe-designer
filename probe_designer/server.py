from flask import Flask, render_template, request
from wtforms import Form, TextField, validators, SelectField

app = Flask(__name__)

# Model
class InputForm(Form):
    genes = TextField(validators=[validators.InputRequired()])
    select = SelectField(u'Programming Language', choices=[('cpp', 'C++'), ('py', 'Python'), ('text', 'Plain Text')])

# View
@app.route('/', methods=['GET', 'POST'])
def index():
    form = InputForm(request.form)
    if request.method == 'POST' and form.validate():
        r = form.genes.data
        selection = form.select.data
        return render_template("view_output.html", form=form, s=r)
    else:
        return render_template("view_input.html", form=form)

if __name__ == '__main__':
    app.run(debug=True)