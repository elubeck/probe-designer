from flask import Flask, render_template, request
from wtforms import Form, TextField, validators, SelectField, DecimalField

app = Flask(__name__)

# Model
class InputForm(Form):
    genes = TextField(validators=[validators.InputRequired()])
    gc_target = DecimalField(default=0.55, validators=[validators.InputRequired()])

# View
@app.route('/', methods=['GET', 'POST'])
def index():
    form = InputForm(request.form)
    if request.method == 'POST' and form.validate():
        genes = form.genes.data
        gc_target = form.gc_target.data
        return render_template("view_output.html", form=form, s=r)
    else:
        return render_template("view_input.html", form=form)

if __name__ == '__main__':
    app.run(host='0.0.0.0', debug=False)