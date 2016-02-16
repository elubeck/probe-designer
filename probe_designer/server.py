from flask import Flask, render_template, request
from flask.ext.wtf import Form
from wtforms import TextField, validators, SelectField, DecimalField, BooleanField, IntegerField, SubmitField

from flask_bootstrap import Bootstrap
from flask_appconfig import AppConfig

app = Flask(__name__)
app.config['SECRET_KEY'] = 'bodole'
app.config['RECAPTCHA_PUBLIC_KEY'] = \
    '6Lfol9cSAAAAADAkodaYldddd22414141'
# AppConfig(app, None)
Bootstrap(app)

# Model
class InputForm(Form):
    std_validation = [validators.InputRequired()]
    non_neg_val = [validators.InputRequired(), validators.NumberRange(min=0)]
    genes = TextField(validators=std_validation)
    gc_target = DecimalField(default=0.55, validators=std_validation)
    gc_min = DecimalField(default=0.35, validators=std_validation)
    gc_max = DecimalField(default=0.75, validators=std_validation)
    cds_only = BooleanField(label='CDS Only', default=True)
    probe_length= IntegerField(default=35, validators=[validators.InputRequired(),
                                                       validators.NumberRange(min=14)])
    spacing = IntegerField(default=1, validators=non_neg_val)
    max_probes = IntegerField(default=48, validators=non_neg_val)
    false_pos_len = IntegerField(default=18, validators=non_neg_val)
    max_off_target = IntegerField(default=50, validators=non_neg_val)
    off_hits = IntegerField(default=6, validators=non_neg_val)
    submit_button = SubmitField('Submit Form')

# View
@app.route('/', methods=['GET', 'POST'])
def index():
    form = InputForm(request.form)
    if request.method == 'POST' and form.validate():
        genes = form.genes.data
        gc_target = form.gc_target.data
    else:
        genes=None
    return render_template("view.html", form=form, genes=genes)

if __name__ == '__main__':
    app.run(host='0.0.0.0', debug=True)