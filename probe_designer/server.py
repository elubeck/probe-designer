from flask import Flask, render_template, request, make_response
from flask.ext.wtf import Form
from wtforms import TextField, validators, SelectField, DecimalField, \
    BooleanField, IntegerField, SubmitField, FloatField
from probe_designer.mRNA_designer import RNARetriever2, design_step_gui
from probe_designer.probe_refiner import ProbeFilter

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
    gc_target = FloatField(default=0.55, validators=std_validation)
    gc_min = FloatField(default=0.35, validators=std_validation)
    gc_max = FloatField(default=0.75, validators=std_validation)
    cds_only = BooleanField(label='CDS Only', default=True)
    length= IntegerField(default=35, validators=[validators.InputRequired(),
                                                       validators.NumberRange(min=14)])
    spacing = IntegerField(default=1, validators=non_neg_val)
    max_probes = IntegerField(default=48, validators=non_neg_val)
    false_pos_len = IntegerField(default=18, validators=non_neg_val)
    max_off_target = IntegerField(default=50, validators=non_neg_val)
    off_hits = IntegerField(default=6, validators=non_neg_val)
    copy_num_db = SelectField(default='brain',
                           validators=std_validation,
                           choices=[('brain', 'brain'), ('embryonic11.5', 'embryonic day 11.5')])
    submit_button = SubmitField('Submit Form')

def probes_2_str(probes, name):
    csv_str = ["{},{}".format(name, probe) for probe in probes]
    return "\n".join(csv_str)

def parse_form(form):
    csv = ""
    for gene_name in form.data['genes'].split(','):
        gene_name = gene_name.strip(' ')
        name, probes1, seq = design_step_gui(gene_name, **form.data)
        filterer = ProbeFilter(db='gencode_tracks_reversed_introns+mRNA',
                               copy_num=form.data['copy_num_db'])
        probes2 = filterer.run(set(probes1), name, **form.data)
        csv = "\n".join([csv, probes_2_str(probes2, gene_name)])
    return csv

# View
@app.route('/', methods=['GET', 'POST'])
def index():
    form = InputForm(request.form)
    if request.method == 'POST' and form.validate():
        csv = parse_form(form)
        response = make_response(csv)
        response.headers["Content-Disposition"] = "attachment; filename=probes.csv"
        return response
    else:
        genes=None
    return render_template("view.html", form=form, genes=genes)

if __name__ == '__main__':
    app.run(host='0.0.0.0', debug=False)