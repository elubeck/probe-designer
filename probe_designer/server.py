from Bio import Entrez
from flask import Flask, render_template, request, make_response
from flask.ext.basicauth import BasicAuth
from flask.ext.wtf import Form
from flask_bootstrap import Bootstrap
from wtforms import validators, SelectField, BooleanField, IntegerField,\
    SubmitField, FloatField, StringField, TextAreaField

from probe_designer.mRNA_designer import design_step_gui
from probe_designer.probe_refiner import ProbeFilter

Entrez.email = 'elubeck@caltech.edu'

app = Flask(__name__)
app.config['SECRET_KEY'] = 'bodole'
app.config['RECAPTCHA_PUBLIC_KEY'] = \
    '6Lfol9cSAAAAADAkodaYldddd22414141'
BasicAuth(app)
app.config['BASIC_AUTH_USERNAME'] = 'cailab'
app.config['BASIC_AUTH_PASSWORD'] = 'all your probes are belong to us'
app.config['BASIC_AUTH_FORCE'] = True
# AppConfig(app, None)
Bootstrap(app)


# Model
class InputForm(Form):
    std_validation = [validators.InputRequired()]
    non_neg_val = [validators.InputRequired(), validators.NumberRange(min=0)]
    genes = TextAreaField(validators=std_validation)
    gc_target = FloatField(label='Target GC', default=0.55, validators=std_validation)
    gc_min = FloatField(label='Minimum GC', default=0.35, validators=std_validation)
    gc_max = FloatField(label='Maximum GC', default=0.75, validators=std_validation)
    length = IntegerField(label='Probe Length', default=35, validators=[validators.InputRequired(),
                                                  validators.NumberRange(min=14)])
    spacing = IntegerField(label='Probe Spacing', default=1, validators=non_neg_val)
    max_probes = IntegerField(label='Maximum Number of Probes', default=48, validators=non_neg_val)
    false_pos_len = IntegerField(label='False Positive Match Length(nt)',
                                 default=18, validators=non_neg_val)
    max_off_target = IntegerField(label='Maximum Off-Target Threshold(AU)',
                                  default=50, validators=non_neg_val)
    off_hits = IntegerField(label='Maximum Off-Target Hits on Any Gene',
                            default=6, validators=non_neg_val)
    copy_num_db = SelectField(label='Copy Number Database',
                              default='brain',
                              validators=std_validation,
                              choices=[('brain', 'brain'), ('embryonic11.5', 'embryonic day 11.5')])
    cds_only = BooleanField(label='CDS Only', default=True)
    output_name = StringField(label='Output CSV Name', default='probes.csv')
    submit_button = SubmitField('Submit Form')


def probes_2_str(probes, name):
    csv_str = ["{},{},{}".format(name, num, probe)
               for num, probe in enumerate(probes)]
    return "\n".join(csv_str)


def parse_form2(form):
    import io
    import csv
    output = io.StringIO()
    writer = csv.writer(output)
    filterer = ProbeFilter(db='gencode_tracks_reversed_introns+mRNA',
                           copy_num=form.data['copy_num_db'])
    for gene_name in form.data['genes'].split(','):
        print(gene_name)
        name = gene_name.strip(' ')
        gene_name = name[0].upper() + name[1:].lower()
        name, probes, seq = design_step_gui(gene_name, **form.data)
        if name == "FAILED":
            esearch = Entrez.read(Entrez.esearch(db='gene',
                                                 term='"{}"[gene] AND "Mus musculus"[orgn]'.format(
                                                     gene_name),
                                                 retmode='xml'))
            esearch2 = Entrez.read(Entrez.esearch(db='gene',
                                                  term='"{}" AND "Mus musculus"[orgn]'.format(
                                                      gene_name),
                                                  retmode='xml'))
            if len(esearch['IdList']) == 1:
                result = Entrez.read(Entrez.efetch(db='gene',
                                                   id=esearch['IdList'][0],
                                                   retmode='xml'))
            elif len(esearch2['IdList']) == 1:
                result = Entrez.read(Entrez.efetch(db='gene',
                                                   id=esearch2['IdList'][0],
                                                   retmode='xml'))
            else:
                continue
            proper_name = result[0]['Entrezgene_gene']['Gene-ref']['Gene-ref_locus']
            name, probes, seq = design_step_gui(proper_name, **form.data)
        if not probes:
            continue
        probes2 = filterer.run(set(probes), name, **form.data)
        for n, probe in enumerate(probes2):
            writer.writerow((gene_name, n + 1, name, probe))
    return output.getvalue()


# View
@app.route('/', methods=['GET', 'POST'])
def index():
    form = InputForm(request.form)
    if request.method == 'POST' and form.validate():
        csv = parse_form2(form)
        response = make_response(csv)
        response.headers["Content-Disposition"] = "attachment; filename={}".format(form.data['output_name'])
        return response
    return render_template("view.html", form=form, )


if __name__ == '__main__':
    app.run(host='0.0.0.0', debug=False)
