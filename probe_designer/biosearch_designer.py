from __future__ import print_function, with_statement, division
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.ui import Select
import pandas as pd
import subprocess
import dataset
import arrow
from random import getrandbits

class Biosearch(object):
    def close(self):
        self.driver.close()
        self.selenium.kill()

    def check_db(self, gene, masking):
        res = self.table.find(Name=gene, Variants=self.variants, Masking=masking)
        for item in res:
            if arrow.get(item['Date']) >= self.date.replace(years = -1):
                print("Found Existing Entry")
                item['Probes'] = pd.DataFrame(list(self.p_table.find(ProbeID=item['ProbeID'])))
                print("Got Probes")
                return item
        else:
            return None


    def write_db(self, probes):
        probes['ProbeID'] = getrandbits(50)
        probes['Date'] = arrow.now().datetime
        probes['Variants'] = self.variants
        probe_df = probes.pop('Probes')
        self.table.insert(probes)
        probe_df['ProbeID'] = probes['ProbeID']
        for k, probe in probe_df.T.to_dict().iteritems():
            self.p_table.insert(probe)

    def design(self, cds, min_probes):
        probes = []
        for gene, data in cds.iteritems():
            probe_set =  [p for p in [self.check_db(gene, mask) for mask in [5,4,3]]
                            if p is not None]
            if probe_set:
                print("Found {} probes for {}".format(len(probe_set), gene))
                for probe in probe_set:
                    probes.append(probe)
                continue
            for n, seq in enumerate(data['CDS List']):
                masking = 5
                key_box = {'ProbeSetName': gene,
                           "SpacingLength": "2",
                           "TargetSequence": seq
                            }
                selection = {'MaskingOrganism': self.organism,
                             "MaskingLevel": masking,
                             }
                if len(seq) < 20:
                    continue
                while True:
                    for k,v in key_box.iteritems():
                        elem = self.driver.find_element_by_name(k)
                        elem.clear()
                        elem.send_keys(v)
                    for k, v in selection.iteritems():
                        if k == 'MaskingLevel':
                            v = str(masking)
                        select = Select(self.driver.find_element_by_name(k))
                        select.select_by_value(v)
                    elem.send_keys(Keys.RETURN)
                    try:
                        element = WebDriverWait(self.driver, 60).until(
                                EC.presence_of_element_located((By.ID, "DesignResultWrapper"))
                            )
                    except:
                        #Incase designer times out
                        print("Error while waiting for response from server")
                        table = []
                        break     
                    table = pd.io.html.read_html(self.driver.page_source, header=0,)[1]
                    table.index = table['Probe #']
                    table = table.drop(table.columns[:2], axis=1)
                    probes.append({'Name': gene, 'Masking': masking,
                                   'Target Seq': seq, 'Probes': table,
                                   'CDS Region #': n, '# isoforms': data['# Isoforms'],
                                   })
                    self.write_db(probes[-1])
                    if masking <= 3:
                        break
                    masking -= 1
                    self.driver.back()
                    self.driver.refresh()
                print("Designed %i probes for %s cds#:%i with masking %i"
                      % (len(table), gene, n, masking))
                self.driver.back()
                self.driver.refresh()
        return probes


    def login(self):
        self.driver.get("https://www.biosearchtech.com/stellarisdesigner/default.aspx")
        elem = self.driver.find_element_by_name("loginuser")
        elem.send_keys('elubeck')
        elem = self.driver.find_element_by_name("loginpassword")
        elem.send_keys('1%$3UuVMvH%e')
        elem.send_keys(Keys.RETURN)

    def __init__(self, organism="mouse", variants=True):
        self.organism = organism
        self.variants = variants
        self.db = dataset.connect("sqlite:///db/biosearch.db")
        self.date = arrow.now()
        self.table = self.db[organism]
        self.p_table = self.db['probes']
        self.selenium = subprocess.Popen("java -jar ../lib/selenium-server-standalone-2.42.2.jar", shell=True) # Will probably only work on nix systems
        self.driver = webdriver.Firefox()
        self.driver.implicitly_wait(30) # seconds
        self.login()
