from __future__ import division, print_function, with_statement

import subprocess
from random import getrandbits

import arrow
import dataset
import pandas as pd
from future.builtins import object, str
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.ui import Select, WebDriverWait


class Biosearch(object):
    def close(self):
        if self.is_open:
            self.driver.close()
            self.selenium.kill()

    def check_db(self, gene, masking):
        res = self.table.find(Name=gene, Variants=self.variants, Masking=masking)
        passed = []
        for item in res:
            if arrow.get(item['Date']) >= self.date.replace(years = -1):
                item['Probes'] = pd.DataFrame(list(self.p_table.find(ProbeID=item['ProbeID'])))
                passed.append(item)
        return passed


    def write_db(self, probes):
        probes['ProbeID'] = getrandbits(50)
        probes['Date'] = arrow.now().datetime
        probes['Variants'] = self.variants
        probe_df = probes.pop('Probes')
        self.table.insert(probes)
        probe_df['ProbeID'] = probes['ProbeID']
        for k, probe in probe_df.T.to_dict().items():
            self.p_table.insert(probe)

    def design(self, cds, min_probes):
        probes = []

        # Check if biosearch can design probes for organism
        organism_options = ["human", "mouse", "rat", "drosophila", "celegans", "other"]
        if self.organism in organism_options:
            chosen_org = self.organism
            mask_range = [5, 4, 3]
        else:
            chosen_org = "other"
            mask_range = [2]

        for gene, data in cds.items():
            probe_set =  [p for p in [self.check_db(gene, mask) for mask in mask_range]
                            if p]
            if probe_set:
                for probe in probe_set:
                    probes.append(probe)
                pl = [pd.concat([p['Probes'] for p in probe], ignore_index=True).shape[0] for probe in probe_set]
                print("Found unblasted probeset for {} with {} probes".format(gene,pl))
                continue
            if self.is_open is False:
                self.open()
            cds_list = []
            for seq in data['CDS List']:
                if len(seq) < 8000:
                    cds_list.append(seq)
                elif len(seq) >= 8000:
                    n_chunks = 1 + len(seq) // 8000
                    for i in range(n_chunks):
                        cds_list.append(seq[i*8000:(i+1)*8000])
            print(len(cds_list), len(data['CDS List']))
            for n, seq in enumerate(cds_list):
                key_box = {'ProbeSetName': gene,
                           "SpacingLength": "2",
                           "ProbesNumber": "200",
                           "TargetSequence": seq
                            }
                selection = {'MaskingOrganism': chosen_org,
                             "MaskingLevel": mask_range[0],
                             }
                if len(seq) < 20:
                    continue
                for masking in mask_range:
                    for k,v in key_box.items():
                        elem = self.driver.find_element_by_name(k)
                        elem.clear()
                        elem.send_keys(v)
                    for k, v in selection.items():
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
                    self.driver.back()
                    self.driver.refresh()
                    print("Designed %i probes for %s cds#:%i with masking %i"
                          % (len(table), gene, n, masking))
        return probes

    def open(self):
        self.selenium = subprocess.Popen("java -jar ../lib/selenium-server-standalone-2.45.0.jar", shell=True) # Will probably only work on nix systems
        self.driver = webdriver.Firefox()
        self.driver.implicitly_wait(30) # seconds
        self.login()
        self.is_open = True

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
        self.is_open = False
