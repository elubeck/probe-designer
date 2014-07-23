from __future__ import print_function, with_statement, division
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.ui import Select
import pandas as pd
import subprocess


class Biosearch(object):
    def close(self):
        self.driver.close()
        self.selenium.close()
    
    def design(self, cds, min_probes):
        probes = []
        for gene, data in cds.iteritems():
            for n, seq in enumerate(data['CDS List']):
                masking = 5
                key_box = {'ProbeSetName': gene,
                           "SpacingLength": "2",
                           "TargetSequence": seq
                            }
                selection = {'MaskingOrganism': "mouse",
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
                    if len(table) >= min_probes or masking <= 3:
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

    def __init__(self):
        self.selenium = subprocess.Popen("java -jar ../lib/selenium-server-standalone-2.42.2.jar", shell=True) # Will probably only work on nix systems
        self.driver = webdriver.Firefox()
        self.driver.implicitly_wait(30) # seconds
        self.login()
