import selenium,sys
from selenium import webdriver


def extract(miRNA_name):
    options = webdriver.ChromeOptions() 
    dir='/Users/quyixiang/Documents/result/'+miRNA_name
    dir=str(dir)
    prefs = {'profile.default_content_settings.popups': 0, 'download.default_directory': dir }
    options.add_experimental_option('prefs', prefs) 

    driver = webdriver.Chrome(chrome_options=options)
    driver.get('http://ophid.utoronto.ca/mirDIP/index.jsp')

    textarea_miRNA = driver.find_element_by_id('idMicroRna')
    textarea_miRNA.send_keys(miRNA_name)

    submit_button = driver.find_element_by_id('submitButton')
    submit_button.click()

    csv_radio = driver.find_element_by_id('csv')
    csv_radio.click()

    download_button = driver.find_element_by_id('submitButton')
    download_button.click()

    link = driver.find_element_by_id('idm_Index')
    link.click()


extract(sys.argv[1])
