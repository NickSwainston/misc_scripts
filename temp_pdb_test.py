import os
import requests
from mwa_pulsar_client.client import *
from requests.auth import HTTPBasicAuth
from mwa_pulsar_client import client

SERVER = 'https://mwa-pawsey-volt01.pawsey.org.au'
#PASS = os.environ['PASS']
#AUTH = HTTPBasicAuth('dpallot', PASS)
#auth = HTTPBasicAuth('mwapulsar','veovys9OUTY=')
auth = HTTPBasicAuth(os.environ['MWA_PULSAR_DB_USER'],os.environ['MWA_PULSAR_DB_PASS'])

try:
    client.detection_create(SERVER,
                            auth,
                            observationid=1227009976,
                            pulsar='J0729-1448',
                            subband=145,
                            coherent=True,
                            observation_type=1,
                            startcchan=130,
                            stopcchan=140,
                            flux=80,
                            width=0.03,
                            dm=130)
except requests.exceptions.RequestException as e:
    import traceback
    traceback.print_exc()
