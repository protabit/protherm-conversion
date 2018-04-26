import requests
import sys
import json
import time
HOST="https://www.protabank.org/"
USERNAME=sys.argv[1]
PASSWORD=sys.argv[2]
r = requests.post(HOST + 'subscriber/api/obtain-token/', data={'username':USERNAME, 'password':PASSWORD})
token = r.json()['token']
headers = {'Authorization': 'Token {}'.format(token), 'content-type':'application/json'}
f=open('submitted_studies.txt','r')
ids = [int(line) for line in f.readlines()]
f.close()
for id in [1]:
    with open('protabank_studies/study'+str(id)+'.json') as json_data:
        d = json.load(json_data)
        json_data.close()
	response = requests.post(HOST + 'input/api/submit/', headers=headers, data=json.dumps(d))
    print "Submitted ", id, " with response ", response.status_code
    if response.status_code != 201:
        print response.content
    else:
        print "submitted", id
    time.sleep(1)
