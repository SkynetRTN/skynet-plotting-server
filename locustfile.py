from locust import FastHttpUser, task, tag, between, events, constant_throughput, constant, log
from factory import faker, DictFactory, fuzzy
from json import dumps
import base64
import sys


class ParamFactory(DictFactory):
    ratioMass = fuzzy.FuzzyFloat(1,10)
    totalMass = fuzzy.FuzzyFloat(2.5,250)
    ratioMassStrain = 1
    totalMassStrain = 1
    phaseStrain = 1
    
gravfilejson = {"default_set": False,
               }

class GravityUser(FastHttpUser):

    # @events.test_start.add_listener
    # def on_test_start(environment, **kwargs):
    #     with open("./gravity-model-data/default-set/GW150914.hdf5", 'rb') as file:
    #         with environment.client.post("/gravprofile", files={'file':file.read()}, data={"default_set": 'false'}, catch_response=True) as response:
                
    #             if response.ok:
    #                 data = response.json()
    #                 with open("locustimage.png", "wb") as fh:
    #                     fh.write(base64.decodebytes(data.image))

# TEST
    @tag("gravityfile")
    @task(1)
    def gravfile(self):
        with open("./gravity-model-data/default-set/GW150914.hdf5", 'rb') as file:
            with self.client.post("/gravfile", files={'file':file.read()}, data={"default_set": 'false','filename':"GW150914.hdf5"}, catch_response=True) as response:
                
                if not response.ok:
                    return
                sessionID = response.json()["sessionID"]

                with self.client.post("/gravitydata", totalMass=1, ratioMass=1, sessionID=sessionID) as retrieved_data:
                    js = retrieved_data.json()
                    
                    if not retrieved_data.ok:
                        print(js["err"])
                        

    @tag("gravityprofile")
    @task(1)
    def gravprofile(self):
        with open("./gravity-model-data/default-set/GW150914.hdf5", 'rb') as file:
            with self.rest("POST", "/gravprofile", files={'file':file.read()}, data={"default_set": 'false','filename':"GW150914.hdf5"}) as response:
                
                if response.ok:
                    # data = response.js
                    # img = base64.b64decode(bytes(data["image"], "utf-8"))
                    # with open("zClientImage.png", "wb", buffering=0) as fh:
                    #     print(fh.write(img))
                    # print("image saved")
                    pass

    host = "http://127.0.0.1:8000/"
    wait_time = constant(10)