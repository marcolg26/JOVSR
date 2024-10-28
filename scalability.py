import os
import time
from model import JOVSR

def scalability_test():

    print("*** Analisi scalabilità iniziata")

    for instance in os.listdir("data/scalability"):
        start = time.time()
        scheduling_model = JOVSR(f'scalability/{instance}')
        if scheduling_model.solve():
            print(f"Istanza {instance}: (2) {round(time.time() - start,2)} s")
        else:
            print(f"Istanza {instance}: (3) {round(time.time() - start,2)} s")

    print("*** Analisi scalabilità completata")

if __name__ == "__main__":
    scalability_test()
