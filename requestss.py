import requests
from main import main

genome_sequencer = "https://genome-sequencer.herokuapp.com/seq"

r = requests.get(genome_sequencer)

data = r.json()["data"]

def alex_genomes():
    for element in data:
        print(element["attributes"]["species"])
        genome_seq = element["attributes"]["sequence"]
        main(genome_seq)

alex_genomes()
