from model import JOVSR

def model_test(instance = "demo"):

    parameters = {
    "Q": 400,           # battery capacity [kWh]
    "c1": 3.9,          # electric consumption rate
    "C": 260,           # diesel capacity [L]
    "c2": 1.3,          # diesel consumption rate [L/km]
    "v": 0.41,          # speed of buses [km/min]
    "r1": 2,            # recharging rate of electric buses
    "r2": 30,           # recharging rate of diesel buses
    "u1": 1,            # time interval of each charging activity
    "u2": 1,            # time interval of each replenishing activity
    "alpha0": 1,        # fixed cost of performing a charge
    "alpha1": 0.35,     # unit electricity cost [€/kWh]
    "alpha2": 1,        # fixed cost of performing a diesel replenishment
    "alpha3": 1.6,      # unit diesel cost [€/L]
    "Qmin": 40,         # minimum charge limit
    "Cmin": 20,         # minimum diesel limit
    "chargers": 7       # number of fast chargers 
    }

    scheduling_model = JOVSR(instance, parameters)

    if scheduling_model.solve():
        scheduling_model.print_results()
        scheduling_model.plot_results()
    else:
        print("Modello non risolvibile")

if __name__ == "__main__":
    model_test()
