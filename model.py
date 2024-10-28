import gurobipy as gb
import pandas as pd
import sys
import matplotlib.pyplot as plt
import os


class dualLogger:  # per stampare sia in console che su file
    def __init__(self, filename):
        self.terminal = sys.stdout
        self.log = open(filename, 'w')

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        pass


# Joint Optimal Vehicle Scheduling and Recharging
class JOVSR:

    def __init__(self, folder="demo", parameters=None):

        self.folder = folder

        default = {
            "Q": 400,  # battery capacity [kWh]
            "c1": 0.9,  # electric consumption rate
            "C": 190,  # diesel capacity [L]
            "c2": 1.3,  # diesel consumption rate [L/km]
            "v": 10000,  # speed of buses [km/h]
            "r1": 8,  # recharging rate of electric buses
            "r2": 10,  # recharging rate of diesel buses
            "u1": 1,  # time interval of each charging activity
            "u2": 1,  # time interval of each replenishing activity
            "alpha0": 1,  # fixed cost of performing a charge
            "alpha1": 0.7,  # unit electricity cost [€/kWh]
            "alpha2": 1,  # fixed cost of performing a diesel replenishment
            "alpha3": 0.3,  # unit diesel cost [€/L]
            "Qmin": 40,  # minimum charge limit
            "Cmin": 60,  # minimum diesel limit
            "chargers" : 12 # number of fast chargers 
        }

        if parameters is not None:
            default.update(parameters)

        self.parameters = default

        self.model = gb.Model("bus")
        self.model.setParam('outputFlag', 0)

        self._load_data()
        self._initialize_model()

    def _load_data(self):

        self.buses = pd.read_csv(f"data/{self.folder}/buses.csv", sep='\t')
        self.trip_data = pd.read_csv(f"data/{self.folder}/trips.csv", sep='\t')
        self.pairs_data = pd.read_csv(
            f"data/{self.folder}/pairs.csv", sep='\t')
        self.trips = list(self.trip_data['code'])

    def _initialize_model(self):

        # SETS

        # set of all trips
        self.I = set(self.trips)

        # set of all eletric buses
        self.E = list(self.buses[self.buses['type'] == 'e']['number'])

        # set of all diesel buses
        self.D = list(self.buses[self.buses['type'] == 'd']['number'])

        # set of all fast chargers
        self.K = range(self.parameters["chargers"])

        # set of all scheduled trips and the trip 0 for all buses to depart from the depot to transit centers to start providing services
        self.I0 = set(self.trips)
        self.I0.add(0)

        self.N = max(self.I0)

        # set of all scheduled trips and the trip N+1 for all buses to return from transit centers to the depot after completing all services
        self.IN1 = set(self.trips)
        self.IN1.add(self.N + 1)

        # I∪{0, N+1}
        self.I0N1 = set(self.trips)
        self.I0N1.add(0)
        self.I0N1.add(self.N + 1)

        # set of all trips pairs
        self.P = set()
        self.P = set(zip(self.pairs_data['i'], self.pairs_data['j']))

        # set of all trips pairs including trips 0 and N+1
        self.Pstar = set(self.P)
        for i in self.I:
            self.Pstar.add((0, i))
            self.Pstar.add((i, self.N + 1))

        # distance of trip i
        self.di = dict(zip(self.trip_data['code'], self.trip_data['dist']))
        self.di[self.N + 1] = 15

        # start time of trip i
        self.ai = dict(zip(self.trip_data['code'], self.trip_data['start']))

        # end time of trip i
        self.bi = dict(zip(self.trip_data['code'], self.trip_data['end']))
        self.bi[0] = 1
        self.max_time = max(self.bi.values())
        self.ai[self.N + 1] = self.max_time + 5
        self.bi[self.N + 1] = self.max_time + 10

        # set of all time slots
        self.T = list(range(self.max_time + 20))

        # distance from the transit center of trip i to the depot to charge electric buses
        self.di_dc = dict(zip(self.trip_data['code'], self.trip_data['de']))
        self.di_dc[0] = 2
        self.di_dc[self.N + 1] = 2

        # distance from the transit center of trip i to the depot to refuel diesel buses
        self.di_dd = dict(zip(self.trip_data['code'], self.trip_data['dd']))
        self.di_dd[0] = 2
        self.di_dd[self.N + 1] = 2

        # set of all feasible charging activities
        self.A = set()
        for (i, j) in self.Pstar:
            for t in self.T:
                if self.bi[i] <= t <= self.ai[j]:
                    if (j < self.N + 1):
                        self.A.add((i, j, t))

        # PARAMETERS

        Q = self.parameters["Q"]  # battery capacity [kWh]
        c1 = self.parameters["c1"]  # eletric consumption rate
        C = self.parameters["C"]  # diesel capacity [L]
        c2 = self.parameters["c2"]  # diesel consumption rate [L/km]

        v = self.parameters["v"]  # speed of buses [km/h]
        r1 = self.parameters["r1"]  # recharging rate of eletric buses
        r2 = self.parameters["r2"]  # recharging rate of diesel buses
        tau = 1  # time duration of each time slot

        u1 = self.parameters["u1"]  # time interval of each charging activity
        u2 = self.parameters["u2"]  # time interval of each replenishing activity
        alpha0 = self.parameters["alpha0"]  # fixed cost of performing a charge
        alpha1 = self.parameters["alpha1"]  # unit electricity cost [€/kWh]
        alpha2 = self.parameters["alpha2"]  # fixed cost of performing a diesel replenishment
        alpha3 = self.parameters["alpha3"]  # unit diesel cost [€/L]

        M = 2000  # infinite costant
        Qmin = self.parameters["Qmin"]  # minimum charge limit
        Cmin = self.parameters["Cmin"]  # minimum diesel limit

        self.n = range(len(self.E))
        self.m = range(len(self.D))
        self.o = range(len(self.I))
        self.o0 = range(len(self.I0))
        self.on1 = range(len(self.IN1))
        self.o0n1 = range(len(self.I0N1))
        self.p = range(len(self.P))
        self.q = range(len(self.K))
        self.r = range(len(self.T))

        # DECISION VARIABLES

        # if the electric bus e ∈ E is scheduled, ze = 1
        self.ze = self.model.addVars(
            [e for e in self.n], vtype=gb.GRB.BINARY, name="ze")

        # if the electric bus d ∈ D is scheduled, zd = 1
        self.zd = self.model.addVars(
            [d for d in self.m], vtype=gb.GRB.BINARY, name="zd")

        # remaining charge of electric bus e ∈ E after completing trip i ∈ I
        self.Qie = self.model.addVars(
            [(i, e) for i in self.I0N1 for e in self.n], lb=0.0, name="Qie")

        # remaining diesel of diesel bus d ∈ D after completing trip i ∈ I
        self.Qid = self.model.addVars(
            [(i, d) for i in self.I0N1 for d in self.m], lb=0.0, name="Qid")

        # If the electric bus e ∈ E goes to the depot to replenish charge after completing the service of trip i ∈ I
        self.Fie = self.model.addVars(
            [(i, e) for i in self.I0N1 for e in self.n], vtype=gb.GRB.BINARY, name="Fie")

        # If the diesel bus d ∈ D goes to the depot to replenish diesel after completing the service of trip i ∈ I
        self.Fid = self.model.addVars(
            [(i, d) for i in self.I0N1 for d in self.m], vtype=gb.GRB.BINARY, name="Fid")

        # If the bus e ∈ E does choose charger k ∈ K to start replenishing charge at time slot t ∈ T after completing the service of trip i ∈ I, Z = 1
        self.Zit_ek = self.model.addVars(
            [(i, e, t, k) for i in self.I0N1 for e in self.n for t in self.T for k in self.q], vtype=gb.GRB.BINARY,
            name="Zitek")

        # If the trip pair (i, j) ∈ P is served sequentially by electric bus e ∈ E, then xi,j,e = 1
        self.xije = self.model.addVars([(i, j, e) for (
            i, j) in self.Pstar for e in self.n], vtype=gb.GRB.BINARY, name="xije")

        # If the trip pair (i, j) ∈ P is served sequentially by diesel bus d ∈ D, then xi,j,d = 1
        self.xijd = self.model.addVars([(i, j, d) for (
            i, j) in self.Pstar for d in self.m], vtype=gb.GRB.BINARY, name="xijd")

        self.xie = self.model.addVars(
            [(i, e) for i in self.I0N1 for e in self.n], vtype=gb.GRB.BINARY, name="xie")

        self.xid = self.model.addVars(
            [(i, d) for i in self.I0N1 for d in self.m], vtype=gb.GRB.BINARY, name="xid")
        
        self.y = self.model.addVars(self.I, vtype=gb.GRB.INTEGER, name="num_bus_per_corsa")

        self.model.ModelSense = gb.GRB.MINIMIZE
        self.model.setObjective(
             gb.quicksum((alpha0 + alpha1 * r1 * u1) * self.Fie[i, e] for i in self.I0N1 for e in self.n) +
             gb.quicksum((alpha2 + alpha3 * r2 * u2) * self.Fid[i, d] for i in self.I0N1 for d in self.m)
        )

        # CONSTRAINTS

        # [2]
        for i in self.I:
            self.model.addConstr(
                gb.quicksum(self.xije[i, j, e] for j in self.IN1 for e in self.n if (i, j) in self.Pstar) +
                gb.quicksum(self.xijd[i, j, d] for j in self.IN1 for d in self.m if (i, j) in self.Pstar) == 1,
                name = f"#2 {i}"
            )
        
        for i in self.I:
            self.model.addConstr(
                self.y[i] == gb.quicksum(self.xije[i, j, e] for j in self.IN1 if (i, j) in self.Pstar for e in self.n) +
               gb.quicksum(self.xije[j, i, e] for j in self.IN1 if (j, i) in self.Pstar for e in self.n) +
               gb.quicksum(self.xijd[i, j, d] for j in self.IN1 if (i, j) in self.Pstar for d in self.m) +
               gb.quicksum(self.xijd[j, i, d] for j in self.IN1 if (j, i) in self.Pstar for d in self.m),
               name=f"bus-trip count_{i}"
               )
        
        for i in self.I:
            self.model.addConstr(
                self.y[i] <=2
                )

        # [3]
        for e in self.n:
            self.model.addConstr(
                gb.quicksum(self.xije[i, j, e] for (i, j) in self.Pstar) <= self.ze[e] * M,
                name=f"#3 {e}"
            )

        # [4]
        for d in self.m:
            self.model.addConstr(
                gb.quicksum(self.xijd[i, j, d] for (i, j) in self.Pstar) <= self.zd[d] * M,
                name=f"#4 {d}"
            )

        # [5]
        self.model.addConstr(gb.quicksum(self.ze[i] for i in self.n) <= len(self.E), name="#5")

        # [6]
        self.model.addConstr(gb.quicksum(self.zd[i] for i in self.m) <= len(self.D), name="#6")

        # [7]
        for e in self.n:
            for i in self.I:
                self.model.addConstr(
                    self.Fie[i, e] <= gb.quicksum(self.xije[i, j, e] for j in self.IN1 if (i, j) in self.Pstar),
                    name=f"#7 {e} {i}"
                )

        # [8]
        for d in self.m:
            for i in self.I:
                self.model.addConstr(
                    self.Fid[i, d] <= gb.quicksum(self.xijd[i, j, d] for j in self.IN1 if (i, j) in self.Pstar),
                    name=f"#8 {d} {i}"
                    )

        # [9]
        for e in self.n:
            for (i, j) in self.Pstar:
                self.model.addConstr(
                    self.xije[i, j, e] * Qmin <= self.Qie[j, e],
                    name=f"#9L ({i},{j}), {e}"
                )
                
                self.model.addConstr(
                    self.Qie[j, e] <= 
                    self.Qie[i, e] + (r1 * u1 - self.di_dc[i] * c1) * self.Fie[i, e] - c1 * self.di[j] *
                    (self.xije[i, j, e]) + Q * (1 - self.xije[i, j, e]),
                    name=f"#9R ({i},{j}), {e}"
                )

        # [10]
        for d in self.m:
            for (i, j) in self.Pstar:
                self.model.addConstr(
                    self.xijd[i, j, d] * Cmin <= self.Qid[j, d],
                    name=f"#10L ({i},{j}), {d}"
                )
                
                self.model.addConstr(
                    self.Qid[j, d] <=
                    self.Qid[i, d] + (r2 * u2 - self.di_dd[i] * c2) * self.Fid[i, d] -c2 * self.di[j] * 
                    self.xijd[i, j, d] + C * (1 - self.xijd[i, j, d]),
                    name=f"#10R ({i},{j}), {d}"
                )

        # [11]
        for i in self.I:
            for e in self.n:
                self.model.addConstr(
                    self.Qie[i, e] + (r1 * u1 - self.di_dc[i] * c1) * self.Fie[i, e] <= Q,
                    name=f"#11 {e} {i}"
                )

        # [12]
        for i in self.I:
            for d in self.m:
                self.model.addConstr(
                    self.Qid[i, d] + (r2 * u2 - self.di_dd[i] * c2) * self.Fid[i, d] <= C,
                    name=f"#12 {e} {i}"
                )

        # [13]
        for i1 in self.I:
            for e in self.n:
                self.model.addConstr(
                    gb.quicksum(self.Zit_ek[i, e, t, k] for k in self.q for (i, j, t) in self.A if i==i1) == self.Fie[i1, e],
                    name=f"#13 {e} {i}"
                )

        # [14]
        for k in self.q:
            self.model.addConstr(
                gb.quicksum(self.Zit_ek[i, e, t1, k] for e in self.n for (i, j, t1) in self.A for t in self.T if t - u1 + 1 <= t1 <= t) <= 1,
                name=f"#14 {k}"
            )

        # [15]
        for d in self.m:
            for (i, j) in self.Pstar:
                self.model.addConstr(
                    self.Fid[i, d] * (self.bi[i] + 2 * self.di_dd[i] / v + u2) <= self.xijd[i, j, d] * self.ai[j] + 10000 * (1 - self.xijd[i, j, d]),
                    name=f"#15 ({i},{j}), {d}"
            )
                
        limit_dep = round((len(self.I0N1))/10, 0)
        if (limit_dep==0):
            limit_dep=1
        
        for e in self.n:
            self.model.addConstr(
                gb.quicksum(self.xije[i, self.N + 1, e] for i in self.I if (i, self.N + 1) in self.Pstar) <= limit_dep)  # dep max e

        for d in self.m:
            self.model.addConstr(
                gb.quicksum(self.xijd[i, self.N + 1, d] for i in self.I if (i, self.N + 1) in self.Pstar) <= limit_dep)  # dep max d

        for i in self.I:
            for e in self.n:
                self.model.addConstr(
                    self.xie[i, e] == gb.quicksum(self.xije[i, j, e] for j in self.IN1 if (i, j) in self.Pstar) +
                    gb.quicksum(self.xije[j, i, e] for j in self.IN1 if (j, i) in self.Pstar)
                )

        for i in self.I:
            for d in self.m:
                self.model.addConstr(
                    self.xid[i, d] == gb.quicksum(self.xijd[i, j, d] for j in self.IN1 if (i, j) in self.Pstar) +
                    gb.quicksum(self.xijd[j, i, d] for j in self.IN1 if (j, i) in self.Pstar)
                )

    def solve(self):
        self.model.optimize()
        return self.model.Status == gb.GRB.OPTIMAL

    def print_results(self):

        if not os.path.exists(f'results/{self.folder}'):
            os.makedirs(f'results/{self.folder}')

        sys.stdout = dualLogger(f'results/{self.folder}/assignments.txt')

        self.fig, self.ax = plt.subplots(figsize=(10, 6))

        bus_schedules = {e: [] for e in self.E}
        bus_schedules.update({d: [] for d in self.D})

        for (i, j) in self.Pstar:
            for e in self.n:
                if (self.xije[i, j, e].x > 0.5):
                    print("Corse (" + str(i) + "," + str(j) + ") effettuate da bus D #" + str(self.E[e]) + " [" + str(self.ai[i]) + "-" + str(self.bi[j]) + "]")
                    if (j == self.N + 1):
                        bus_schedules[self.E[e]].append((self.ai[i], self.bi[i], 'limegreen', str(i)))
                    else:
                        bus_schedules[self.E[e]].append((self.ai[i], self.bi[i], 'limegreen', str(i)))
                        bus_schedules[self.E[e]].append((self.bi[i], self.ai[j], 'gray', "x"))
                        bus_schedules[self.E[e]].append((self.ai[j], self.bi[j], 'limegreen', str(j)))

        for (i, j) in self.Pstar:
            for d in self.m:
                if (self.xijd[i, j, d].x > 0.5):
                    print("Corse (" + str(i) + "," + str(j) + ") effettuate da bus E #" + str(self.D[d]) + " [" + str(self.ai[i]) + "-" + str(self.bi[j]) + "]")
                    if (j == self.N + 1):
                        bus_schedules[self.D[d]].append((self.ai[i], self.bi[i], 'orange', str(i)))
                    else:
                        bus_schedules[self.D[d]].append((self.ai[i], self.bi[i], 'orange', str(i)))
                        bus_schedules[self.D[d]].append((self.bi[i], self.ai[j], 'gray', "x"))
                        bus_schedules[self.D[d]].append((self.ai[j], self.bi[j], 'orange', str(j)))

        for d in self.m:
            for i in self.I0N1:
                if (self.Fid[i, d].x > 0):
                    print(f"Bus {self.D[d]} effettua rifornimento dopo la corsa {i}")
                    if (i != self.N + 1 and i > 0):
                        self.ax.scatter(self.bi[i] + 1, self.D[d], color='yellow', s=100, zorder=2, marker='^')

        for e in self.n:
            for i in self.I0N1:
                for t in self.T:
                    for k in self.q:
                        if (self.Zit_ek[i, e, t, k].x > 0):
                            print(f"Bus {self.E[e]} in ricarica dopo corsa {i} (charger {k}, minuto {t})")
                            self.ax.scatter(t + 1, self.E[e], color='lightskyblue', s=100, zorder=2, marker=f"${k}$")

        for bus_id, schedule in bus_schedules.items():
            for start, end, color, text in schedule:
                self.ax.barh(bus_id, end - start, left=start,height=0.4, align='center', color=color)
                self.ax.text((start + end) / 2, bus_id, text, color='black', va='center', ha='center')

        self.ax.set_xlabel('Orario')
        self.ax.set_ylabel('Macchina')
        self.ax.set_title('Assegnazioni autobus')
        self.ax.grid(True)

    def plot_results(self):
        plt.savefig(f'results/{self.folder}/assignments.png')
        plt.show()
