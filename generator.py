import os
import sys
import pandas as pd
import numpy as np

def trips_generator(I=200, t=1000, e=20, d=30):

    np.random.seed(0)

    trip_data = {
        'code': range(1, I + 1),
        'start': np.random.randint(0, t, I),
        'end': 0, 
        'dist': np.random.randint(5, 80, I),
        'de': np.random.randint(1, 15, I),
        'dd': np.random.randint(1, 15, I)
    }

    trip_data['start'] = np.sort(trip_data['start'])
    trip_data['end'] = trip_data['start'] + np.random.randint(10, 100, I)

    df1 = pd.DataFrame(trip_data)

    code = df1['code']
    start = df1['start']
    end = df1['end']

    pairs = []
    for i in range(len(code)):
        for j in range(len(code)):
            if end[i] < start[j] and code[i] != code[j]:
                pairs.append((code[i], code[j]))
    
    df2 = pd.DataFrame(pairs, columns=['i', 'j'])

    buses = []
    for n in range(e):
        buses.append((n+1, 'e'))
    for m in range(d):
        buses.append((m+n+2, 'd'))

    df3 = pd.DataFrame(buses, columns=['number', 'type'])

    if not os.path.exists(f'data/{I}_{t}'):
        os.makedirs(f'data/{I}_{t}')

    path1 = f'data/{I}_{t}/trips.csv'
    df1.to_csv(path1, sep='\t', index=False)

    path2 = f'data/{I}_{t}/pairs.csv'
    df2.to_csv(path2, sep='\t', index=False)

    path3 = f'data/{I}_{t}/buses.csv'
    df3.to_csv(path3, sep='\t', index=False)

if __name__ == "__main__":

    if len(sys.argv) > 1:
        N = int(sys.argv[1])
        T = int(sys.argv[2])
        E = int(sys.argv[3])
        D = int(sys.argv[4])
        trips_generator(N, T, E, D)

    else:
        trips_generator()