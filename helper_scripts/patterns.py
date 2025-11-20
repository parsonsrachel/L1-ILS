import pandas as pd
import warnings
warnings.filterwarnings('ignore')

def patterns(filename):
    df = pd.read_csv(filename)

    print(f'Starting number of TEs: {len(df["TE label"].unique())}')

    # remove rows that do not pass the orthologous check
    rm_te_lab = df.loc[df['orthologous'] == False]['TE label'].unique()
    for l in rm_te_lab:
        df.drop(df[df['TE label'] == l].index, inplace=True)

    print(f'Number of orthologous TEs: {len(df["TE label"].unique())}')
    print(f'Number of orthologous TE with presence/absence determined: {len(df["TE label"].unique()) - len(df.loc[df["presence/absence"] == "?"]["TE label"].unique())}')

    # remove rows that were not able to make a presence/absence call for
    for l in df.loc[df['presence/absence'] == '?']['TE label'].unique():
        df.drop(df[df['TE label'] == l].index, inplace=True)

    return df

def id_with_pat(df,pat):
    # get the TEs who have the pattern pat (SO, BO, Gorilla, Human, Chimp, Bonobo)
    # df[['TE label','species','presence/absence']]
    df['presence/absence'] = pd.to_numeric(df['presence/absence'], errors='coerce')

    pivot = df.pivot_table(index='TE label', columns='species', values='presence/absence', fill_value=0)

    desired = dict()
    ind = 0
    species = ["S.Orangutan", "B.Orangutan", "Gorilla", "Human", "Chimp", "Bonobo"]
    for i in pat:
        if i == 1:
            desired[species[ind]] = 1
        ind += 1
    del desired['Human']

    # All other species must be 0
    required_species = set(pivot.columns)
    desired_species = set(desired.keys())
    other_species = list(required_species - desired_species)

    if desired:
        present_spec = (pivot[list(desired.keys())[0]] == 1)
        for s in desired:
            present_spec &= (pivot[s] == 1)
    else:
        present_spec = pd.Series(True, index=pivot.index)

    matching_tes = pivot[
        present_spec &
        (pivot[other_species].sum(axis=1) == 0)  # All others must be zero
    ]
    return matching_tes.index.tolist()
