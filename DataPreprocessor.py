from tdc.multi_pred import PPI
import numpy as np
import pandas as pd
from itertools import islice


class Preprocessor:

    def __init__(self):
        self.df = None
        self.aa_group1 = ['A', 'G', 'V']
        self.aa_group2 = ['I', 'L', 'F', 'P', 'J']
        self.aa_group3 = ['Y', 'M', 'T', 'S']
        self.aa_group4 = ['H', 'N', 'Q', 'W']
        self.aa_group5 = ['R', 'K']
        self.aa_group6 = ['D', 'E']
        self.aa_group7 = ['C']
        self.any = ['X', 'U', 'O']
        all_aa_letters = ['A', 'G', 'V', 'I', 'L', 'F', 'P', 'Y', 'M', 'T', 'S', 'H', 'N', 'Q', 'W', 'R', 'K', 'D', 'E', 'C']
        aa_groups = ['g1', 'g2', 'g3', 'g4', 'g5', 'g6', 'g7', 'ANY']
        self.feature_list = []
        one_triple = ""
        for i in aa_groups:
            for j in aa_groups:
                for k in aa_groups:
                    one_triple = i + "_" + j + "_" + k
                    self.feature_list.append(one_triple)

    def get_triple_feature_type(self, triple):
        groups = []
        for aa in triple:
            if aa in self.aa_group1:
                groups.append("g1")
            if aa in self.aa_group2:
                groups.append("g2")
            if aa in self.aa_group3:
                groups.append("g3")
            if aa in self.aa_group4:
                groups.append("g4")
            if aa in self.aa_group5:
                groups.append("g5")
            if aa in self.aa_group6:
                groups.append("g6")
            if aa in self.aa_group7:
                groups.append("g7")
            if aa in self.any:
                groups.append("ANY")
        # print(triple, groups)
        return groups[0] + "_" + groups[1] + "_" + groups[2]

    @staticmethod
    def window(sequence, n):
        it = iter(sequence)
        result = tuple(islice(it, n))
        if len(result) == n:
            yield result
        for elem in it:
            result = result[1:] + (elem,)
            yield result

    def get_aa_triples(self, sequence):
        return ["".join(x) for x in self.window(sequence, 3)]

    def get_sequence_vector(self, sequence):
        seq_vector = [0] * len(self.feature_list)
        count = 0
        triples = self.get_aa_triples(sequence)
        for triple in triples:
            count = seq_vector[self.feature_list.index(self.get_triple_feature_type(triple))]
            seq_vector[self.feature_list.index(self.get_triple_feature_type(triple))] = count + 1

        for i in range(len(seq_vector)):
            seq_vector[i] = seq_vector[i] / len(triples)

        min_frequency = min(seq_vector)
        max_frequency = max(seq_vector)

        for i in range(len(seq_vector)):
            seq_vector[i] = (seq_vector[i] - min_frequency) / max_frequency

        return seq_vector

    def set_dataframe(self, df):
        df = df.replace('\*', '', regex=True)

        # Protein 1
        all_sequences = np.array(df['Protein1'])
        vector = []
        data = []
        for seq in all_sequences:
            vector = self.get_sequence_vector(seq)
            data.append((*vector,))
        protein1_features = ["P1_" + feature for feature in self.feature_list]
        print("protein1_features", protein1_features)
        new_df_1 = pd.DataFrame.from_records(data, columns=protein1_features)

        # Protein 2
        all_sequences = np.array(df['Protein2'])
        vector = []
        data = []
        for seq in all_sequences:
            vector = self.get_sequence_vector(seq)
            data.append((*vector,))
        protein2_features = ["P2_" + feature for feature in self.feature_list]
        print("protein2_features", protein2_features)
        new_df_2 = pd.DataFrame.from_records(data, columns=protein2_features)
        return pd.concat([new_df_1, new_df_2], axis=1)


    def get_protein_pair_sequences(self):
        data = PPI(name='HuRI').neg_sample(frac=1)
        self.df = data.get_data()
        # df = df.drop(['Protein1', 'Protein2'], axis=1)
        # stacked = df[['Protein1_ID', 'Protein2_ID']].stack()
        # df[['Protein1_ID', 'Protein2_ID']] = pd.Series(stacked.factorize()[0]+1, index=stacked.index).unstack()
        # DF1 = df.loc[df['Y'] == 1]
        # DF0 = df.loc[df['Y'] == 0]
        # DF1.to_csv('hi-iii.csv', header=False, index=False)
        # DF1_test = DF1.sample(frac=0.2)
        # DF1_train = DF1.drop(DF1_test.index).sample(frac=0.875)
        # DF1_val = DF1.drop(DF1_test.index).drop(DF1_train.index)
        # DF0_test = DF0.sample(frac=0.2)
        # DF0_train = DF0.drop(DF0_test.index).sample(frac=0.875)
        # DF0_val = DF0.drop(DF0_test.index).drop(DF0_train.index)
        # train = pd.concat([DF1_train, DF0_train], axis=0)
        # test = pd.concat([DF1_test, DF0_test], axis=0).sample(frac=1)
        # val = pd.concat([DF1_val, DF0_val], axis=0).sample(frac=1)
        # protein_df1 = df.drop(['Protein2', 'Protein2_ID', 'Y'], axis=1)
        # protein_df1 = protein_df1.rename(columns={"Protein1": "Protein", "Protein1_ID": "Protein_ID"})
        # protein_df2 = df.drop(['Protein1', 'Protein1_ID', 'Y'], axis=1)
        # protein_df2 = protein_df2.rename(columns={"Protein2": "Protein", "Protein2_ID": "Protein_ID"})
        # protein_df = pd.concat([protein_df1, protein_df2]).drop_duplicates()
        # protein_df.drop(['Protein_ID'], axis=1).to_csv('protein_seq_list.csv', header=True, index=False)
        # protein_df.drop(['Protein'], axis=1).to_csv('protein_id_list.csv', header=False, index=False)
        target = self.df['Y']
        return self.df.drop(['Protein1_ID', 'Protein2_ID', 'Y'], axis=1)

    def get_protein_pair_target(self):
        return self.df['Y']


