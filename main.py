from tdc.multi_pred import PPI
import numpy as np
import pandas as pd
import DataPreprocessor

from tdc.single_pred import ADME
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score
import pickle


def random_forest_classifier(data, target):
    folds = StratifiedKFold(n_splits=5)
    scores_rfc = []
    counter = 0
    for train_index, test_index in folds.split(np.array(data), np.array(target)):
        X_train, X_test, y_train, y_test = np.array(data)[train_index], np.array(data)[test_index], np.array(target)[train_index], np.array(target)[test_index]
        counter = counter + 1
        file_name1 = 'RandomForestClassifier_model_'
        file_name2 = 'RFC_predictions_'
        model = RandomForestClassifier(n_estimators=50)
        print("model:", model)
        model.fit(X_train, y_train)
        test_predictions_np = np.array(model.predict(X_test))
        scores_rfc.append(model.score(X_test, y_test))
        file_name2 = file_name2 + str(counter) + '.csv'
        # save the model to disk
        file_name1 = file_name1 + str(counter) + '.sav'
        pickle.dump(model, open(file_name1, 'wb'))
    print("scores_rfc:", scores_rfc)
    print("Average Score:", sum(scores_rfc) / len(scores_rfc))


def svc_classifier(data, target):
    folds = StratifiedKFold(n_splits=3, shuffle=True)
    scores_svc = []
    counter = 0
    for train_index, test_index in folds.split(np.array(data), np.array(target)):
        X_train, X_test, y_train, y_test = np.array(data)[train_index], np.array(data)[test_index], np.array(target)[train_index], np.array(target)[test_index]
        counter = counter + 1
        file_name1 = 'SVC_model_'

        model = SVC()
        print("model:", model)
        model.fit(X_train, y_train)
        test_predictions_np = np.array(model.predict(X_test))
        scores_svc.append(model.score(X_test, y_test))
        # file_name2 = file_name2 + str(counter) + '.csv'
        # np.savetxt(file_name2, test_predictions_np, fmt="%d")

        # save the model to disk
        file_name1 = file_name1 + str(counter) + '.sav'
        pickle.dump(model, open(file_name1, 'wb'))
        print("scores_svc:", scores_svc)
    print("scores_svc:", scores_svc)
    print("Average Score:", sum(scores_svc) / len(scores_svc))

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    preprocessor = DataPreprocessor.Preprocessor()
    # protein_pairs_df = preprocessor.get_protein_pair_sequences()
    # processed_protein_pairs_df = preprocessor.set_dataframe(protein_pairs_df)
    # processed_protein_pairs_df.to_csv('processed_protein_pairs_df.csv', header=True, index=False)
    # print(processed_protein_pairs_df.head())
    # df = pd.read_csv('protein_seq_list.csv')
    # preprocessor.set_dataframe(df).to_csv('processed_sequences.csv', header=True, index=False)
    # print(preprocessor.get_sequence_vector('AAAAIIIAA'))

    # X = pd.read_csv('processed_protein_pairs_df.csv')
    # y = preprocessor.get_protein_pair_target()
    # svc_classifier(X, y)
    # X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)
    # model = RandomForestClassifier(n_estimators=50)
    # model.fit(X_train, y_train)
    # print(model.score(X_test, y_test))
    # y_predicted = model.predict(X_test)

    # load the model from disk
    # loaded_model = pickle.load(open(filename, 'rb'))
    # result = loaded_model.score(X_test, Y_test)
    # print(result)

    # data = PPI(name='HuRI').neg_sample(frac=1)
    # df = data.get_data()
    # df.to_csv('all_df.csv', header=True, index=False)
    sample_df = pd.read_csv("all_df.csv")
    sample_df["Y"].to_csv('not_shuffled_targets.csv', header=True, index=False)
    print("helloooo")
