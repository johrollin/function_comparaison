# /usr/bin/python3

"""
This script compare annotation method results  

* 1: Removing stop word with nltk  
* 2: For each gene compare each annotation result and if ratio > treeshold keep the best annotation  
"""
import argparse
import yaml
from nltk.corpus import stopwords
from nltk.tokenize import word_tokenize
from fuzzywuzzy import fuzz
import pandas as pd
from nltk.tokenize import RegexpTokenizer
import operator
import Levenshtein # make fuzzy faster to use levenshtein than difflib

def annotation_launch():
    """ 
    Main function of the script

    """

    parser = argparse.ArgumentParser()
    parser.add_argument('-yaml', help="yaml containing all needed parameter", required=True)
    args = parser.parse_args()
    yaml_file = args.yaml
    try:
        stream = open(yaml_file, 'r')
    except OSError as e:
        raise('cannot open yaml file')

    # extract all data from yaml parameter file
    df_all_data, output_file, column_to_compare, list_column_to_filter, \
    list_value_of_filter, how_to_compare, list_column_to_keep, threshold_ratio, \
    column_to_group = clean_args_data(stream)    

    # clean column_to_compare before the comparaison and filter the data
    df_clean_filter, df_filter = clean_data(\
        df_all_data, list_column_to_keep, column_to_compare,\
        list_column_to_filter, list_value_of_filter, how_to_compare)

    # compare text of column_to_compare
    dict_annot_cpd = compare_annotation(df_clean_filter, \
    threshold_ratio, column_to_group, column_to_compare, list_column_to_filter)


    print(df_clean_filter)

def clean_args_data(stream):
    """ 
    Read yaml file to extract arguments to use  
    Open data_file and extract data into a panda dataframe  

    :type stream: file (from yaml)  
    :param stream: parameter file  

    """
    # read yaml file
    try:
        data = yaml.load(stream)
    except:
        raise('yaml file format is not valid')
    data_file_path = data['data_file']
    # read data file put the data in data_frame
    try:
        data_file = open(data_file_path, 'r')
        df_all_data = pd.read_csv(data_file, sep='\t', header=0, dtype='str')
    except OSError as e:
        msg = 'cannot open data_file: ' + str(data_file)
        raise(msg)
    data_file.close()

    output_file = data['output_file']
    column_to_group = data['column_to_group']
    column_to_compare = data['column_to_compare']
    list_column_to_filter = data['list_column_to_filter']
    list_value_of_filter = data['list_value_of_filter']
    how_to_compare = data['how_to_compare']
    list_column_to_keep = data['list_column_to_keep']
    threshold_ratio = data['threshold_ratio']

    return df_all_data, output_file, column_to_compare, list_column_to_filter, \
        list_value_of_filter, how_to_compare, list_column_to_keep, threshold_ratio, \
        column_to_group

def make_comparaison(comparator, data_value, filter_value):
    """
    make filter comparaison to choose if the data must be keeped or not  

    :type comparator: string  
    :param comparator: how_to_compare[]  

    :type data_value: float  
    :param data_value:   

    :type filter_value: float  
    :param filter_value: list_value_of_filter[]  

    """
    is_valid = False

    if comparator=='>=':
        if data_value >= filter_value:
            is_valid = True
    elif comparator=='>':
        if data_value > filter_value:
            is_valid = True
    elif comparator=='<=':
        if data_value <= filter_value:
            is_valid = True
    elif comparator=='<':
        if data_value < filter_value:
            is_valid = True
    elif comparator=='=':
        if data_value >= filter_value:
            is_valid = True
    return is_valid
    
def clean_data(df, list_column_to_keep, column_to_compare, list_column_to_filter, list_value_of_filter, how_to_compare):
    """
    Return two DataFrame with only the column to keep that pass the filter values  
    The punctuation and stop world will be removed from the column to compare  

    :param df: DataFrame to clean  
    :type df: DataFrame  
    :param list_column_to_keep: list of column name to keep  
    :type list_column_to_keep: list ['str']  
    :param column_to_compare: column name to clean  
    :type column_to_compare: str  
    :param list_column_to_filter: list of column name to filter  
    :type list_column_to_filter: list ['str']  
    :param list_value_of_filter: list of filter values treeshold data >= values  
    :type list_value_of_filter: list ['float']  
    :param how_to_compare: (>, >=, <, <=, =)  
    :type how_to_compare: list ['']  
    """

    # Remove NULL values
    df = df.dropna(subset=list_column_to_keep)
    # reset row index
    df.reset_index(inplace=True)

    # create dataframe with only considered column
    new_df = pd.DataFrame()
    df_filter = pd.DataFrame()
    df_clean_filter = pd.DataFrame()

    for i in range(0, len(list_column_to_keep)):
        col = list_column_to_keep[i]
        new_df.insert(i, col, df[col], True)
        df_clean_filter[col] = ""
        df_filter[col] = ""

    # list of stop world
    stop_words = set(stopwords.words('english'))
    # ponctuation
    tokenizer = RegexpTokenizer(r'\w+')

    count = 0
    # enum on dataframe
    for i, row in new_df.iterrows():
        new_desc = ''
        result = []
        is_valid = True
        # make filter using list_column_to_filter 
        if list_column_to_filter != []:
            for j in range(0, len(list_column_to_filter)):
                data_value = float(row[list_column_to_filter[j]])
                filter_value = float(list_value_of_filter[j])
                comparator = how_to_compare[j]
                is_valid = make_comparaison(comparator, data_value, filter_value)
        # data value pass the given filter requirement 
        if is_valid:
            # fill df filtrer with column_to_compare data
            df_filter = df_filter.append(row)
            count += 1
            # select stop word + ponctuation of column_to_compare
            word_tokens = word_tokenize(row[column_to_compare])
            for word in word_tokens:
                # keep only meaning word
                if not word in stop_words:
                    new_desc += word + ' '
            # remove ponctuation
            result = tokenizer.tokenize(new_desc)
            # put the clean column_to_compare in df
            row[column_to_compare] = ''
            for element in result:
                row[column_to_compare] += element + ' '
            # remplis le df clean_filtrer avec column_to_compare nettoye
            df_clean_filter = df_clean_filter.append(row)

    # both index must be the same
    # df_clean_filter will be use to make fuzzy comparaison
    # df_filter will be use to make the ouptpu file
    df_filter.reset_index(inplace=True, drop=True)
    df_clean_filter.reset_index(inplace=True, drop=True)
    return df_clean_filter, df_filter


def compare_annotation(df, threshold_ratio, column_to_group, column_to_compare, list_column_to_filter): 
    # TODO prendre result avec + de hit > ratio et pas le meilleur ratio
    # TODO use list_column_to_filter to choos better result
    """ 
    For each gene compare each annotation result and if ratio > treeshold keep the best annotation  

    :type df: 
    :param df:

    :type threshold_ratio:
    :param threshold_ratio:

    :type column_to_group:
    :param column_to_group: 

    :type column_to_compare:
    :param column_to_compare: 

    :type list_column_to_filter:
    :param list_column_to_filter:     
    """


    dict_annot_cpd = {}
    for i, row in df.iterrows():
        current_GO_id = row[column_to_group]
        # On fait des sous groupe de resultat pour le meme GO_id
        subgroup = df.loc[df[column_to_group] == current_GO_id]
        # si le Go_id a plusieurs resultat est n'a pas deja ete traite
        if len(subgroup) > 1 and current_GO_id not in dict_annot_cpd.keys():
            dict_annot_cpd[current_GO_id] = ""
            dict_tmp_annot = {}
            # calcul l'index sur dataframe du sous groupe de resultat sur le meme gene
            max_index = i+len(subgroup)-1
            # on itere sur les resulats pour le meme GO_id
            for j, el in subgroup.iterrows():
                list_distance = []
                list_df_index = []
                if max_index > j:
                    # Pour chaque resultat, on le compare au autre si pas deja fait
                    # le dernier a deja ete comparer au autre
                    k = j
                    max_index = i+len(subgroup)-1
                    for k in range(j, max_index):
                        k += 1
                        # j / k correspond a l'index dans la dataframe du resultat en train d'etre compare
                        field1 = df.loc[j, column_to_compare]
                        field2 = df.loc[k, column_to_compare]
                        # il faut conserver les index pour retrouver les resultat plus tard
                        list_df_index.append([j, k])
                        ratio = fuzz.token_set_ratio(field1, field2)
                        print(str(j) + '-' + str(k))
                        print(ratio)
                        # list_distance.append('toto')
                        list_distance.append(ratio)
                        k += 1
                    # on conserve le meilleur ratio entre les IP_name
                    index_max_ratio = list_distance.index(max(list_distance))
                    # on stocke pour pouvoir ensuite choisir le meilleur pour le GO_id
                    dict_tmp_annot[j] = [
                        list_df_index[index_max_ratio], list_distance[index_max_ratio]]
            # trouver la meilleur assocation de resultat pour un gene
            max_val = 0
            # for key, value in d.iteritems(): # python 2
            for key, value in dict_tmp_annot.items():  # python 3
                new_val = int(value[1])
                if new_val > max_val:
                    max_val = new_val
                    max_val_key = key
            # choisir laquel des deux annotation on prend
            # on prend le hit parmis les deux restant qui a les meilleure metrique
            # si ex-eaquo, le premier au hasard ?
            if max_val > 70:
                index_hit_1 = dict_tmp_annot[max_val_key][0][0]
                index_hit_2 = dict_tmp_annot[max_val_key][0][1]
                metric_1 = df.loc[index_hit_1, 'lrap']
                metric_2 = df.loc[index_hit_2, 'lrap']
                if float(metric_1) >= float(metric_2):
                    dict_annot_cpd[current_GO_id] = [index_hit_1, max_val]
                elif float(metric_1) < float(metric_2):
                    dict_annot_cpd[current_GO_id] = [index_hit_2, max_val]
            else:
                dict_annot_cpd.pop(current_GO_id, None)
    return dict_annot_cpd


def make_annotation(df_filter, dict_annot_cpd, list_column_to_report):
    """
    """
    list_annotation_final = []
    # for key, value in dict_annot_cpd.iteritems(): # python 2
    for key, value in dict_annot_cpd.items():  # python 3
        list_annotation_tmp = []
        list_annotation_tmp.append(key)
        for column_to_use in list_column_to_report:
            list_annotation_tmp.append(
                df_filter.loc[int(value[0])][column_to_use])
        list_annotation_tmp.append(str(value[1]))
        list_annotation_final.append(list_annotation_tmp)

    return list_annotation_final


def write_db_file(list_annotation_final, output_name):
    """
    """

    with open(output_name, 'w') as f:
        for el in list_annotation_final:
            f.write('\t'.join(el) + '\n')


if __name__ == "__main__":

    annotation_launch()
