# /usr/bin/python3

# comparing annotation method results
# 1 removing stop word with nltk
# 2 for each gene comparing each annotation result and if ratio > treeshold keep one annotation

from nltk.corpus import stopwords 
from nltk.tokenize import word_tokenize 
from fuzzywuzzy import fuzz
import pandas as pd
from nltk.tokenize import RegexpTokenizer
import operator

def example1():
    """
    """
    example_sent = "This is a sample sentence, showing off the stop words filtration."

    stop_words = set(stopwords.words('english')) 

    word_tokens = word_tokenize(example_sent)

    filtered_sentence = [w for w in word_tokens if not w in stop_words] 

    filtered_sentence = [] 
    for w in word_tokens: 
        if w not in stop_words: 
            filtered_sentence.append(w) 

    print(word_tokens) 
    print(filtered_sentence) 


def example2():
    """
    """
    a=fuzz.token_set_ratio('Deluxe Room, 1 King Bed', 'Deluxe King Room')
    print(a)
#def get_ratio(row):
#    name = row['Expedia']
#    name1 = row['Booking.com']
#    return fuzz.token_set_ratio(name, name1)
#len(df[df.apply(get_ratio, axis=1) > 70]) / len(df)

def open_file(sid):
    """
    """   
    filepath = str(sid)+'/IPRSCAN'+str(sid)
    filepath = str(sid)+'/testIPR'
    file = open(filepath, 'r')
    df_iprscan = pd.read_csv(file, sep='\t', header=0, dtype='str')
    file.close()

    filepath = str(sid)+'/UNIFIRE'+str(sid)
    file = open(filepath, 'r')
    df_unifire = pd.read_csv(file, sep='\t', header=0, dtype='str')
    file.close()

    # filepath = str(sid)+'/Swissprot'+str(sid)
    # file = open(filepath, 'r')
    # df_swissprot = pd.read_csv(file, sep='\t', header=0, dtype='str')
    # file.close()

    # filepath = str(sid)+'/TrEMBL'+str(sid)
    # file = open(filepath, 'r')
    # df_trembl = pd.read_csv(file, sep='\t', header=0, dtype='str')
    # file.close()

    # filepath = str(sid)+'/annot'+str(sid)
    # file = open(filepath, 'r')
    # df_annot = pd.read_csv(file, sep='\t', header=0, dtype='str')
    # file.close()
    df_swissprot=""
    df_trembl=""
    df_annot=""


    return df_iprscan, df_unifire, df_swissprot, df_trembl, df_annot


def clean_data(df, list_column_to_keep, column_to_clean, list_column_to_filter, list_value_of_filter):
    """
    Return a DataFrame with only the column to keep that pass the filters values,
    the punctuation and stop world will be removed fromthe column to clean
    :param df: DataFrame to clean
    :type df: DataFrame
    :param list_column_to_keep: list of column name to keep
    :type list_column_to_keep: list ['str']
    :param column_to_clean: column name to clean
    :type column_to_clean: str
    :param list_column_to_filter: list of column name to filter
    :type list_column_to_filter: list ['str']
    :param list_value_of_filter: list of filter values treeshold data >= values
    :type list_value_of_filter: list ['float']
    """
    
    # enleve les NULL
    df = df.dropna(subset = list_column_to_keep)
    # remet l'index des row 
    df.reset_index(inplace=True)

    #cree le dataframe de travail juste avec les colonnes voulues
    new_df = pd.DataFrame()
    df_filter = pd.DataFrame()
    df_clean_filter = pd.DataFrame()


    for i in range(0, len(list_column_to_keep)):
        col = list_column_to_keep[i]
        new_df.insert(i, col, df[col], True)
        df_clean_filter[col]=""
        df_filter[col]=""

    # liste des stop world = determinant, mot de liaison ...
    stop_words = set(stopwords.words('english')) 
    # pour la ponctuation
    tokenizer = RegexpTokenizer(r'\w+')

    count=0
    # itere sur dataframe
    for i, row in new_df.iterrows():
        new_desc = ''
        result = []
        is_valid = True
        # filtre la qualite de resultat a garder
        if list_column_to_filter != []:
            for j in range(0, len(list_column_to_filter)):
                filter_name = list_column_to_filter[j]
                filter_value = float(list_value_of_filter[j])
                if not float(row[filter_name]) >= filter_value:
                    #new_df = new_df.drop(i)
                    is_valid = False
        if is_valid:
            # remplis le df filtrer avec column_to_clean d'origine
            df_filter = df_filter.append(row)
            count+=1
            # selection de tout les mots + ponctuation de column_to_clean
            word_tokens = word_tokenize(row[column_to_clean])
            for word in word_tokens:
                # si le mots de column_to_clean n'est pas un stop word, on le conserve
                if not word in stop_words:
                    new_desc += word + ' '
            # enleve ponctuation
            result = tokenizer.tokenize(new_desc)
            # remet GON_description nettoye
            row[column_to_clean] = ''
            for element in result:
                row[column_to_clean] += element + ' '
            # remplis le df clean_filtrer avec column_to_clean nettoye
            df_clean_filter = df_clean_filter.append(row)

    df_filter.reset_index(inplace=True, drop=True)
    df_clean_filter.reset_index(inplace=True, drop=True)
    return df_clean_filter, df_filter

def compare_annotation(df):
    """
    """

    dict_annot_cpd={}
    for i, row in df.iterrows():
        current_GO_id = row['GO_id']
        # On fait des sous groupe de resultat pour le meme GO_id
        subgroup = df.loc[df['GO_id'] == current_GO_id]
        # si le Go_id a plusieurs resultat est n'a pas deja ete traite
        if len(subgroup) > 1 and current_GO_id not in dict_annot_cpd.keys():
            dict_annot_cpd[current_GO_id] = ""
            dict_tmp_annot={}
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
                        k+=1
                        # j / k correspond a l'index dans la dataframe du resultat en train d'etre compare
                        field1 = df.loc[j, 'IP_name']
                        field2 = df.loc[k, 'IP_name']
                        # il faut conserver les index pour retrouver les resultat plus tard
                        list_df_index.append([j, k])
                        ratio = fuzz.token_set_ratio(field1, field2)
                        print(str(j) + '-' + str(k))
                        print(ratio)
                        #list_distance.append('toto')
                        list_distance.append(ratio)
                        k+=1
                    # on conserve le meilleur ratio entre les IP_name
                    index_max_ratio = list_distance.index(max(list_distance))
                    # on stocke pour pouvoir ensuite choisir le meilleur pour le GO_id
                    dict_tmp_annot[j] = [list_df_index[index_max_ratio], list_distance[index_max_ratio]]
            # trouver la meilleur assocation de resultat pour un gene
            max_val=0
            #for key, value in d.iteritems(): # python 2
            for key, value in dict_tmp_annot.items(): # python 3
                new_val = int(value[1])
                if  new_val > max_val:
                    max_val=new_val
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
    #for key, value in dict_annot_cpd.iteritems(): # python 2
    for key, value in dict_annot_cpd.items(): # python 3
        list_annotation_tmp = []
        list_annotation_tmp.append(key)
        for column_to_use in list_column_to_report:
            list_annotation_tmp.append(df_filter.loc[int(value[0])][column_to_use])
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
    # sid = 10057
    sid = 1389

    df_iprscan, df_unifire, df_swissprot, df_trembl, df_annot = open_file(sid)

    ## IPRSCAN 'GO_id', 'lrap', 'IP_id', 'GON_id', 'IP_name', 'GON_description'
    list_column_to_keep = ['GO_id', 'lrap', 'IP_id', 'CPD_id', 'IP_name']
    column_to_clean = 'IP_name'
    list_column_to_filter = ['lrap']
    list_value_of_filter = [0.8]
    df_clean_filter_iprscan, df_filter_iprscan = clean_data(df_iprscan, list_column_to_keep, column_to_clean, \
    list_column_to_filter, list_value_of_filter)
    print(df_filter_iprscan)
    dict_annot_cpd = compare_annotation(df_clean_filter_iprscan)
    list_column_to_report = ['CPD_id', 'lrap', 'IP_name']
    list_annotation_final = make_annotation(df_filter_iprscan, dict_annot_cpd, list_column_to_report)
    write_db_file(list_annotation_final, 'bla_bla')

    # ## UNIFIRE GO_id, UNI_id, UNI_value
    # list_column_to_keep = ['GO_id', 'UNI_id', 'UNI_value']
    # column_to_clean = 'UNI_value'
    # list_column_to_filter = ['lrap']
    # list_value_of_filter = [0.8]
    # df_clean_filter_unifire, df_filter_unifire= clean_data(df_unifire, list_column_to_keep, column_to_clean, \
    # list_column_to_filter, list_value_of_filter)
    # print(df_clean_unifire.head(10))
    

    
    # example2()
