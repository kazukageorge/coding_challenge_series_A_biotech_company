'''

Description of the code:
    1. Given a gbff file, compute the most common coding sequence (CDS) features
        !python3 solution.py -g <directory_to_gbff_file> -c
    2. Given a gbff file and feature(s), output a csv with rows=CDS and column=features. There are default features.
        !python3 solution.py -g <directory_to_gbff_file> -a default
          or
        !python3 solution.py -g <directory_to_gbff_file> -a <feature1> -a <feature2> ... -a <featureN>
          or
        !python3 solution.py -g <directory_to_gbff_file> -a default -a <feature1> -a <feature2> ... -a <featureN>
    3. Given a gbff file and keyword, output csv that has a non-case sensitive match with the keyword (in the given
        assignment, it was cas9) in each feature. keyword must be cas\d, ex. "cas9"
        !python3 solution.py -g <directory_to_gbff_file> -k <keyword>
    4. Given a gbff file and keyword, output a textfile sending each CDS's translation (Amino acid) to NCBI's BLAST
        database. If the keyword is "cas", it will only look into the features that has repeated_region
        !python3 solution.py -g <directory_to_gbff_file> -k <keyword> -b

'''

import sys, getopt, os
import pandas as pd
import re
from datetime import datetime
from Bio.Blast import NCBIWWW
from Bio import SeqIO

def get_args(argv=None):
    '''
    Description:
        Check and organize the input arguments.

    Input:
        command line input arguments -> argv

    Output:
        variables       -> verbose, check_common_feat, gbff, add_feat, filt, blast

    Format:
        -h --help:                      optional, help command
        -b --blast:                     optional, run ncbi's BLAST database. Used to check the filter for discovering new filter (Cas9) sequences
                                                  takes foerever to run
        -g --gbff <directory>:          required, enter gbff file directory for analysis
        -c --check_common_features:     optional, used to check the most common features of the CDS.
        -a --add_features               optional, used to add features to add to the output csv. If entering multiple
                                                  features, please use the format "-a feature1 -a feature2 ...".
        -k --keyword:                   optional, used for determining a set of criteria

    :param argv:
    :return:
    '''
    if argv == None:
        print('solution.py --help  --blast --gbff <directory/file.gbff> --check_common_features -a feature1 -a feature2 --keyword <cas\d>')
        print('I used pyCharm with the arguments\n--verbose -c -b -g /Users/george/Desktop/GCF_002014815.1_ASM201481v1_genomic.gbff -a note -a gene -f cas9')
        print('')
        sys.exit(2)

    check_common_feat = False
    gbff = None
    keyword = None
    blast = False

    try:
        opts, args = getopt.getopt(argv, "hvcbg:a:k:",
                                   [
                                       "help",
                                       "verbose",
                                       "check_common_features",
                                       "blast",
                                       "gbff=",
                                       'add_features='
                                       "keyword="
                                   ])
        add_feat =[]
    except getopt.GetoptError as err:
        print('Invalid input arguments..')
        print('Please type \npython3 solution.py --help')
        print(str(err))
        sys.exit(2)

    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print('\
        -h --help:                      optional, help command  \n\
        -b --blast:                     optional, run NCBI\'s BLAST database. Used to check the filter for discovering new keyword (cas) candidate sequences\n\
        -g --gbff <directory>:          required, enter gbff file directory for analysis\n\
        -c --check_common_features:     optional, used to check the most common features of the CDS.\n\
        -a --add_features <feature1>    optional, used to add features to add to the output csv. If entering multiple features, please use the format "-a feature1 -a feature2 ...".\n\
        -k --keyword:                   optional, used for determining a set of criteria')

            print('\n\
        Description of the code:\n\
             1. Given a gbff file, compute the most common coding sequence (CDS) features\n\
                    !python3 solution.py -g <directory_to_gbff_file> -c\n\
             2. Given a gbff file and feature(s), output a csv with rows=CDS and column=features. There are default features.\n\
                    !python3 solution.py -g <directory_to_gbff_file> -a default\n\
                       or\n\
                    !python3 solution.py -g <directory_to_gbff_file> -a <feature1> -a <feature2> ... -a <featureN>\n\
                       or\n\
                    !python3 solution.py -g <directory_to_gbff_file> -a default -a <feature1> -a <feature2> ... -a <featureN>\n\
             3. Given a gbff file and keyword, output csv that has a non-case sensitive match with the keyword (in the given assignment, it was cas9) in each feature. keyword must be cas\d, ex. "cas9"\n\
                    !python3 solution.py -g <directory_to_gbff_file> -k <keyword>\n\
             4. Given a gbff file and keyword, output a textfile sending each CDS\'s translation (Amino acid) to NCBI\'s BLAST database. If the keyword is "cas", it will only look into the features that has repeated_region\n\
                    !python3 solution.py -g <directory_to_gbff_file> -k <keyword> -b\n\n')


            print('python3 solution.py --help --blast --gbff <directory/file.gbff> --check_common_features -a feature1 -a feature2 --keyword <cas\d>')
            print('Example: python3 solution.py -c -b -g /Users/george/Desktop/GCF_002014815.1_ASM201481v1_genomic.gbff -a note -a gene -k cas9')
            sys.exit()
        elif opt in ("-c", "--check_common_features"):
            check_common_feat = True
        elif opt in ("-b", "--blast"):
            blast = True
        elif opt in ("-g", "--gbff"):
            if "-" in arg:
                print('Invalid arguments. Cannot have "-gbff {}" '.format(arg))
                sys.exit(2)
            gbff = arg
        elif opt in ("-a", "--add_features"):
            if "-" in arg:
                print('Invalid arguments. Cannot have "-a {}" '.format(arg))
                sys.exit(2)
            add_feat.append(arg)
        elif opt in ("-k", "--keyword"):
            if "-" in arg:
                print('Invalid arguments. Cannot have "-k {}" '.format(arg))
                sys.exit(2)
            keyword = arg
            flag = re.match("cas\d", keyword.lower())
            if flag is None:
                print("-" * 100)
                print("keyword must be cas[1,2,3....], ex. 'cas9'")
                print("You entered '{}'".format(keyword))
                print("Exiting...")
                sys.exit(2)
    if not add_feat:
        add_feat = None
    if gbff is None:
        print('Invalid input arguments..')
        print('Please type \npython3 solution.py --help')
        sys.exit(2)

    print('-'*100)
    print('Input Arguments: \n')
    print('GBFF:                  {}'.format(gbff))
    print('CHECK_COMMON_FEATURES: {}'.format(check_common_feat))
    print('ADDITIONAL FEATURE(S): {}'.format(add_feat))
    print('KEYWORD:               {}'.format(keyword))
    print('RUN BLAST              {}'.format(blast))
    print('-'*100)



    return check_common_feat, gbff, add_feat, keyword, blast


def check_common_features(gbff):
    '''
    Display all features in the gbff file, add each feature into a HashMap (key=feature, value=frequency) and output all
    features in a list based on its frequency

    :param gbff:
    :return:
    '''
    print('Checking common features in the CDS of {}: \n'.format(gbff))
    hashMap = {}
    records = SeqIO.parse(gbff, format='genbank')
    # iterate over the records, find the features in feature.qualifiers dictionary
    for record in records:
        for feature in record.features:
            if feature.type == 'CDS':
                for key in feature.qualifiers:
                    if key not in hashMap:
                        hashMap[key] = 1
                    else:
                        hashMap[key] += 1

    # # sort the hashMap
    # {k: v for k, v in sorted(hashMap.items(), key=lambda item: item[1])}
    sort_features = sorted(hashMap.items(), key=lambda x: x[1], reverse=True)

    for i in sort_features:
        print(i[0], i[1])

    print('-'*100)

def get_nearest_region(record_locs,list_repeat_region, record_id ):
    '''
    Description:
        Computes the nearest_repeat_region_nucleotide. nearest_repeat_region_nucleotide is defined as the shortest distance
        from any part of the subject sequence to any part of the target sequence

    Input:
        record_locs             -> List
        list_repeat_region      -> List
        record_id               -> str
    Output:
        repeat_dist_nucleotide  -> List

    :param record_locs:
    :param list_repeat_region:
    :param record_id:
    :return:
    '''

    if record_id not in list_repeat_region:
        return ["infinite"]*len(record_locs)
    # get the distance to the nearest repeating region (either head of string or tail of the string)
    # from the head or the tail of the subject string
    repeat_region_nuc = list_repeat_region[record_id]

    repeat_dist_nucleotide = []
    # Find the nearest location. This will be seq1-tail to head-seq2 or seq2-tail to head-seq1.
    for record_loc  in record_locs:
        if max(record_loc) < min(repeat_region_nuc):
            repeat_dist_nucleotide.append(min(repeat_region_nuc) - max(record_loc))
        else:
            repeat_dist_nucleotide.append(min(record_loc) - max(repeat_region_nuc))
    return repeat_dist_nucleotide

def get_features_to_csv(gbff, add_feat):
    '''
     Description:
        Computes the features (default and additional features) given CDS. Iterate over gbff file using Bio, genebank unwrapper
        and (1) append feature values if it is CDS section (2) separately append repeat_region

        The following are default features:
            1. organism
            2. record_id
            3. feature_location
            4. locus_tag
            5. product
            6. protein_id
            7. protein_size
            8. translation (Amino acid sequence)
            9. nearest_repeat_region_nucleotide
            - - - -
            10 onwards. added features
            - - - -

    Input:
        gbff                    -> directory
        add_feat                -> List
        verbose                 -> Bool
    Output:
        df_features             -> DataFrame

    :param gbff:
    :return: list of features
    '''
    # print selected features
    save = 1
    if not add_feat:
        add_feat= ['default']
        save = 0

    if save:
        print('Obtaining CSV for the features ', end=" "),
    else:
        print('Obtain DataFrame only')
    default_marker = 0
    if add_feat:
        for i, feat in enumerate(add_feat):
            if feat == "default":
                feat = "default_values"
                default_marker =1
            print('{}.{} '.format(i+1,feat), end=" "),
        print('\n')

        # remove default from the features.
        if "default" in add_feat:
            add_feat.remove("default")

    default_vals = ["organism", "record_id", "feature_location","locus_tag","protein_id", "gene", "protein_size",
                    "distance_to_repeat_region", "translation"]
    if default_marker:
        print("default_values are: ", end="")
        for i, default_val in enumerate(default_vals):
            print('{}.{} '.format(i + 1, default_val), end=" "),
        print("\n")
    if add_feat:
        # remove duplicate features
        print('Removing duplicated features..')
        add_feat = list(set(add_feat))
        for af in add_feat:
            if af in default_vals:
                add_feat.remove(af)

    # key = record_id
    # value = [start, end]
    list_repeat_region = {}

    # [start, end]
    list_nearest_region = []

    # print('Obtaining features given cds for {}'.format(gbff))
    find_features = []
    if default_marker:
        find_features = ['locus_tag', 'product','protein_id','gene','protein_size','translation']
        find_features += add_feat
    records = SeqIO.parse(gbff, format='genbank')
    # iterate over the records, find the features in feature.qualifiers dictionary

    # column used for data frame. More will be appended later
    list_feature_header = []
    if default_marker:
        list_feature_header = ['organism', 'record_id','feature_location']
        list_feature_header += find_features
    list_features = []
    print('Organizing to DataFrame and searching for repeat region (for CRISPR identifier..)')
    count_rec = 0
    count_feat = 0
    # Iterate over records, append the features to the approriate index on the list
    for record in records:
        # if count_rec %10 == 0:
        #     print("\nAppending Record #{}".format(count_rec),end="")
        temp_list_record = []
        record_organism = record.annotations['organism']
        record_id = record.id

        temp_list_record.append(record_organism)
        temp_list_record.append(record_id)
        record_loc = []

        for feature in record.features:
            temp_list_features = []
            # Only keep log for 'CDS' and 'repeat_region'
            if feature.type == 'CDS':
                feature_location = feature.location
                temp_list_features.append(str(feature_location))
                for find_feature in find_features:
                    if find_feature == 'protein_size':
                        continue
                    if find_feature in feature.qualifiers:
                        if find_feature == 'translation':
                            n = str(len(feature.qualifiers[find_feature][0]))
                            temp_list_features.append(n)
                        temp_list_features.append(feature.qualifiers[find_feature][0])
                    else:
                        if find_feature == 'translation':
                            temp_list_features.append(' ')
                        temp_list_features.append(' ')
                list_features.append((temp_list_record+temp_list_features))
                record_loc.append((int(feature.location.start), int(feature.location.end)))

            elif feature.type == 'repeat_region':
                print("")
                print('Found repeat region! Used for identifying Cas candidate sequences')
                print("")

                if list_repeat_region.get(record.id) is None:
                    list_repeat_region[record.id] = [
                        int(feature.location.start), int(feature.location.end)
                    ]
                else:
                    list_repeat_region[record.id].append(feature.location.start, feature.location.end)
            # if count_feat % 1000 == 0:
            #     print("Appending feature #{}".format(count_feat),end = "")
            count_feat +=1
        count_rec +=1
        list_nearest_region+=get_nearest_region(record_loc, list_repeat_region, record_id)


    print("Appeneded a total of {} records, {} features".format(count_rec, count_feat))
    # convert the list to dataframe. Easier to append two DF than two lists
    df_features = pd.DataFrame(list_features, columns=list_feature_header)
    df_nearest_region = pd.DataFrame(list_nearest_region, columns = ['nearest_repeat_region_nucleotide'])
    df_features = pd.concat([df_features, df_nearest_region], axis=1)
    cols = ['organism', 'locus_tag', 'protein_id', 'locus_tag', 'record_id',
                                'feature_location', 'protein_size' , 'gene',  'nearest_repeat_region_nucleotide']
    cols = cols + add_feat + ['translation']

    # Organize the dataframe
    df_features = df_features[cols]
    if save:
        name_feat = ''
        if default_marker:
            name_feats = ['default'] + add_feat
            for name in name_feats:
                name_feat += "_" + name

        else:
            name_feats = add_feat
            for name in name_feats:
                name_feat +="_" + name

        now = datetime.now()
        current_time = now.strftime("%d%m%Y_%H%M%S")

        pattern =  'GCF' + '(?P<name>.*)\.gbff$'
        # pattern =  os.path.sep + '(?P<name>.*)\.gbff$'

        name_re = re.findall(pattern , gbff)
        name_gbff = 'GCF' +name_re[0] + '_gbff'
        save_name = ('output/CDS_features/CDS_features{}_{}_{}.csv'.format(name_feat,name_gbff, current_time))
        # name =
        # match = pattern.search(gbff)    # gbff_name = re.match(pattern, gbff)
        print('Saving the DataFrame to "{}"'.format(save_name))
        # create output folder
        if not os.path.exists('output/CDS_features'):
            os.makedirs('output/CDS_features')
            print('output folder created in {}'.format(os.getcwd()))
        df_features.to_csv(save_name)

    return df_features

def filter_given_protein_features(features_df, filt):
    '''
    Description:
        Used the filter keyword to find the target of interest. If found, the row that the target resides in the dataframe
        will be written into a csv

        Iterate through the row of the dataframe, only looking at the features "gene" and "product"
        If there is a match (regular expression) then find the row index, make csv its row with other features

    Input:
        features_df             -> dataframe, features and its info
        filt                    -> string, filter keyword
    Output:
        None

    :param features_df:
    :param filters:
    :return:
    '''
    print("-"*100)
    print('Searching for keyword {} in DataFrame..\n'.format(keyword))

    idx_gene = idx_product = -1
    row_gene = row_product = []
    if "gene" in features_df:
        # print('Using feature "gene"')
        for gene in features_df.gene.tolist():
            if gene == filt:
                row_gene = features_df.loc[features_df['gene'] == filt]
                idx_gene = row_gene.index[0]
                print("{} is present at row index {}. Protein ID is {} \n".format(filt, idx_gene,row_gene['protein_id'].tolist()[0]))
    if "product" in features_df:
        # print('Using feature "product"')
        #product is a default feature
        for sentence in features_df['product'].tolist():
            if re.findall(filt.lower(), sentence.lower()):
                row_product = features_df.loc[features_df['product'] == sentence]
                idx_product = row_product.index[0]
                print("{} is present at row index {}. Protein ID is {} \n".format(filt, idx_product,row_product['protein_id'].tolist()[0]))

    # If there is a hit
    if idx_gene > -1 or idx_product > -1:
        print("Saving DataFrame with the filter/keyword: {}".format(keyword))
        # create output folder
        if not os.path.exists('output/filter_{}'.format(keyword)):
            os.makedirs('output/filter_{}'.format(keyword))
            print('Creating and saving csv to "output/filter_{}"'.format(keyword))
    now = datetime.now()
    current_time = now.strftime("%d%m%Y_%H%M%S")
    pattern = 'GCF' + '(?P<name>.*)\.gbff$'
    # pattern =  os.path.sep + '(?P<name>.*)\.gbff$'

    name_re = re.findall(pattern, gbff)
    name_gbff = 'GCF' + name_re[0] + '_gbff'
    save_name = ('output/filter_{}/{}_{}.csv'.format(filt, keyword, current_time))

    # save appropriate dataframe column
    if idx_gene ==-1 and idx_product ==-1:
        print('Filter "{}" not found'.format(filt))
    if idx_gene > -1 and idx_product == -1:
        print('Filter "{}" found from feature "gene"'.format(filt))
        print('Saving filtered DataFrame with keyword {} to "{}"'.format(keyword, save_name))

        row_gene.to_csv('output/filter_{}/{}_{}.csv'.format(filt, keyword, current_time))
    if idx_gene == -1 and idx_product > -1:
        print('Filter "{}" found from feature "product"'.format(filt))
        print('Saving filtered DataFrame with keyword {} to "{}"'.format(keyword, save_name))

        row_product.to_csv('output/filter_{}/{}_{}.csv'.format(filt, keyword, current_time))
    if idx_gene > -1 and idx_product > -1:
        print('Filter "{}" found from features "gene" and "product"'.format(filt))
        print('Saving filtered DataFrame with keyword {} to "{}"'.format(keyword, save_name))

        if idx_gene == idx_product:
            row_gene.to_csv('output/filter_{}/{}_{}.csv'.format(filt,keyword,current_time))
        else:
            # combine both dataframe
            row_gene = pd.concat([row_gene, row_product])
            row_gene.to_csv('output/filter_{}/{}_{}.csv'.format(filt,keyword, current_time))

    pass

def run_blast_internet(features_df, filt):
    '''
     Description:
        Used the filter keyword to find the target of interest. If found, the row that the target resides in the dataframe
        will be written into a textfile.

        Use NCBI's BLAST to demystify if the subject of interest is cas9 based on E-value and percent match
        Running BLAST over the internet takes forever. Better to run it locally, but
        for the purpose of this problem, we keep everything simple and just run it over the net.
        Much faster to just copy and paste the Amino acid sequence to the database manually
        <https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome>

    Input:
        features_df             -> dataframe, features and its info
        filt                    -> string, filter keyword
    Output:
        None


    :param features_df:
    :param filt:
    :return:
    '''
    print("-"*100)
    print('Searching for potential {} candidate sequences via BLAST\n'.format(filt))
    print('Only consider CDS near repeated region if any..')

    # make the dataframe into variables
    nearest_repeat_region_nuc = features_df['nearest_repeat_region_nucleotide']
    near_rep_nucs = features_df.nearest_repeat_region_nucleotide.tolist()
    amino_acids = features_df.translation.tolist()
    # locus_tags = features_df.locus_tag.tolist()
    protein_ids = features_df.protein_id.tolist()

    # iterate over the variables, if there is a repeated region, then send the sequence to the database
    # use filter=command_line_argument_Filter (cas9) to just filter out useful info
    print('Sending the Amino acid sequence to BLAST and retrieving the results...')
    for near_rep_nuc, amino_acid, protein_id in zip(near_rep_nucs, amino_acids, protein_ids):
        if near_rep_nuc is not "infinite" and len(amino_acid) > 1:
            if amino_acid:
                print('Near repeated region sequence found!')
                # Put in fasta format for searching
                s = '>{}.\n{}'.format(protein_id, amino_acid)
                print('Blasting protein_id: "{}", Amino acid: "{}"\nThis may take a while'.format(protein_id, amino_acid))
                result_handle = NCBIWWW.qblast(program="blastp", database='nr', sequence=s, format_type="Text",
                                               filter=filt)

                print("Saving the result from BLAST for protein ID: {} with the filter/keyword: {}".format(protein_id,filt))
                # create output folder
                if not os.path.exists('output/BLAST/'):
                    os.makedirs('output/BLAST/')
                    print('Creating and saving csv to "output/BLAST/raw/{}_{}.txt"'.format(protein_id, filt))

                save_file = open('output/BLAST/{}_{}.txt'.format(protein_id, filt), 'w')
                blast_results = result_handle.read()
                save_file.write(blast_results)
                save_file.close()

    print('BLAST complete')
    print('')

    pass

if __name__ == '__main__':
    # Check if input arguments are valid, and parse them
    check_feat, gbff, add_feat, keyword, blast = get_args(sys.argv[1:])


    # check the most common features (features.qualifiers)
    if check_feat:
        check_common_features(gbff)

    # outputs a dataframe. Also generates a csv file
    if add_feat or keyword:
        features_df = get_features_to_csv(gbff, add_feat)

    if keyword:
        filter_given_protein_features(features_df, keyword)
        if blast:
            run_blast_internet(features_df, keyword)

    print('solution.py complete!')
