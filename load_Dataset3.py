def dataset3():
    lncrna_name = []
    disease_name = []
    mirna_name = []

    input_lncrna_disease = [[], []]
    input_lncrna_mirna = [[], []]
    input_mirna_disease = [[], []]

    f = open("./Data/Dataset3/lncRNA_240.txt", 'r', encoding='utf-8')
    contents = f.readlines()
    for content in contents: lncrna_name.append(content.lower().strip('\n'))
    f.close()

    f = open("./Data/Dataset3/disease_412.txt", 'r', encoding='utf-16')
    contents = f.readlines()
    for content in contents: disease_name.append(content.lower().strip('\n'))
    f.close()

    f = open("./Data/Dataset3/miRNA_495.csv", 'r', encoding='utf-8')
    contents = f.readlines()
    for content in contents: mirna_name = content.lower().strip('\n').split(',')
    f.close()

    f = open("./Data/Dataset3/lncRNA_Disease.txt", 'r+', encoding='utf-8')
    contents = f.readlines()
    for content_i in range(len(contents)):
        value = contents[content_i].split('\t')
        for i in range(len(value)):
            if int(value[i].strip('\n')) == 1:
                input_lncrna_disease[0].append(int(content_i))
                input_lncrna_disease[1].append(int(i))
    f.close()

    f = open("./Data/Dataset3/lncRNA_miRNA.csv", 'r', encoding='utf-8')
    contents = f.readlines()
    for i in range(len(contents)):
        value = contents[i].lower().strip('\n').split(',')
        for j in range(1, len(value)):
            if value[j] == '1':
                input_lncrna_mirna[0].append(i)
                input_lncrna_mirna[1].append(j - 1)
    f.close()

    f = open("./Data/Dataset3/miRNA_Disease.csv", 'r', encoding='utf-8')
    contents = f.readlines()
    for i in range(len(contents)):
        value = contents[i].lower().strip('\n').split(',')
        for j in range(1, len(value)):
            if value[j] == '1':
                input_mirna_disease[0].append(i)
                input_mirna_disease[1].append(j - 1)
    f.close()

    lncrna_num = len(lncrna_name)
    disease_num = len(disease_name)
    mirna_num = len(mirna_name)

    print("\nDataset3 loading completed!")

    return lncrna_num, mirna_num, disease_num, input_lncrna_disease, input_mirna_disease, input_lncrna_mirna
