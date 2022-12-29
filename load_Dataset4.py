def dataset4():
    lncrna_name = []
    disease_name = []
    mirna_name = []

    dataset_lncrna_disease = []
    dataset_mirna_disease = []
    dataset_lncrna_mirna = []

    f = open("./Data/Dataset4/lncRNA_Disease.txt", 'r', encoding='utf-8')
    contents = f.readlines()
    disease_name = contents[0].split('\t')
    del disease_name[0]
    disease_name = [data.strip('\n').lower() for data in disease_name]

    for content in contents[1:]:
        value = content.split('\t')
        lncrna_name.append(value[0].strip().lower().strip('\n'))
        for i in range(len(value[1:])):
            if int(value[i + 1]) == 1:
                link = []
                link.append(value[0].strip().lower().strip('\n'))
                link.append(disease_name[i])
                dataset_lncrna_disease.append(link)
    f.close()

    f = open("./Data/Dataset4/miRNA_Disease.txt", 'r', encoding='utf-8')
    contents = f.readlines()
    for content in contents[1:]:
        value = content.split('\t')
        mirna_name.append(value[0].strip().lower().strip('\n'))
        for i in range(len(value[1:])):
            if int(value[i + 1]) == 1:
                link = []
                link.append(value[0].strip().lower().strip('\n'))
                link.append(disease_name[i])
                dataset_mirna_disease.append(link)
    f.close()

    f = open("./Data/Dataset4/miRNA_lncRNA.txt", 'r', encoding='utf-8')
    contents = f.readlines()
    for content in contents[1:]:
        value = content.split('\t')
        for i in range(len(value[1:])):
            if int(value[i + 1]) == 1:
                link = []
                link.append(lncrna_name[i])
                link.append(value[0].strip().lower().strip('\n'))
                dataset_lncrna_mirna.append(link)
    f.close()

    lncrna_num = len(lncrna_name)
    disease_num = len(disease_name)
    mirna_num = len(mirna_name)

    lncrna_index = dict(zip(lncrna_name, range(0, lncrna_num)))
    mirna_index = dict(zip(mirna_name, range(0, mirna_num)))
    disease_index = dict(zip(disease_name, range(0, disease_num)))

    input_lncrna_disease = [[], []]
    for i in range(len(dataset_lncrna_disease)):
        input_lncrna_disease[0].append(lncrna_index.get(dataset_lncrna_disease[i][0]))
        input_lncrna_disease[1].append(disease_index.get(dataset_lncrna_disease[i][1]))

    input_lncrna_mirna = [[], []]
    for i in range(len(dataset_lncrna_mirna)):
        input_lncrna_mirna[0].append(lncrna_index.get(dataset_lncrna_mirna[i][0]))
        input_lncrna_mirna[1].append(mirna_index.get(dataset_lncrna_mirna[i][1]))

    input_mirna_disease = [[], []]
    for i in range(len(dataset_mirna_disease)):
        input_mirna_disease[0].append(mirna_index.get(dataset_mirna_disease[i][0]))
        input_mirna_disease[1].append(disease_index.get(dataset_mirna_disease[i][1]))

    print("\nDataset4 loading completed!")

    return lncrna_num, mirna_num, disease_num, input_lncrna_disease, input_mirna_disease, input_lncrna_mirna
