import xlrd


def dataset1():
    lncrna_name = []
    disease_name = []
    mirna_name = []

    lncrna_disease = []
    mirna_disease = []
    lncrna_mirna = []

    f = open('./Data/Dataset1/lncRNA_Disease.txt', 'r', encoding='utf-8')
    contents = f.readlines()
    for content in contents:
        value = content.split('\t')
        value[0] = value[0].lower()
        if value[0] not in lncrna_name: lncrna_name.append(value[0])
        value[1] = value[1].lower().strip('\n')
        if value[1] not in disease_name: disease_name.append(value[1])
        lncrna_disease.append(value)
    f.close()

    f = open("./Data/Dataset1/miRNA_Disease.txt", 'r', encoding='utf-8')
    contents = f.readlines()
    for content in contents:
        value = content.split(',')
        value[0] = value[0].lower()
        if value[0] not in mirna_name: mirna_name.append(value[0])
        value[1] = value[1].lower().strip('\n')
        if value[1] not in disease_name: disease_name.append(value[1])
        mirna_disease.append(value)
    f.close()

    f = "./Data/Dataset1/starBaseV2.0.xlsx"
    book_l_m = xlrd.open_workbook(f)
    sheet_l_m = book_l_m.sheet_by_name('Sheet1')
    for r in range(sheet_l_m.nrows):
        col = []
        for c in range(sheet_l_m.ncols):
            col.insert(0, sheet_l_m.cell(r, c).value.lower())
        if col[0] in lncrna_name:
            lncrna_mirna.append(col)
            if col[1] not in mirna_name: mirna_name.append(col[1])

    lncrna_num = len(lncrna_name)
    disease_num = len(disease_name)
    mirna_num = len(mirna_name)

    lncrna_index = dict(zip(lncrna_name, range(0, lncrna_num)))
    mirna_index = dict(zip(mirna_name, range(0, mirna_num)))
    disease_index = dict(zip(disease_name, range(0, disease_num)))

    input_lncrna_disease = [[], []]
    for i in range(len(lncrna_disease)):
        input_lncrna_disease[0].append(lncrna_index.get(lncrna_disease[i][0]))
        input_lncrna_disease[1].append(disease_index.get(lncrna_disease[i][1]))

    input_lncrna_mirna = [[], []]
    for i in range(len(lncrna_mirna)):
        input_lncrna_mirna[0].append(lncrna_index.get(lncrna_mirna[i][0]))
        input_lncrna_mirna[1].append(mirna_index.get(lncrna_mirna[i][1]))

    input_mirna_disease = [[], []]
    for i in range(len(mirna_disease)):
        input_mirna_disease[0].append(mirna_index.get(mirna_disease[i][0]))
        input_mirna_disease[1].append(disease_index.get(mirna_disease[i][1]))

    print("\nDataset1 loading completed!")

    return lncrna_num, mirna_num, disease_num, input_lncrna_disease, input_mirna_disease, input_lncrna_mirna