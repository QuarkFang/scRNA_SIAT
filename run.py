import dataset


if __name__ == '__main__':
    # dataset.add(6, 'Human', 'Kidney', 'UBERON_0002113', 'Normal', 'Normal cell',
    #             'Proximal tubular cell', 'NA', 'AA1',
    #             '9263997')

    a = dataset.query('AAA')
    print(a)