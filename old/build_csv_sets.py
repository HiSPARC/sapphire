import tables
import csv


DATAFILE = 'analysis-e14.h5'


def write_csv_set(data, tablename, columns):
    with open(tablename + '-e14.csv', 'w') as file:
        writer = csv.writer(file, delimiter='\t')

        table = getattr(data.root, tablename)
        for row in table:
            writer.writerow([row[x] for x in columns])


if __name__ == '__main__':
    data = tables.openFile(DATAFILE, 'r')

    columns = ['PID', 'R', 'PHI', 'T', 'E', 'pid', 'r', 'phi', 't', 'e']
    write_csv_set(data, 'g50_60_0', columns)
    write_csv_set(data, 'g50_60_6', columns)
    write_csv_set(data, 'g50_60_10', columns)
    write_csv_set(data, 'g20_25_0', columns)
    write_csv_set(data, 'g20_25_6', columns)
    write_csv_set(data, 'g20_25_10', columns)
    write_csv_set(data, 'g10_15_0', columns)
    write_csv_set(data, 'g10_15_6', columns)
    write_csv_set(data, 'g10_15_10', columns)
