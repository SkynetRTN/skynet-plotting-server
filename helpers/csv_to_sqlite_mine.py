import os

import csv_to_sqlite
import csv
from cluster_pro_scraper import hms2d, dms2d

raw_file_path = "clusters.csv"
new_file_path = "MWSC.csv"
sqlite_path = "../MWSC.sqlite"

# os.remove(new_file_path)

with open(raw_file_path, 'r', encoding='utf-8') as raw_file:
    reader = csv.DictReader(raw_file, delimiter=',', quotechar='\"')
    writefile = open(new_file_path, 'a', newline='', encoding='utf-8')
    writer = csv.writer(writefile, delimiter=',', quotechar='\'')
    header = False
    for row in reader:
        if not header:
            writer.writerow(row.keys())
            header = True
        new_row = row
        new_row['ra'] = str(hms2d(row['ra']))
        new_row['dec'] = str(dms2d(row['dec']))
        print(list(new_row.values()))
        writer.writerow(list(new_row.values()))


writefile.close()

options = csv_to_sqlite.CsvOptions(encoding="utf-8")
try:
    os.remove(sqlite_path)
finally:
    csv_to_sqlite.write_csv([new_file_path], sqlite_path, options)
