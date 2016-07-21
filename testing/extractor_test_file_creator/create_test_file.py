# coding: utf-8

import os
import sys
import csv

import Oligotyping.lib.fastalib as u

def main():
    filelist = os.listdir(os.getcwd())
    filelist.remove("create_test_file.py")
    fieldnames = ["ID", "Anticodon"]

    with open("../extractor_test_file", "w") as writefile:
        writer = csv.DictWriter(writefile, fieldnames=fieldnames,
            delimiter="\t")
        writer.writeheader()
    
        for file in filelist:
            cur_fasta = u.SequenceSource(file)
            while cur_fasta.next():
                cur_dict = {}
                cur_list = cur_fasta.id.split(" ")
                cur_dict["ID"] = cur_list[0]
                cur_dict["Anticodon"] = cur_list[5].strip("(").strip(")")
                writer.writerow(cur_dict)

if __name__ == "__main__":
    main()
